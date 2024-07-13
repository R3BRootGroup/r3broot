/******************************************************************************
 *   Copyright (C) 2019 GSI Helmholtzzentrum für Schwerionenforschung GmbH    *
 *   Copyright (C) 2019-2023 Members of R3B Collaboration                     *
 *                                                                            *
 *             This software is distributed under the terms of the            *
 *                 GNU General Public Licence (GPL) version 3,                *
 *                    copied verbatim in the file "LICENSE".                  *
 *                                                                            *
 * In applying this license GSI does not waive the privileges and immunities  *
 * granted to it by virtue of its status as an Intergovernmental Organization *
 * or submit itself to any jurisdiction.                                      *
 ******************************************************************************/

#include "R3BNeulandMillepede.h"
#include "Utilities.h"
#include <R3BException.h>
#include <R3BNeulandCalToHitParTask.h>
#include <R3BNeulandCommon.h>
#include <SteerWriter.h>
#include <TGraphErrors.h>
#include <fmt/ranges.h>
#include <optional>
#include <range/v3/algorithm.hpp>
#include <range/v3/numeric.hpp>
#include <range/v3/view.hpp>

namespace rng = ranges;

constexpr auto SCALE_FACTOR = 10.F;
constexpr auto REFERENCE_BAR_NUM = 25;

namespace R3B::Neuland::Calibration
{
    void MillepedeEngine::Init()
    {
        cal_to_hit_par_ = GetTask()->GetCal2HitPar();
        const auto max_module_num = std::min(GetMaxModuleNum(), GetModuleSize());

        position_cal_.set_max_number_of_modules(max_module_num);

        mille_handler_.init();
        mille_handler_.set_max_number_of_modules(max_module_num);
        mille_handler_.set_scale_factor(SCALE_FACTOR);

        pede_handler_.init();
        pede_handler_.set_max_number_of_modules(max_module_num);
        pede_handler_.set_scale_factor(SCALE_FACTOR);
        pede_handler_.set_data_filename(mille_handler_.get_filename());
    }

    void MillepedeEngine::set_minimum_values(const CalData& signals)
    {
        // make sure only one hit exists in one bar
        auto t_sum = 0.;
        for (const auto& [plane_num, plane_cal_data] : signals)
        {
            for (const auto& [bar_num, bar_signal] : plane_cal_data.bar_cal_data)
            {
                if (not has_only_one_signal(bar_signal))
                {
                    continue;
                }
                const auto& left_signal = bar_signal.left.front();
                const auto& right_signal = bar_signal.right.front();
                t_sum += (left_signal.leading_time - left_signal.trigger_time + right_signal.leading_time -
                          right_signal.trigger_time)
                             .value;
            }
        }

        mille_handler_.set_avarage_t_sum(static_cast<float>(t_sum / GetBarCalDataSize(signals)));
        R3BLOG(debug, fmt::format("Average t_sum is calculated to be {}", average_t_sum_.value()));
    }

    auto MillepedeEngine::SignalFilter(const CalData& signals) -> bool
    {
        // select out rays with few hits
        if (GetBarCalDataSize(signals) < minimum_hit_)
        {
            return false;
        }

        // Make sure the track has large zenith angle
        auto is_too_vertical = rng::any_of(
            signals, [this](const auto& plane_cal) { return plane_cal.second.bar_cal_data.size() > plane_max_hit_; });
        if (is_too_vertical)
        {
            return false;
        }

        set_minimum_values(signals);
        return true;
    }

    auto MillepedeEngine::get_bar_cal_residual(const BarCalData& signal, double pred) -> double
    {
        const auto position = calculate_position_along_bar(signal, cal_to_hit_par_);
        return std::abs(pred - position);
    }

    void MillepedeEngine::AddSignals(const CalData& signals)
    {
        // all bar signal must have one signal on both sides

        switch (current_state_)
        {
            case State::pre_calibration:
                position_cal_.fill_time_differences(signals);
                break;
            case State::millepede_calibration:
                fill_data_to_mille(signals);
                break;
        }
    }

    template <typename FillAction>
    void MillepedeEngine::fill_filtered_bar_cal_data(const CalData& signals, FillAction&& fill_action)
    {
        auto max_module_num = std::min(GetMaxModuleNum(), GetModuleSize());
        for (const auto& [_, plane_cal_data] : signals)
        {
            auto plane_num = plane_cal_data.plane_num;

            // only look for the bar data with 1 signal from the left and 1 from the right
            auto filtered_signals = rng::filter_view(
                plane_cal_data.bar_cal_data | rng::views::all,
                [max_module_num](const auto& bar_signal)
                { return has_only_one_signal(bar_signal.second) and bar_signal.second.module_num <= max_module_num; });

            if (filtered_signals.empty())
            {
                return;
            }

            // check if bar signals are together
            const auto min_bar_num = rng::min(filtered_signals | rng::views::keys);
            const auto max_bar_num = rng::max(filtered_signals | rng::views::keys);
            const auto filtered_size =
                static_cast<int>(rng::distance(filtered_signals.begin(), filtered_signals.end()));
            auto is_together = (max_bar_num - min_bar_num + 1 == filtered_size);
            if (not is_together)
            {
                // fmt::print("Plane signals has seperate bar signals: {}\n", plane_cal_data);
                return;
            }

            fill_action(plane_num, filtered_signals);
        }
    }

    void MillepedeEngine::fill_data_to_mille(const CalData& signals)
    {
        auto fill_bar_position_action = [this](auto plane_num, auto& filtered_bar_signals)
        {
            const auto filtered_size =
                static_cast<double>(rng::distance(filtered_bar_signals.begin(), filtered_bar_signals.end()));
            auto displacement = std::accumulate(filtered_bar_signals.begin(),
                                                filtered_bar_signals.end(),
                                                double{},
                                                [](double init, const auto& bar_signal)
                                                { return init + GetBarVerticalDisplacement(bar_signal.first); }) /
                                filtered_size;
            mille_handler_.add_spacial_local_constraint(plane_num, displacement, BarSize_XY * filtered_size / 2.);
        };

        fill_filtered_bar_cal_data(signals, fill_bar_position_action);
        if (not mille_handler_.fit_local())
        {
            return;
        }

        auto fill_time_action = [this](auto plane_num, auto& filtered_bar_signals)
        {
            const auto filtered_size =
                static_cast<int>(rng::distance(filtered_bar_signals.begin(), filtered_bar_signals.end()));
            if (filtered_size == 1)
            {
                mille_handler_.add_signal_t(filtered_bar_signals.begin()->second, cal_to_hit_par_);
                return;
            }
            const auto is_horizontal = IsPlaneIDHorizontal(plane_num - 1);
            const auto pos_z = PlaneID2ZPos<double>(plane_num - 1);
            const auto pred = mille_handler_.get_bar_positions().get_prediction(is_horizontal, pos_z);

            // TODO: need to optimized here. get_bar_cal_residual get called too many times

            // choose the bar signal closest to the regression line
            auto bar_iter = rng::min_element(filtered_bar_signals,
                                             [this, pred](const auto& first_bar, const auto& second_bar) {
                                                 return get_bar_cal_residual(second_bar.second, pred) >
                                                        get_bar_cal_residual(first_bar.second, pred);
                                             });
            if (get_bar_cal_residual(bar_iter->second, pred) < pos_residual_threshold)
            {
                mille_handler_.add_signal_t(bar_iter->second, cal_to_hit_par_);
            }
            return;
        };

        fill_filtered_bar_cal_data(signals, fill_time_action);
    }

    void MillepedeEngine::Calibrate(Cal2HitPar& hit_par)
    {
        switch (current_state_)
        {
            case State::pre_calibration:
                pre_calibrate(hit_par);
                break;
            case State::millepede_calibration:
                millepede_calibrate(hit_par);
                break;
        }
    }

    void MillepedeEngine::pre_calibrate(Cal2HitPar& cal_to_hit_par)
    {
        R3BLOG(info, "Begin histogram calibration");
        position_cal_.calibrate(cal_to_hit_par);
        current_state_ = State::millepede_calibration;
        pede_handler_.init_steer_writer(cal_to_hit_par);
        R3BLOG(info, "Histogram calibration finished. Ready for millepede calibration!");
    }

    void MillepedeEngine::millepede_calibrate(Cal2HitPar& cal_to_hit_par)
    {
        R3BLOG(info, "Starting millepede calibration...");
        mille_handler_.close();
        pede_handler_.calibrate();
        pede_handler_.fill_module_parameters(cal_to_hit_par);
        R3BLOG(info, "Millepede calibration finished.");
    }

    void MillepedeEngine::EndOfEvent(unsigned int /*event_num*/) { mille_handler_.do_end_of_event(); }

    void MillepedeEngine::EndOfTask()
    {
        average_t_sum_.reset();
        mille_handler_.buffer_clear();
    }
} // namespace R3B::Neuland::Calibration
