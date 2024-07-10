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
#include <TFitResult.h>
#include <TGraphErrors.h>
#include <fmt/ranges.h>
#include <optional>
#include <range/v3/algorithm.hpp>
#include <range/v3/numeric.hpp>
#include <range/v3/view.hpp>

namespace rng = ranges;

constexpr auto DEFAULT_RES_FILENAME = "millepede.res";
constexpr auto SCALE_FACTOR = 10.F;
constexpr auto REFERENCE_BAR_NUM = 25;
constexpr auto MILLE_BUFFER_SIZE = std::size_t{ 100000 };

namespace R3B::Neuland::Calibration
{
    namespace
    {
        void calculate_time_offset(R3B::Neuland::Cal2HitPar& cal_to_hit_par)
        {
            auto& module_pars = cal_to_hit_par.GetListOfModuleParRef();
            for (auto& [module_num, module_par] : module_pars)
            {
                if (module_par.effectiveSpeed.value != 0)
                {
                    module_par.tDiff = module_par.offset_effective_c / module_par.effectiveSpeed;
                }
            }
        }

        auto calculate_time_difference(const R3B::Neuland::BarCalData& bar_cal_data) -> R3B::ValueError<double>
        {
            const auto& left_signal = bar_cal_data.left.front();
            const auto& right_signal = bar_cal_data.right.front();
            const auto t_diff = (right_signal.leading_time - right_signal.trigger_time) -
                                (left_signal.leading_time - left_signal.trigger_time);
            return t_diff;
        }

        inline auto has_only_one_signal(const BarCalData& bar_signal) -> bool
        {
            return bar_signal.left.size() == 1 and bar_signal.right.size() == 1;
        }

    } // namespace

    void MillepedeEngine::fill_time_differences(const CalData& cal_data, TH2I* hist_time_offsets)
    {
        auto max_module_num = std::min(GetMaxModuleNum(), GetModuleSize());
        for (const auto& [plane_num, plane_cal_data] : cal_data)
        {
            for (const auto& [bar_num, bar_signal] : plane_cal_data.bar_cal_data)
            {
                if (bar_signal.module_num > max_module_num or not has_only_one_signal(bar_signal))
                {
                    continue;
                }
                auto time_difference = calculate_time_difference(bar_signal);
                hist_time_offsets->Fill(bar_signal.module_num, time_difference.value);
            }
        }
    }

    void MillepedeEngine::Init()
    {
        cal_to_hit_par_ = GetTask()->GetCal2HitPar();
        mille_handler_.init();
        mille_handler_.set_max_number_of_modules(std::min(GetMaxModuleNum(), GetModuleSize()));

        par_result_.set_filename(DEFAULT_RES_FILENAME);
        pede_launcher_.set_steer_filename(pede_steer_filename_);
        pede_launcher_.set_parameter_filename(parameter_filename_);
        // binary_data_writer_.set_buffer_size(MILLE_BUFFER_SIZE);
    }

    void MillepedeEngine::fill_module_parameters(const Millepede::ResultReader& result,
                                                 Neuland::Cal2HitPar& cal_to_hit_par)
    {
        const auto& pars = result.get_pars();
        const auto max_module_num = std::min(GetModuleSize(), GetMaxModuleNum());
        for (const auto& [par_id, par] : pars)
        {
            const auto [module_num, global_label] = to_module_num_label(par_id, max_module_num);
            auto& module_pars = cal_to_hit_par.GetListOfModuleParRef();

            auto& par_ref = module_pars.emplace(module_num, HitModulePar{}).first->second;
            switch (global_label)
            {
                case GlobalLabel::tsync:
                    par_ref.tSync += ValueErrorD{ par.value, par.error } * SCALE_FACTOR;
                    break;
                case GlobalLabel::offset_effective_c:
                    // The value here is the product of tDiff and effectiveSped. Real tDiff will be calculated later
                    par_ref.offset_effective_c += (ValueErrorD{ par.value, par.error } * SCALE_FACTOR);
                    break;
                case GlobalLabel::effective_c:
                    par_ref.effectiveSpeed += ValueErrorD{ par.value, par.error };
                    break;
                default:
                    throw std::runtime_error("An error occured with unrecognized global tag");
            }
        }
        calculate_time_offset(*cal_to_hit_par_);
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
            case State::histogram_calibration:
                fill_time_differences(signals, hist_time_offsets_);
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
        // for (const auto& [plane_num, plane_signals] : signals)
        // {
        //     for (const auto& [bar_num, bar_signal] : plane_signals.bar_cal_data)
        //     {
        //         add_spacial_local_constraint(plane_num, GetBarVerticalDisplacement(bar_num), BarSize_XY / 2.);
        //         if (has_only_one_signal(bar_signal))
        //         {
        //             add_signal_t(bar_signal);
        //         }
        //     }
        // }
        // write_local_constraint(bar_positions_);

        // return;

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
            case State::histogram_calibration:
                histogram_calibrate(hit_par);
                break;
            case State::millepede_calibration:
                millepede_calibrate(hit_par);
                break;
        }
    }

    void MillepedeEngine::histogram_calibrate(Cal2HitPar& cal_to_hit_par)
    {
        R3BLOG(info, "Begin histogram calibration");
        auto max_module_num = std::min(GetModuleSize(), GetMaxModuleNum());
        auto& module_pars = cal_to_hit_par.GetListOfModuleParRef();
        for (int module_num{ 1 }; module_num <= max_module_num; ++module_num)
        {
            const auto [time_offset, effective_c] = calculate_time_offset_effective_speed(module_num);

            auto& par_ref = module_pars.emplace(module_num, HitModulePar{}).first->second;
            par_ref.tDiff.value = time_offset.value;
            par_ref.tDiff.error = time_offset.error;
            par_ref.effectiveSpeed.value = effective_c.value;
            par_ref.effectiveSpeed.error = effective_c.error;

            const auto offset_effective_c = time_offset * effective_c;
            par_ref.offset_effective_c.value = offset_effective_c.value;
            par_ref.offset_effective_c.error = offset_effective_c.error;
        }

        current_state_ = State::millepede_calibration;
        init_steer_writer(cal_to_hit_par);
        R3BLOG(info, "Histogram calibration finished. Ready for millepede calibration!");
    }

    void MillepedeEngine::millepede_calibrate(Cal2HitPar& cal_to_hit_par)
    {
        R3BLOG(info, "Starting millepede calibration...");
        mille_handler_.close();
        // binary_data_writer_.close();

        R3BLOG(info, "Launching pede algorithm...");
        pede_launcher_.sync_launch();
        pede_launcher_.end();

        par_result_.read();
        fill_module_parameters(par_result_, cal_to_hit_par);
        R3BLOG(info, "Millepede calibration finished.");
    }

    void MillepedeEngine::EndOfEvent(unsigned int /*event_num*/)
    {
        // TODO: could be an empty event
        mille_handler_.do_end_of_event();
        // binary_data_writer_.end();
    }

    void MillepedeEngine::EventReset() {}

    void MillepedeEngine::HistInit(DataMonitor& histograms)
    {
        const auto module_size = GetModuleSize();

        constexpr auto hist_time_offsets_y_size = 1000;
        constexpr auto hist_time_offsets_y_max = 500;
        hist_time_offsets_ = histograms.add_hist<TH2I>("hist_time_offsets",
                                                       "hist_time_offsets",
                                                       module_size,
                                                       0.5,
                                                       module_size + 0.5,
                                                       hist_time_offsets_y_size,
                                                       -hist_time_offsets_y_max,
                                                       hist_time_offsets_y_max);
    }

    void MillepedeEngine::write_to_buffer() { mille_handler_.write_to_buffer(); }

    void MillepedeEngine::EndOfTask()
    {
        average_t_sum_.reset();
        mille_handler_.buffer_clear();
    }

    void MillepedeEngine::init_steer_writer(const Cal2HitPar& /*cal_to_hit_par*/)
    {
        auto steer_writer = SteerWriter{};
        steer_writer.set_filepath(pede_steer_filename_);
        steer_writer.set_parameter_file(parameter_filename_);
        steer_writer.set_data_filepath(mille_handler_.get_filename());
        steer_writer.add_method(SteerWriter::Method::inversion,
                                std::make_pair(static_cast<float>(pede_interartion_number_), pede_error_threshold_));
        steer_writer.add_other_options(std::vector<std::string>{ "outlierdownweighting", "4" });
        // steer_writer.add_other_options(
        //     std::vector<std::string>{ "scaleerrors", fmt::format("{}", error_scale_factor_) });

        // const auto& module_pars = cal_to_hit_par.GetModulePars();
        // for (const auto& [module_num, module_par] : module_pars)
        // {
        //     steer_writer.add_parameter_default(
        //         to_global_label_id(static_cast<int>(module_num), GlobalLabel::effective_c),
        //         std::make_pair(static_cast<float>(module_par.effectiveSpeed.value), 0.F));
        //     steer_writer.add_parameter_default(
        //         to_global_label_id(static_cast<int>(module_num), GlobalLabel::offset_effective_c),
        //         std::make_pair(
        //             static_cast<float>(module_par.effectiveSpeed.value * module_par.tDiff.value / SCALE_FACTOR),
        //             0.F));
        // }
        steer_writer.write();
    }

    auto MillepedeEngine::calculate_time_offset_effective_speed(int module_num) -> std::pair<ValueErrorD, ValueErrorD>
    {
        auto bin_num = hist_time_offsets_->GetXaxis()->FindBin(module_num);
        auto* bar_dist_time_diff = hist_time_offsets_->ProjectionY("_y", bin_num, bin_num);
        if (bar_dist_time_diff == nullptr)
        {
            throw R3B::runtime_error(
                fmt::format("Cannot project histogram in y direction at module num {}", module_num));
        }
        // Normalization:
        bar_dist_time_diff->Scale(1. / bar_dist_time_diff->GetEntries());

        auto* bar_dist_time_diff_CDF = bar_dist_time_diff->GetCumulative();

        // Range determination
        constexpr auto low_value_cut = 0.05;
        constexpr auto high_value_cut = 0.95;

        auto low_thres = bar_dist_time_diff_CDF->GetBinCenter(bar_dist_time_diff_CDF->FindFirstBinAbove(low_value_cut));
        auto high_thres =
            bar_dist_time_diff_CDF->GetBinCenter(bar_dist_time_diff_CDF->FindFirstBinAbove(high_value_cut));

        // fit
        // options: S: return the result Q: quite mode 0: not canvas drawing
        auto fit_result = bar_dist_time_diff_CDF->Fit("pol1", "SQ0", "", low_thres, high_thres);

        // TODO: error analysis better needed here
        auto offset = ValueErrorD{ fit_result->Parameter(0), fit_result->Error(0) };
        auto slope = ValueErrorD{ fit_result->Parameter(1), fit_result->Error(1) };

        auto effective_speed = slope * BarLength * 2.;

        // TODO: the error calculation here is wrong as two values are not independent
        auto time_offset = (ValueErrorD{ 0.5, 0 } - offset) / slope;

        return std::make_pair(time_offset, effective_speed);
    }
} // namespace R3B::Neuland::Calibration
