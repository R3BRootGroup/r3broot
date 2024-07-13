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

#include "R3BNeulandPositionCal.h"
#include "Utilities.h"
#include <R3BDataMonitor.h>
#include <R3BException.h>
#include <R3BNeulandCommon.h>
#include <TFitResult.h>
#include <TH2I.h>
#include <fmt/format.h>

namespace R3B::Neuland::Calibration
{
    namespace
    {
        auto calculate_time_difference(const Neuland::BarCalData& bar_cal_data) -> ValueError<double>
        {
            const auto& left_signal = bar_cal_data.left.front();
            const auto& right_signal = bar_cal_data.right.front();
            const auto t_diff = (right_signal.leading_time - right_signal.trigger_time) -
                                (left_signal.leading_time - left_signal.trigger_time);
            return t_diff;
        }
    } // namespace

    void PositionCalibrator::HistInit(DataMonitor& histograms)
    {
        constexpr auto hist_time_offsets_y_size = 1000;
        constexpr auto hist_time_offsets_y_max = 500;
        hist_time_offsets_ = histograms.add_hist<TH2I>("hist_time_offsets",
                                                       "hist_time_offsets",
                                                       max_number_of_modules_,
                                                       0.5,
                                                       max_number_of_modules_ + 0.5,
                                                       hist_time_offsets_y_size,
                                                       -hist_time_offsets_y_max,
                                                       hist_time_offsets_y_max);
    }

    void PositionCalibrator::fill_time_differences(const CalData& cal_data)
    {
        for (const auto& [plane_num, plane_cal_data] : cal_data)
        {
            for (const auto& [bar_num, bar_signal] : plane_cal_data.bar_cal_data)
            {
                if (bar_signal.module_num > max_number_of_modules_ or not has_only_one_signal(bar_signal))
                {
                    continue;
                }
                auto time_difference = calculate_time_difference(bar_signal);
                hist_time_offsets_->Fill(bar_signal.module_num, time_difference.value);
            }
        }
    }

    void PositionCalibrator::calibrate(Cal2HitPar& cal_to_hit_par)
    {
        auto& module_pars = cal_to_hit_par.GetListOfModuleParRef();
        for (int module_num{ 1 }; module_num <= max_number_of_modules_; ++module_num)
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
    }

    auto PositionCalibrator::calculate_time_offset_effective_speed(int module_num)
        -> std::pair<ValueErrorD, ValueErrorD>
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
        // options: S: return the result; Q: quite mode; 0: no canvas drawing
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
