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
#pragma once

#include "R3BNeulandCalData2.h"
#include "R3BNeulandCalToHitPar.h"

namespace R3B::Neuland::Calibration
{
    enum class GlobalLabel
    {
        tsync,              // tsync
        offset_effective_c, // offset times effective_C
        effective_c         // effective speed of light
    };
    inline auto calculate_position_along_bar(const BarCalData& signal, Cal2HitPar* cal_to_hit_par) -> double
    {
        const auto module_num = static_cast<int>(signal.module_num);
        auto effective_c = cal_to_hit_par->GetModuleParAt(module_num).effectiveSpeed.value;
        auto time_offset = cal_to_hit_par->GetModuleParAt(module_num).tDiff.value;
        auto offset_effective_c = effective_c * time_offset;

        const auto& left_signal = signal.left.front();
        const auto& right_signal = signal.right.front();
        const auto t_diff = (right_signal.leading_time - right_signal.trigger_time) -
                            (left_signal.leading_time - left_signal.trigger_time);

        return (offset_effective_c - t_diff.value * effective_c) / 2.;
    }

    inline auto to_module_num_label(int par_num, int max_module_num) -> std::pair<int, GlobalLabel>
    {
        auto res = std::pair<int, GlobalLabel>{};
        const auto factor = (par_num - 1) / max_module_num;
        res.first = (par_num - 1) % max_module_num + 1;

        switch (factor)
        {
            // case 0:
            //     res.second = GlobalLabel::tsync;
            //     break;
            // case 1:
            //     res.second = GlobalLabel::offset_effective_c;
            //     break;
            // case 2:
            //     res.second = GlobalLabel::effective_c;
            //     break;
            case 0:
                res.second = GlobalLabel::offset_effective_c;
                break;
            case 1:
                res.second = GlobalLabel::effective_c;
                break;
            default:
                throw R3B::logic_error(fmt::format("An error occured with unrecognized global par id: {}", par_num));
        }

        return res;
    }

    inline auto to_global_label_id(int module_num, GlobalLabel label, int max_module_num) -> int
    {
        switch (label)
        {
            // case GlobalLabel::tsync:
            //     return module_num;
            // case GlobalLabel::offset_effective_c:
            //     return module_num + num_of_module;
            // case GlobalLabel::effective_c:
            //     return module_num + 2 * num_of_module;
            case GlobalLabel::offset_effective_c:
                return module_num;
            case GlobalLabel::effective_c:
                return module_num + max_module_num;
            case GlobalLabel::tsync:
                return 0;
            default:
                throw std::runtime_error("An error occured with unrecognized global tag");
        }
    }
} // namespace R3B::Neuland::Calibration
