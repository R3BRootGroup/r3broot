#pragma once

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

#include <R3BFormatters.h>
#include <R3BShared.h>
#include <TObject.h>
#include <numeric>
#include <range/v3/view.hpp>
#include <vector>

namespace R3B::Neuland
{

    // TODO: Provide both flattened data and the composite
    struct CalDataSignal
    {
        ValueError<double> leading_time;        // ns
        ValueError<double> time_over_threshold; // ns
        ValueError<double> trigger_time;        // ns
        ClassDefNV(CalDataSignal, 1)
    };

    struct BarCalData
    {
      public:
        BarCalData() = default;
        explicit BarCalData(int mod_num)
            : module_num{ mod_num }
        {
        }
        int module_num = 0; // 1 based bar num
        std::vector<CalDataSignal> left;
        std::vector<CalDataSignal> right;
        ClassDefNV(BarCalData, 1)
    };

    struct PlaneCalData
    {
      public:
        PlaneCalData() = default;
        explicit PlaneCalData(int num)
            : plane_num{ num }
        {
        }

        int plane_num = 0;
        std::unordered_map<int, BarCalData> bar_cal_data;
        ClassDefNV(PlaneCalData, 1)
    };

    using CalData = std::unordered_map<int, PlaneCalData>;

    inline auto GetBarCalDataSize(const CalData& cal_data) -> int
    {
        return std::accumulate(cal_data.begin(),
                               cal_data.end(),
                               int{},
                               [](int init, const auto& plane_cal_data)
                               { return init + plane_cal_data.second.bar_cal_data.size(); });
    }

} // namespace R3B::Neuland

template <>
class fmt::formatter<R3B::Neuland::CalDataSignal>
{
  public:
    static constexpr auto parse(format_parse_context& ctx) { return ctx.end(); }
    template <typename FmtContent>
    constexpr auto format(const R3B::Neuland::CalDataSignal& signal, FmtContent& ctn) const
    {
        return format_to(ctn.out(),
                         "{{leadT: {}, tot: {}, trigT: {} }}",
                         signal.leading_time,
                         signal.time_over_threshold,
                         signal.trigger_time);
    }
};

template <>
class fmt::formatter<R3B::Neuland::BarCalData>
{
  public:
    static constexpr auto parse(format_parse_context& ctx) { return ctx.end(); }
    template <typename FmtContent>
    constexpr auto format(const R3B::Neuland::BarCalData& signal, FmtContent& ctn) const
    {
        return format_to(ctn.out(),
                         "ModuleNum: {}, left bar: [{}], right bar: [{}]",
                         signal.module_num,
                         fmt::join(signal.left, ", "),
                         fmt::join(signal.right, ", "));
    }
};

template <>
class fmt::formatter<R3B::Neuland::PlaneCalData>
{
  public:
    static constexpr auto parse(format_parse_context& ctx) { return ctx.end(); }
    template <typename FmtContent>
    constexpr auto format(const R3B::Neuland::PlaneCalData& signal, FmtContent& ctn) const
    {
        return format_to(ctn.out(),
                         "PlaneNum: {}, bar signals: \n\t{}",
                         signal.plane_num,
                         fmt::join(signal.bar_cal_data | ranges::views::values, "\n\t"));
    }
};
