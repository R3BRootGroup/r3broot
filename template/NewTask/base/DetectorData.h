/******************************************************************************
 *   Copyright (C) 2019 GSI Helmholtzzentrum für Schwerionenforschung GmbH    *
 *   Copyright (C) 2019-2025 Members of R3B Collaboration                     *
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

#include <TNamed.h>

// if needed
#include <fmt/format.h>

namespace R3B
{
    struct DetectorData : public TNamed
    {
      public:
        explicit DetectorData(std::string_view name = "DefaultName")
            : TNamed(name.data(), ""){};

        int module_num = 0; // starts from 1
        int module_id = 0;  // starts from 0
        double energy = 0.;
        uint32_t tpat = 0U;

        void reset()
        {
            module_num = 0;
            module_id = 0;
            energy = 0.;
            tpat = 0U;
        }
        ClassDefOverride(DetectorData, 1);
    };

    struct DetectorCalData : public TObject
    {
      public:
        DetectorCalData() = default;
        int module_num = 0; // starts from 1
        double left_time = 0.;
        double right_time = 0.;
        ClassDefOverride(DetectorCalData, 1);
    };

    struct DetectorHitData : public TObject
    {
      public:
        DetectorHitData() = default;
        int module_num = 0; // starts from 1
        double energy = 0.;
        double time = 0.;
        double position = 0.;
        ClassDefOverride(DetectorHitData, 1);
    };
} // namespace R3B

// optional for the nice printing
template <>
class fmt::formatter<R3B::DetectorData>
{
  public:
    static constexpr auto parse(format_parse_context& ctx) { return ctx.end(); }
    template <typename FmtContent>
    constexpr auto format(const R3B::DetectorData& data, FmtContent& ctn) const
    {
        return format_to(
            ctn.out(), "{{module_id: {}, energy: {}, tpat: {:x}}}", data.module_id, data.energy, data.tpat);
    }
};
