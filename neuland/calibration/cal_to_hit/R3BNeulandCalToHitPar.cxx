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

#include "R3BNeulandCalToHitPar.h"
#include <R3BNeulandCommon.h>
#include <R3BShared.h>
#include <algorithm>
#include <cmath>
#include <range/v3/all.hpp>

namespace
{
    auto calculate_threshold(const R3B::Neuland::HitModulePar& modulePar, R3B::Side side)
    {
        const auto value = modulePar.PMTThreshold.get(side).value *
                           std::exp(R3B::Neuland::TotalBarLength / modulePar.lightAttenuationLength.value);
        const auto err = 0.;
        return R3B::ValueError<double>{ value, err };
    }

} // namespace

namespace R3B::Neuland
{
    auto Cal2HitPar::GetModuleParAt(unsigned int module_num) const -> const HitModulePar&
    {
        if (auto module_par = module_pars_.find(module_num); module_par == module_pars_.end())
        {
            throw R3B::runtime_error(fmt::format("Cannot find the parameter with module num {}", module_num));
        }
        return module_pars_.at(module_num);
    }

    Cal2HitPar::Cal2HitPar(std::string_view name, std::string_view title, std::string_view context, Bool_t own)
        : ParSet(name, title, context, own)
    {
    }

    void Cal2HitPar::CalculateEnergyCutoff()
    {
        auto threshold_views =
            ranges::views::transform(module_pars_ | ranges::views::values,
                                     [](const auto& ele)
                                     {
                                         return std::max(calculate_threshold(ele, R3B::Side::left).value,
                                                         calculate_threshold(ele, R3B::Side::right).value);
                                     });

        auto max_value = ranges::max_element(threshold_views);
        if (max_value != threshold_views.end())
        {
            energy_cut_ = *max_value;
        }
    }

} // namespace R3B::Neuland

ClassImp(R3B::Neuland::Cal2HitPar);
