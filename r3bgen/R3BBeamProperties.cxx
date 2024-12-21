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

#include "R3BBeamProperties.h"
#include "R3BDistribution2D.h"
#include "R3BDistribution3D.h"

#include "TGeoPhysicalConstants.h" // for  amu_c2

constexpr auto AMU_MeV = TGeoUnit::amu_c2 / TGeoUnit::MeV; // MeV

R3BBeamProperties::R3BBeamProperties()
    : fVertex_cm(R3BDistribution3D::Delta(0, 0, 0))
    , fSpread_mRad(R3BDistribution2D::Delta(0, 0))
    , fBeta(R3BDistribution1D::Delta(0.5))
{
}

void R3BBeamProperties::SetEnergyDistribution(R3BDistribution<1> energyDistribution_AMeV,
                                              int nNucleons,
                                              double mass_MeV)
{
    const auto nucleonMass = (mass_MeV == 0. ? nNucleons * AMU_MeV : mass_MeV);

    fBeta = R3BDistribution<1>(
        [energyDistribution_AMeV, nucleonMass, nNucleons](std::array<double, 1> values) mutable -> std::array<double, 1>
        {
            const auto energy = energyDistribution_AMeV.GetRandomValues(values)[0];
            const auto gamma = (nNucleons * energy) / nucleonMass + 1.;
            const auto beta = std::sqrt(1. - 1. / (gamma * gamma));
            return { beta };
        });
}
