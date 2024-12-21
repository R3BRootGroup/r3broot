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

#ifndef NEULANDCOMMON_H
#define NEULANDCOMMON_H

#include <cmath>
#include <limits>
#include <vector>

namespace Neuland
{
    // Constants

    constexpr auto __sqrt12 = 3.464101615;
    constexpr auto NaN = std::numeric_limits<double>::quiet_NaN();
    constexpr auto Inf = std::numeric_limits<double>::infinity();

    // Usefull functions I do not want to redefine in every file

    template <typename T>
    constexpr T Sqr(const T val)
    {
        return val * val;
    }

    // Initialize variables from Birk' s Law

    constexpr double BirkdP = 1.032;
    constexpr double BirkC1 = 0.013 / BirkdP;
    constexpr double BirkC2 = 9.6e-6 / (BirkdP * BirkdP);

    constexpr double GetLightYieldAfterBirk(const int charge, const double length_mm, const double edep_MeV)
    {
        return (charge * length_mm == 0.
                    ? edep_MeV
                    : edep_MeV / (1. + BirkC1 * (Sqr(charge) > 1 ? 7.2 / 12.6 : 1.) * (edep_MeV / length_mm) +
                                  BirkC2 * Sqr(edep_MeV / length_mm)));
    }

    // Physical Constants

    constexpr auto CLight = 29.9792458;     // Speed of light [cm/ns]
    constexpr auto InvCLight = 1. / CLight; // Speed of light [cm/ns]>

    // Electronics Constans

    constexpr auto MaxCalTime = 5. * 2048;

    // Geometry & Material Constants

    constexpr auto BarSize_XY = 5.0;                                  // cm NeuLAND parameter
    constexpr auto BarUncertainty_XY = BarSize_XY / __sqrt12;         // cm NeuLAND parameter
    constexpr auto BarSize_Z = 5.0;                                   // cm NeuLAND parameter
    constexpr auto BarUncertainty_Z = BarSize_Z / __sqrt12;           // cm NeuLAND parameter
    constexpr auto BarLength = 250.0;                                 // cm NeuLAND parameter
    constexpr auto LightGuideLength = 10.0;                           // cm NeuLAND parameter
    constexpr auto TotalBarLength = BarLength + 2 * LightGuideLength; // cm NeuLAND parameter, Bar including Light Guide

    constexpr auto ScintillatorDensity = 1.032;        // g / cm^3
    constexpr auto MIPStoppingPowerPerDensity = 1.956; // MeV cm^2 / g
    constexpr auto MIPStoppingPower = 1.73;            // MeV / cm

    constexpr auto FirstHorizontalPlane = 0;
    constexpr auto BarsPerPlane = 50;
    constexpr auto MaxNumberOfPlanes = 60;
    constexpr auto MaxNumberOfBars = MaxNumberOfPlanes * BarsPerPlane;

    constexpr bool IsPlaneHorizontal(const int plane) { return (plane % 2 == FirstHorizontalPlane); }
    constexpr bool IsPlaneVertical(const int plane) { return !IsPlaneHorizontal(plane); }
    constexpr int GetPlaneNumber(const int barID) { return barID / BarsPerPlane; }

    // Average Parameters

    constexpr auto AvgTimeResolution = 0.150;   // ns
    constexpr auto AvgEffectiveCLight = -7.95;  // cm / ns
    constexpr auto AvgGain = 15;                // MeV / ns
    constexpr auto AvgThreshold = 1.75;         // MeV
    constexpr auto AvgAttenuationLength = 400.; // cm

    constexpr auto SaturationCoefficient = 1.75e-3; // 1 / ns

} // namespace Neuland

#endif
