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

// FIXME Root Problem remove this as soon as Fairsoft is using Root >=6.20.00
namespace Neuland
{
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
} // namespace Neuland

constexpr Int_t nPlanes = 2;
constexpr Int_t nEvents = 1024;

void testNeulandTSyncer()
{
    TRandom3 rnd(23);

    std::array<std::array<Double_t, Neuland::BarsPerPlane>, nPlanes> tsync;
    Neuland::Calibration::TSyncer* syncer = new Neuland::Calibration::TSyncer();

    for (auto plane = 0; plane < nPlanes; ++plane)
    {
        for (auto bar = 0; bar < Neuland::BarsPerPlane; ++bar)
        {
            tsync[plane][bar] = rnd.Uniform(-100., 100.);
        }
    }

    for (auto event = 0; event < nEvents; ++event)
    {
        const auto offset = rnd.Uniform(-100., 100.);
        for (auto plane = 0; plane < nPlanes; ++plane)
        {
            for (auto bar = 0; bar < Neuland::BarsPerPlane; ++bar)
            {
                syncer->AddBarData(plane * Neuland::BarsPerPlane + bar, rnd.Gaus(tsync[plane][bar], 0.05) + offset);
            }
        }
        syncer->DoEvent();
    }

    const auto result = syncer->GetTSync(nPlanes);
    const auto offset = result[0].Value - tsync[0][0];
    for (auto plane = 0; plane < nPlanes; ++plane)
    {
        for (auto bar = 0; bar < Neuland::BarsPerPlane; ++bar)
        {
            if (!std::isfinite(result[plane * Neuland::BarsPerPlane + bar].Value) ||
                fabs(result[plane * Neuland::BarsPerPlane + bar].Value - tsync[plane][bar] - offset) > 0.05)
            {
                std::cout << "FAILED" << std::endl;
                return;
            }
        }
    }

    std::cout << "SUCCESS" << std::endl;
}
