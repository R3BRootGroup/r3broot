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
    constexpr auto __sqrt12 = 3.464101615;

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

constexpr UInt_t nEvents = (1U << 16);

void testNeulandBar()
{
    TRandom3 rnd(0);

    const double tdiff = rnd.Uniform(-1000., 1000.);
    const double veff = rnd.Gaus(Neuland::AvgEffectiveCLight, 0.2);
    const double invveff = 1. / veff;
    const std::array<double, 2> gain = { rnd.Gaus(Neuland::AvgGain, 0.2), rnd.Gaus(Neuland::AvgGain, 0.2) };
    const double attlen = rnd.Gaus(Neuland::AvgAttenuationLength, 25.);
    const double invattlen = 1. / attlen;

    auto bar = new Neuland::Calibration::HitCalibrationBar(0);

    for (auto event = 0; event < nEvents; ++event)
    {
        const auto position = rnd.Uniform(-0.5 * Neuland::BarLength, 0.5 * Neuland::BarLength);
        const std::array<double, 2> distance = { 0.5 * Neuland::TotalBarLength + position,
                                                 0.5 * Neuland::TotalBarLength - position };
        const std::array<double, 2> time = { distance[0] * 0.5 * fabs(invveff) + tdiff * 0.5,
                                             distance[1] * 0.5 * fabs(invveff) - tdiff * 0.5 };

        const auto eloss = rnd.Uniform(1., 20.);
        const std::array<double, 2> eatpmt = { eloss * std::exp(-invattlen * distance[0]),
                                               eloss * std::exp(-invattlen * distance[1]) };

        const auto pos = rnd.Gaus(position, 1.);
        const std::array<double, 2> tdc = { rnd.Gaus(time[0], 0.05), rnd.Gaus(time[1], 0.05) };
        const std::array<double, 2> qdc = { gain[0] * eatpmt[0] / (1. + Neuland::AvgSaturationConstant * eatpmt[0]),
                                            gain[1] * eatpmt[1] / (1. + Neuland::AvgSaturationConstant * eatpmt[1]) };
        const auto energy = rnd.Gaus(eloss, eloss * 0.05);

        if (energy <= 0.)
            continue;

        for (auto side = 0; side < 2; ++side)
            bar->Set(side, tdc[side], qdc[side]);

        bar->Add(0., pos, pos, energy);
    }

    bar->Calibrate();

    const auto pars = bar->GetParameters();
    std::cout << endl << pars.GetTimeOffset(2) - pars.GetTimeOffset(1) << "  " << tdiff << std::endl;
    std::cout << pars.GetEffectiveSpeed() << "  " << veff << std::endl;
    std::cout << pars.GetEnergyGain(1) << "  " << gain[0] << std::endl;
    std::cout << pars.GetEnergyGain(2) << "  " << gain[1] << std::endl;
    std::cout << pars.GetLightAttenuationLength() << "  " << attlen << std::endl;

    if (fabs(pars.GetTimeOffset(2) - pars.GetTimeOffset(1) - tdiff) < 0.01 &&
        fabs(pars.GetEffectiveSpeed() - veff) < 0.1 && fabs(pars.GetEnergyGain(1) - gain[0]) < 0.1 &&
        fabs(pars.GetEnergyGain(2) - gain[1]) < 0.01 && fabs(pars.GetLightAttenuationLength() - attlen) < 25)
    {
        std::cout << "SUCCESS" << std::endl;
    }
    else
    {
        std::cout << "FAILED" << std::endl;
    }
}
