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
#include "R3BDigitizingEngine.h"

/**
 * Simulation of NeuLAND Bar/Paddle
 * @author Yanzhao Wang
 * @since July 15, 2022
 */

namespace R3B::Digitizing::Neuland
{
    class NeulandPaddle : public Digitizing::Paddle
    {
      public:
        explicit NeulandPaddle(uint16_t paddleID);

      private:
        [[nodiscard]] auto ComputeTime(const Channel::Signal& firstSignal, const Channel::Signal& secondSignal) const
            -> double override;
        [[nodiscard]] auto ComputeEnergy(const Channel::Signal& firstSignal, const Channel::Signal& secondSignal) const
            -> double override;
        [[nodiscard]] auto ComputePosition(const Channel::Signal& leftSignal, const Channel::Signal& rightSignal) const
            -> double override;
        auto ComputeChannelHits(const Hit& hit) const -> Pair<Channel::Hit> override;

      public:
        static constexpr double gHalfLength = 135.;   // [cm]
        static constexpr double gCMedium = 14.;       // speed of light in material in [cm/ns]
        static constexpr double gAttenuation = 0.008; // light attenuation of plastic scintillator [1/cm]
        static constexpr double gLambda = 1. / 2.1;
        static const double ReverseAttenFac;

        static auto MatchSignals(const Channel::Signal& firstSignal, const Channel::Signal& secondSignal) -> float;
        [[nodiscard]] static auto SignalCouplingNeuland(const Channel::Signals& firstSignals,
                                                        const Channel::Signals& secondSignals)
            -> std::vector<ChannelSignalPair>;
        static auto GenerateChannelHit(Double_t mcTime, Double_t mcLight, Double_t dist) -> Channel::Hit;
    };
} // namespace R3B::Digitizing::Neuland
