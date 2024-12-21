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
#include "R3BDigitizingPaddle.h"
#include <gmock/gmock.h>

/**
 * Simulation of Mock Bar/Paddle
 * @author Yanzhao Wang
 * @since April 12, 2023
 */
namespace R3B::testing::Neuland
{
    using Digitizing::Channel;
    using Digitizing::ChannelSide;
    using Digitizing::Paddle;

    struct GMockNeulandPaddle : public Paddle
    {
        explicit GMockNeulandPaddle(uint16_t paddleID)
            : Paddle{ paddleID }
        {
        }
        MOCK_METHOD(double,
                    ComputeTime,
                    (const Channel::Signal& firstSignal, const Channel::Signal& secondSignal),
                    (const, override));
        MOCK_METHOD(double,
                    ComputeEnergy,
                    (const Channel::Signal& firstSignal, const Channel::Signal& secondSignal),
                    (const, override));
        MOCK_METHOD(double,
                    ComputePosition,
                    (const Channel::Signal& firstSignal, const Channel::Signal& secondSignal),
                    (const, override));
        MOCK_METHOD(Pair<Channel::Hit>, ComputeChannelHits, (const Hit& hit), (const, override));
    };

    struct GMockChannel : public Channel
    {
        explicit GMockChannel(ChannelSide side)
            : Digitizing::Channel{ side } {};
        void AddHit(Hit hit) override { hits_.push_back(hit); }

        MOCK_METHOD(void, AttachToPaddle, (Paddle*), (override));
        MOCK_METHOD(Signals, ConstructSignals, (), (override));
        std::vector<Hit> hits_;
    };
} // namespace R3B::testing::Neuland
