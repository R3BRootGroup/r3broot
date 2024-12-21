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

#include "R3BDigitizingPaddle.h"
#include "FairLogger.h"

namespace R3B::Digitizing
{
    Paddle::Paddle(int paddleID, SignalCouplingStrategy strategy)
        : fPaddleID{ paddleID }
        , fSignalCouplingStrategy{ std::move(strategy) }
    {
    }

    void Paddle::SetChannel(std::unique_ptr<Channel> channel)
    {
        channel->SetPaddle(this);
        channel->AttachToPaddle(this);
        if (channel->GetSide() == ChannelSide::left)
        {
            fLeftChannel = std::move(channel);
        }
        else if (channel->GetSide() == ChannelSide::right)
        {
            fRightChannel = std::move(channel);
        }
        else
        {
            LOG(error) << "Channel side is invalid!";
        }
    }

    void Paddle::DepositLight(const Hit& hit)
    {
        auto channelHits = ComputeChannelHits(hit);
        fLeftChannel->AddHit(channelHits.left);
        fRightChannel->AddHit(channelHits.right);
    }

    auto Paddle::HasFired() const -> bool
    {
        if (!fLeftChannel || !fRightChannel)
        {
            LOG(fatal) << "channels failed to be constructed when checking fire! ";
            return false;
        }
        return (fLeftChannel->HasFired() && fRightChannel->HasFired());
    }

    auto Paddle::HasHalfFired() const -> bool
    {
        return (fLeftChannel->HasFired() && !fRightChannel->HasFired()) ||
               (!fLeftChannel->HasFired() && fRightChannel->HasFired());
    }

    auto Paddle::ConstructPaddelSignals(const Channel::Signals& firstSignals,
                                        const Channel::Signals& secondSignals) const -> Signals
    {
        auto channelSignalPairs = fSignalCouplingStrategy(firstSignals, secondSignals);

        auto paddleSignals = std::vector<Signal>();
        paddleSignals.reserve(channelSignalPairs.size());

        for (auto& it : channelSignalPairs)
        {
            // warning if two channel signals are from the same side
            if (it.left.get().side == it.right.get().side)
            {
                LOG(warn) << "DigitizingEngine: ConstructPaddleSignals(): cannot construct the paddle signal from the "
                             "channels on the same side!";
            }

            // swap the channel signals if necessary
            if (it.left.get().side != ChannelSide::left)
            {
                std::swap(it.left, it.right);
            }

            auto paddleSignal = Signal{ { it.left, it.right } };
            paddleSignal.energy = ComputeEnergy(it.left, it.right);
            paddleSignal.time = ComputeTime(it.left, it.right);
            paddleSignal.position = ComputePosition(it.left, it.right);
            paddleSignals.push_back(paddleSignal);
        }
        return paddleSignals;
    }

    auto Paddle::GetSignals() const -> const Signals&
    {
        if (!fSignals.valid())
        {
            if (HasFired())
            {
                auto signals = ConstructPaddelSignals(fLeftChannel->GetSignals(), fRightChannel->GetSignals());
                fSignals.set(std::move(signals));
            }
            else
            {
                fSignals.set({});
            }
        }

        return fSignals.getRef();
    }

    auto Paddle::SignalCouplingByTime(const Channel::Signals& firstSignals, const Channel::Signals& secondSignals)
        -> std::vector<ChannelSignalPair>
    {
        auto firstSignalRefs =
            std::vector<std::reference_wrapper<const Channel::Signal>>(firstSignals.begin(), firstSignals.end());
        auto secondSignalRefs =
            std::vector<std::reference_wrapper<const Channel::Signal>>(secondSignals.begin(), secondSignals.end());
        auto channelSignalPairs = std::vector<ChannelSignalPair>{};
        channelSignalPairs.reserve(std::min(firstSignals.size(), secondSignals.size()));

        auto sortSignalByTime = [](auto signals)
        {
            std::sort(signals.begin(),
                      signals.end(),
                      [](auto& first, auto& second) { return first.get().tdc < second.get().tdc; });
        };
        sortSignalByTime(firstSignalRefs);
        sortSignalByTime(secondSignalRefs);

        for (auto first_it = firstSignalRefs.begin(), second_it = secondSignalRefs.begin();
             first_it < firstSignalRefs.end() && second_it < secondSignalRefs.end();
             ++first_it, ++second_it)
        {
            channelSignalPairs.emplace_back(*first_it, *second_it);
        }
        return channelSignalPairs;
    }

    auto Paddle::GetTrigTime() const -> double
    {
        return std::min(fRightChannel->GetTrigTime(), fLeftChannel->GetTrigTime());
    }

} // namespace R3B::Digitizing
