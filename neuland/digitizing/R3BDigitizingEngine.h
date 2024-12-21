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

#ifndef NEULAND_DIGITIZING_ENGINE_H
#define NEULAND_DIGITIZING_ENGINE_H

#include "FairLogger.h"
#include "R3BDigitizingChannel.h"
#include "R3BDigitizingPaddle.h"
#include "Rtypes.h"
#include <cmath>
#include <functional>
#include <map>
#include <memory>

namespace R3B::Digitizing
{
    class DigitizingEngineInterface
    {
      public:
        DigitizingEngineInterface() = default;
        // rule of 5
        virtual ~DigitizingEngineInterface() = default;
        DigitizingEngineInterface(const DigitizingEngineInterface& other) = delete;
        auto operator=(const DigitizingEngineInterface& other) -> DigitizingEngineInterface& = delete;
        DigitizingEngineInterface(DigitizingEngineInterface&& other) = default;
        auto operator=(DigitizingEngineInterface&& other) -> DigitizingEngineInterface& = delete;

        virtual void DepositLight(int paddle_id, double time, double light, double dist) = 0;
        [[nodiscard]] virtual auto GetTriggerTime() const -> double = 0;
        virtual auto ExtractPaddles() -> std::map<int, std::unique_ptr<Paddle>> = 0;
        virtual void Init() = 0;
    };

    // factory classes for paddle and channel:
    template <typename ChannelClass,
              typename = typename std::enable_if<std::is_base_of<Channel, ChannelClass>::value>::type>
    struct UseChannel
    {
        template <typename... Args>
        explicit UseChannel(Args&&... args)
            : BuildChannel([&](ChannelSide side)
                           { return std::make_unique<ChannelClass>(side, std::forward<Args>(args)...); })
        {
        }
        std::function<std::unique_ptr<ChannelClass>(ChannelSide)> BuildChannel;
    };

    template <typename PaddleClass,
              typename = typename std::enable_if<std::is_base_of<Paddle, PaddleClass>::value>::type>
    struct UsePaddle
    {
        template <typename... Args>
        explicit UsePaddle(Args&&... args)
            : BuildPaddle([&](int paddleID)
                          { return std::make_unique<PaddleClass>(paddleID, std::forward<Args>(args)...); })
        {
        }
        std::function<std::unique_ptr<PaddleClass>(int)> BuildPaddle;
    };

    template <typename PaddleClass, typename ChannelClass, typename InitFunc = std::function<void()>>
    class DigitizingEngine : public DigitizingEngineInterface
    {
      private:
        UsePaddle<PaddleClass> paddleClass_;
        UseChannel<ChannelClass> channelClass_;
        std::map<int, std::unique_ptr<Paddle>> paddles;
        InitFunc initFunc_;

      public:
        DigitizingEngine(
            const UsePaddle<PaddleClass>& p_paddleClass,
            const UseChannel<ChannelClass>& p_channelClass,
            InitFunc initFunc = []() {})
            : paddleClass_{ p_paddleClass }
            , channelClass_{ p_channelClass }
            , DigitizingEngineInterface()
            , initFunc_{ initFunc }
        {
        }

        void DepositLight(int paddle_id, double time, double light, double dist) override
        {
            if (paddles.find(paddle_id) == paddles.end())
            {
                auto newPaddle = paddleClass_.BuildPaddle(paddle_id);
                newPaddle->SetChannel(channelClass_.BuildChannel(Digitizing::ChannelSide::left));
                newPaddle->SetChannel(channelClass_.BuildChannel(Digitizing::ChannelSide::right));
                paddles[paddle_id] = std::move(newPaddle);
            }
            paddles.at(paddle_id)->DepositLight({ time, light, dist });
        }

        [[nodiscard]] auto GetTriggerTime() const -> double override
        {
            auto min_element = std::min_element(paddles.begin(),
                                                paddles.end(),
                                                [](const auto& left, const auto& right)
                                                { return left.second->GetTrigTime() < right.second->GetTrigTime(); });
            return (min_element == paddles.end()) ? NAN : min_element->second->GetTrigTime();
        }

        [[nodiscard]] auto ExtractPaddles() -> std::map<int, std::unique_ptr<Paddle>> override
        {
            return std::move(paddles);
        }

        void Init() override { initFunc_(); }
        void SetInit(const InitFunc& initFunc) { initFunc_ = initFunc; }
    };

    // helper to create owning digitizingEngine:
    template <typename... Args>
    [[nodiscard]] auto CreateEngine(Args&&... args)
        -> std::unique_ptr<decltype(DigitizingEngine{ std::forward<Args>(args)... })>
    {
        return std::make_unique<decltype(DigitizingEngine{ std::forward<Args>(args)... })>(std::forward<Args>(args)...);
    }

} // namespace R3B::Digitizing

#endif // NEULAND_DIGITIZING_ENGINE_H
