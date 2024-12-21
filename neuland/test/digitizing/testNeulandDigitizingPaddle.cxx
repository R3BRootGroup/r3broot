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

#include "MockModels.h"
#include "R3BDigitizingPaddleNeuland.h"
#include "gtest/gtest.h"
#include <vector>

/**
 * gtest for neuland digitizing modules
 * @author Yanzhao Wang
 * @since April 12, 2023
 */
namespace
{
    namespace Digitizing = R3B::Digitizing;
    using Digitizing::Channel;
    using Digitizing::ChannelSide;
    using Digitizing::PaddleHit;
    using Digitizing::Neuland::NeulandPaddle;
    using R3B::testing::Neuland::GMockChannel;
    using ::testing::Return;
    template <typename Type>
    using Pair = Digitizing::LRPair<Type>;

    class testNeulandPaddle : public ::testing::Test
    {
      protected:
        testNeulandPaddle() = default;

        void SetUp() override
        {
            auto leftChannel = std::make_unique<GMockChannel>(ChannelSide::left);
            auto rightChannel = std::make_unique<GMockChannel>(ChannelSide::right);
            channels_.left = leftChannel.get();
            channels_.right = rightChannel.get();
            EXPECT_CALL(*channels_.left, AttachToPaddle).Times(1);
            EXPECT_CALL(*channels_.right, AttachToPaddle).Times(1);
            paddle_->SetChannel(std::move(leftChannel));
            paddle_->SetChannel(std::move(rightChannel));
        }

        // void SetExpectCall()
        // {
        //     EXPECT_CALL(*channels_.left,
        //     ConstructSignals()).Times(1).WillOnce(Return(std::move(channelSignals_.left)));
        //     EXPECT_CALL(*channels_.right, ConstructSignals())
        //         .Times(1)
        //         .WillOnce(Return(std::move(channelSignals_.right)));
        // }

        template <ChannelSide side>
        void AddChannelSignal(const double& time, const double& energy)
        {
            auto signal = Channel::Signal{};
            signal.tdc = time;
            signal.qdc = energy;
            signal.qdcUnSat = energy;
            signal.side = side;
            if constexpr (side == ChannelSide::right)
            {
                channelSignals_.right.push_back(signal);
            }
            else
            {
                channelSignals_.left.push_back(signal);
            }
        }

        template <ChannelSide side>
        auto ExtractChannelSignals()
        {
            if constexpr (side == ChannelSide::right)
            {
                return std::move(channelSignals_.right);
            }
            else
            {
                return std::move(channelSignals_.left);
            }
        }
        [[nodiscard]] const auto& GetChannels() const { return channels_; }
        [[nodiscard]] auto* GetPaddle() const { return paddle_.get(); }

      private:
        std::unique_ptr<NeulandPaddle> paddle_ = std::make_unique<NeulandPaddle>(0);
        Pair<std::vector<Channel::Signal>> channelSignals_;
        Pair<GMockChannel*> channels_ = { nullptr, nullptr };
    };

    TEST_F(testNeulandPaddle, check_inputN_to_channel) // NOLINT
    {
        auto* paddle = GetPaddle();
        paddle->DepositLight(PaddleHit{ 10., 20., 0. });
        paddle->DepositLight(PaddleHit{ 10., 20., 0. });
        paddle->DepositLight(PaddleHit{ 20., 0., 0. });
        const auto& leftHits = GetChannels().left->hits_;
        const auto& rightHits = GetChannels().right->hits_;
        ASSERT_EQ(leftHits.size(), rightHits.size()) << "Numbers of hits of both PMT channels don't match!";
        ASSERT_EQ(leftHits.size(), 3) << "Number of hits for a PMT channel doesn't match with paddle inputs";
    }

    TEST_F(testNeulandPaddle, check_light_attenuation) // NOLINT
    {
        auto* paddle = GetPaddle();
        auto const inputE = 20.;
        auto paddleHit = PaddleHit{ 10., inputE, 0 };
        paddle->DepositLight(paddleHit);
        const auto& leftHits = GetChannels().left->hits_;
        const auto& rightHits = GetChannels().right->hits_;
        ASSERT_LT(leftHits[0].light, inputE) << "No light attenuation for the left channel!";
        ASSERT_LT(rightHits[0].light, inputE) << "No light attenuation for the right channel!";
    }

    TEST_F(testNeulandPaddle, check_positive_direction) // NOLINT
    {
        auto* paddle = GetPaddle();
        auto const inputE = 20.;
        auto const initPos = 0.8 * NeulandPaddle::gHalfLength;
        auto paddleHit = PaddleHit{ 10., inputE, initPos };
        paddle->DepositLight(paddleHit);
        const auto& leftHits = GetChannels().left->hits_;
        const auto& rightHits = GetChannels().right->hits_;
        ASSERT_LT(leftHits[0].light, rightHits[0].light) << "left channel has larger light input with longer distance!";
        ASSERT_GT(leftHits[0].time, rightHits[0].time) << "light hit left channel ealier with longer distance!";
    }

    TEST_F(testNeulandPaddle, check_coupling_counts) // NOLINT
    {
        auto* paddle = GetPaddle();
        paddle->DepositLight(PaddleHit{ 10., 20., 0.5 * NeulandPaddle::gHalfLength });
        paddle->DepositLight(PaddleHit{ 60., 10., -0.3 * NeulandPaddle::gHalfLength });
        const auto& leftHits = GetChannels().left->hits_;
        const auto& rightHits = GetChannels().right->hits_;
        AddChannelSignal<ChannelSide::left>(leftHits[0].time, leftHits[0].light);
        AddChannelSignal<ChannelSide::left>(leftHits[1].time, leftHits[1].light);
        AddChannelSignal<ChannelSide::right>(rightHits[1].time, rightHits[1].light);
        AddChannelSignal<ChannelSide::right>(rightHits[0].time, rightHits[0].light);
        EXPECT_CALL(*GetChannels().left, ConstructSignals())
            .Times(1)
            .WillOnce(Return(ExtractChannelSignals<ChannelSide::left>()));
        EXPECT_CALL(*GetChannels().right, ConstructSignals())
            .Times(1)
            .WillOnce(Return(ExtractChannelSignals<ChannelSide::right>()));
        const auto& signals = paddle->GetSignals();
        ASSERT_EQ(signals.size(), 2) << "Fail to reconstruct the same number of signals!";
    }

    TEST_F(testNeulandPaddle, check_output_singal) // NOLINT
    {
        auto* paddle = GetPaddle();
        paddle->DepositLight(PaddleHit{ 10., 20., 0.5 * NeulandPaddle::gHalfLength });
        const auto& leftHits = GetChannels().left->hits_;
        const auto& rightHits = GetChannels().right->hits_;
        AddChannelSignal<ChannelSide::left>(leftHits[0].time, leftHits[0].light);
        AddChannelSignal<ChannelSide::right>(rightHits[0].time, rightHits[0].light);
        EXPECT_CALL(*GetChannels().left, ConstructSignals())
            .Times(1)
            .WillOnce(Return(ExtractChannelSignals<ChannelSide::left>()));
        EXPECT_CALL(*GetChannels().right, ConstructSignals())
            .Times(1)
            .WillOnce(Return(ExtractChannelSignals<ChannelSide::right>()));
        auto const& signals = paddle->GetSignals();
        ASSERT_EQ(signals.size(), 1) << "Fail to reconstruct the same number of signals!";
        const double err = 0.1 * 0.5 * NeulandPaddle::gHalfLength;
        ASSERT_NEAR(0.5 * NeulandPaddle::gHalfLength, signals.front().position, err)
            << "reconstructed position is way different to predefined position!";
    }

    // TEST_F(testNeulandPaddle, check_coupling_NoEQ_counts) // NOLINT
    // {
    //     auto* paddle = GetPaddle();
    //     paddle->DepositLight(PaddleHit{ 10., 20., 0.5 * NeulandPaddle::gHalfLength });
    //     paddle->DepositLight(PaddleHit{ 60., 10., -0.3 * NeulandPaddle::gHalfLength });
    //     const auto& leftHits = GetChannels().left->hits_;
    //     const auto& rightHits = GetChannels().right->hits_;
    //     AddChannelSignal<ChannelSide::left>(leftHits[0].time, leftHits[0].light);
    //     AddChannelSignal<ChannelSide::left>(leftHits[1].time, leftHits[1].light);
    //     AddChannelSignal<ChannelSide::right>(rightHits[1].time, rightHits[1].light);
    //     SetExpectCall();
    //     const auto& signals = paddle->GetSignals();
    //     ASSERT_EQ(signals.size(), 1) << "inequal number of channel signals are not paired correctly!";
    // }

    // TEST_F(testNeulandPaddle, check_coupling_badParing) // NOLINT
    // {
    //     FairLogger::GetLogger()->SetLogScreenLevel("error");
    //     auto* paddle = GetPaddle();
    //     paddle->DepositLight(PaddleHit{ 10., 20., 0.5 * NeulandPaddle::gHalfLength });
    //     paddle->DepositLight(PaddleHit{ 60., 10., -0.3 * NeulandPaddle::gHalfLength });
    //     const auto& leftHits = GetChannels().left->hits_;
    //     const auto& rightHits = GetChannels().right->hits_;
    //     AddChannelSignal<ChannelSide::left>(leftHits[0].time, leftHits[0].light);
    //     AddChannelSignal<ChannelSide::left>(leftHits[0].time, leftHits[0].light);
    //     AddChannelSignal<ChannelSide::right>(rightHits[0].time, rightHits[0].light);
    //     AddChannelSignal<ChannelSide::right>(rightHits[1].time, rightHits[1].light);
    //     SetExpectCall();
    //     const auto& signals = paddle->GetSignals();
    //     ASSERT_EQ(signals.size(), 1) << "recursive paring is not prohibited!";
    // }
} // namespace
