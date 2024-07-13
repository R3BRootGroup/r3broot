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
#pragma once

#include "R3BNeulandCosmicEngine.h"
#include "R3BNeulandMille.h"
#include "R3BNeulandPede.h"
#include "R3BNeulandPositionCal.h"
#include <optional>

namespace R3B::Neuland::Calibration
{
    class MillepedeEngine : public CosmicEngineInterface
    {
      public:
        enum class State
        {
            pre_calibration,
            millepede_calibration
        };

        MillepedeEngine() = default;

        // Setters:
        void set_pede_interation_number(int number) { pede_handler_.set_pede_interation_number(number); }
        void set_pede_error_threshold(float thres) { pede_handler_.set_pede_error_threshold(thres); }
        void set_error_factor(float scale) { pede_handler_.set_error_factor(scale); }

      private:
        State current_state_ = State::pre_calibration;
        int minimum_hit_ = 1;
        int plane_max_hit_ = 3;
        std::optional<float> average_t_sum_;
        double pos_residual_threshold = BarSize_XY;
        MilleHandler mille_handler_;
        PedeHandler pede_handler_;
        PositionCalibrator position_cal_;

        // parameter:
        Cal2HitPar* cal_to_hit_par_ = nullptr;

        // private virtual methods:
        void Init() override;
        void AddSignals(const CalData& signals) override;
        void Calibrate(Cal2HitPar& hit_par) override;
        void EndOfEvent(unsigned int event_num = 0) override;
        void EventReset() override {}
        auto SignalFilter(const CalData& signals) -> bool override;
        void EndOfTask() override;
        void HistInit(DataMonitor& histograms) override { position_cal_.HistInit(histograms); }
        void SetMinStat(int min) override
        {
            minimum_hit_ = min;
            R3BLOG(info, fmt::format("Minimum number of hits is set to {}", minimum_hit_));
        }

        // private non-virtual methods:
        void set_minimum_values(const CalData& signals);
        void fill_data_to_mille(const CalData& signals);
        void fill_time_differences(const CalData& cal_data, TH2I* hist_time_offsets);
        template <typename FillAction>
        void fill_filtered_bar_cal_data(const CalData& signals, FillAction&& fill_action);
        void pre_calibrate(Cal2HitPar& cal_to_hit_par);
        void millepede_calibrate(Cal2HitPar& cal_to_hit_par);
        auto get_bar_cal_residual(const BarCalData& signal, double pred) -> double;
    };
} // namespace R3B::Neuland::Calibration
