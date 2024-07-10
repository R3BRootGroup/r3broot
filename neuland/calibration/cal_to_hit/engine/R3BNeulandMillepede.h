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

#include "BarPositionFilter.h"
#include "R3BNeulandCosmicEngine.h"
#include "R3BNeulandMille.h"
// #include "Utilities.h"
#include <Mille.h>
#include <ParResultReader.h>
#include <PedeLauncher.h>
#include <R3BNeulandCommon.h>
#include <optional>
// #include <RankChecker.h>

namespace R3B::Neuland::Calibration
{
    constexpr auto DEFAULT_PEDE_ERROR_THRES = 0.1F;

    class MillepedeEngine : public CosmicEngineInterface
    {
      public:
        enum class State
        {
            histogram_calibration,
            millepede_calibration
        };

        MillepedeEngine() = default;

        // Setters:
        void set_pede_interation_number(int number) { pede_interartion_number_ = number; }
        void set_pede_error_threshold(float thres) { pede_error_threshold_ = thres; }
        void set_error_factor(float scale) { error_scale_factor_ = scale; }

      private:
        State current_state_ = State::histogram_calibration;
        int minimum_hit_ = 1;
        int plane_max_hit_ = 3;
        int pede_interartion_number_ = 3;
        float pede_error_threshold_ = DEFAULT_PEDE_ERROR_THRES;
        float error_scale_factor_ = 1.F;
        double pos_residual_threshold = BarSize_XY;
        // float minimum_pos_z_ = 0;
        // float smallest_time_sum_ = 0.;
        std::optional<float> average_t_sum_;
        float init_effective_c_ = DEFAULT_EFFECTIVE_C;
        // MilleDataPoint input_data_buffer_;
        // std::string input_data_filename_ = "neuland_cosmic_mille.bin";
        std::string pede_steer_filename_ = "neuland_steer.txt";
        std::string parameter_filename_ = "neuland_pars.txt";
        // Mille binary_data_writer_{ input_data_filename_ };
        MilleHandler mille_handler_;
        Millepede::ResultReader par_result_;
        Millepede::Launcher pede_launcher_;
        // BarPositionFitter bar_positions_;

        // histograms:
        TH2I* hist_time_offsets_ = nullptr;

        // parameter:
        Cal2HitPar* cal_to_hit_par_ = nullptr;

        // private virtual methods:
        void Init() override;
        void AddSignals(const CalData& signals) override;
        void Calibrate(Cal2HitPar& hit_par) override;
        void EndOfEvent(unsigned int event_num = 0) override;
        void EventReset() override;
        auto SignalFilter(const CalData& signals) -> bool override;
        void EndOfTask() override;
        void HistInit(DataMonitor& histograms) override;
        void SetMinStat(int min) override
        {
            minimum_hit_ = min;
            R3BLOG(info, fmt::format("Minimum number of hits is set to {}", minimum_hit_));
        }

        // private non-virtual methods:
        void write_to_buffer();
        void add_signal_t_sum(const BarCalData& signal);
        void add_signal_t_diff(const BarCalData& signal);
        void add_spacial_local_constraint(int plane_num, double displacement, double error);
        void write_local_constraint(const BarPositionFitter& bar_positions);
        void set_minimum_values(const CalData& signals);
        void fill_module_parameters(const Millepede::ResultReader& result, Neuland::Cal2HitPar& cal_to_hit_par);
        void fill_data_to_mille(const CalData& signals);
        void fill_time_differences(const CalData& cal_data, TH2I* hist_time_offsets);

        template <typename FillAction>
        void fill_filtered_bar_cal_data(const CalData& signals, FillAction&& fill_action);

        void init_parameter();
        void init_steer_writer(const Cal2HitPar& cal_to_hit_par);

        void histogram_calibrate(Cal2HitPar& cal_to_hit_par);
        void millepede_calibrate(Cal2HitPar& cal_to_hit_par);

        auto calculate_time_offset_effective_speed(int module_num) -> std::pair<ValueErrorD, ValueErrorD>;

        auto get_bar_cal_residual(const BarCalData& signal, double pred) -> double;

        // inline void add_signal_t(const BarCalData& signal)
        // {
        //     // add_signal_t_sum(signal);
        //     add_signal_t_diff(signal);
        // }
    };

} // namespace R3B::Neuland::Calibration
