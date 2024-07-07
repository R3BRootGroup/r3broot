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
#include <Fit/Fitter.h>
#include <Mille.h>
#include <ParResultReader.h>
#include <PedeLauncher.h>
#include <R3BNeulandCommon.h>
#include <TF1.h>
#include <optional>
// #include <RankChecker.h>

class TGraphErrors;

namespace R3B::Neuland::Calibration
{
    enum class GlobalLabel
    {
        tsync,              // tsync
        offset_effective_c, // offset times effective_C
        effective_c         // effective speed of light
    };

    class BarPositionFitter
    {
      public:
        struct BarPosition
        {
            bool is_horizontal = true;
            double displacement = 0.;
            double pos_z = 0.;
            double displacement_error = 0.;
        };

        BarPositionFitter();
        void clear()
        {
            horizontal_bars_.clear();
            vertical_bars_.clear();
            y_pos_fitting_function_.SetParameter(0, 0.);
            y_pos_fitting_function_.SetParameter(1, 0.);
            x_pos_fitting_function_.SetParameter(0, 0.);
            x_pos_fitting_function_.SetParameter(1, 0.);
        }

        void add_one_bar_signal(const BarPosition& bar_position)
        {
            auto& bar_positions = bar_position.is_horizontal ? horizontal_bars_ : vertical_bars_;
            bar_positions.displacements.push_back(bar_position.displacement);
            bar_positions.displacement_errors.push_back(bar_position.displacement_error);
            bar_positions.pos_z.push_back(bar_position.pos_z);
        }

        auto fit_positions() -> bool;
        auto get_prediction(bool is_horizontal, double pos_z) -> double;

        template <typename FillAction>
        void fill_to_mille_writer(FillAction&& do_action) const
        {
            // TODO: reformat duplications here
            namespace rng = ranges;
            auto bar_position = BarPosition{};
            for (const auto& [displacement, error, pos_z] : rng::views::zip(
                     horizontal_bars_.displacements, horizontal_bars_.displacement_errors, horizontal_bars_.pos_z))
            {
                bar_position.is_horizontal = true;
                bar_position.pos_z = pos_z;
                bar_position.displacement = displacement;
                bar_position.displacement_error = error;
                do_action(bar_position);
            }
            for (const auto& [displacement, error, pos_z] : rng::views::zip(
                     vertical_bars_.displacements, vertical_bars_.displacement_errors, vertical_bars_.pos_z))
            {
                bar_position.is_horizontal = false;
                bar_position.pos_z = pos_z;
                bar_position.displacement = displacement;
                bar_position.displacement_error = error;
                do_action(bar_position);
            }
        }

      private:
        static constexpr auto fitting_function_range = 200.; // cm
        struct BarPositions
        {
            std::vector<double> displacements;
            std::vector<double> displacement_errors;
            std::vector<double> pos_z;
            void clear()
            {
                displacements.clear();
                displacement_errors.clear();
                pos_z.clear();
            }
            [[nodiscard]] auto size() const -> unsigned int { return displacements.size(); }
        };

        ROOT::Fit::Fitter fitter_;
        TF1 y_pos_fitting_function_{ "hfun", "pol1", -fitting_function_range, fitting_function_range };
        TF1 x_pos_fitting_function_{ "vfun", "pol1", -fitting_function_range, fitting_function_range };
        BarPositions horizontal_bars_;
        BarPositions vertical_bars_;

        auto fit_with(const BarPositions& positions, TF1& func) -> bool;
    };

    class MillepedeEngine : public CosmicEngineInterface
    {
      public:
        enum class State
        {
            histogram_calibration,
            millepede_calibration
        };

        MillepedeEngine() = default;
        void enable_rank_check(bool rank_check = true) { has_rank_check_ = rank_check; }

      private:
        bool has_rank_check_ = false;
        State current_state_ = State::histogram_calibration;
        int minimum_hit_ = 1;
        int plane_max_hit_ = 3;
        double pos_residual_threshold = BarSize_XY;
        float error_scale_factor_ = 1000.F;
        // float minimum_pos_z_ = 0;
        // float smallest_time_sum_ = 0.;
        std::optional<float> average_t_sum_;
        float init_effective_c_ = DEFAULT_EFFECTIVE_C;
        MilleDataPoint input_data_buffer_;
        std::string input_data_filename_ = "neuland_cosmic_mille.bin";
        std::string pede_steer_filename_ = "neuland_steer.txt";
        std::string parameter_filename_ = "neuland_pars.txt";
        Mille binary_data_writer_{ input_data_filename_ };
        Millepede::ResultReader par_result_;
        Millepede::Launcher pede_launcher_;
        BarPositionFitter bar_positions_;

        // histograms:
        TGraphErrors* graph_time_offset_ = nullptr;
        TGraphErrors* graph_time_sync_ = nullptr;
        TGraphErrors* graph_effective_c_ = nullptr;
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
        void SetErrorScale(float scale) override { error_scale_factor_ = scale; }

        // private non-virtual methods:
        void buffer_clear();
        void write_to_buffer();
        void add_signal_t_sum(const BarCalData& signal);
        void add_signal_t_diff(const BarCalData& signal);
        void add_spacial_local_constraint(int plane_num, double displacement, double error);
        void write_local_constraint(const BarPositionFitter& bar_positions);
        void set_minimum_values(const CalData& signals);
        inline auto to_global_label_id(int module_num, GlobalLabel label) -> int;
        inline auto to_module_num_label(int par_num) -> std::pair<int, GlobalLabel>;
        void fill_module_parameters(const Millepede::ResultReader& result, Neuland::Cal2HitPar& cal_to_hit_par);
        void fill_data_to_figure(Cal2HitPar& cal_to_hit_par);
        void fill_data_to_mille(const CalData& signals);
        void fill_time_differences(const CalData& cal_data, TH2I* hist_time_offsets);

        template <typename FillAction>
        void fill_filtered_bar_cal_data(const CalData& signals, FillAction&& fill_action);

        void init_parameter();
        void init_steer_writer(const Cal2HitPar& cal_to_hit_par);

        void histogram_calibrate(Cal2HitPar& cal_to_hit_par);
        void millepede_calibrate(Cal2HitPar& cal_to_hit_par);

        auto calculate_time_offset_effective_speed(int module_num) -> std::pair<ValueErrorD, ValueErrorD>;
        auto calculate_position(const BarCalData& signal) -> double;

        auto get_bar_cal_residual(const BarCalData& signal, double pred) -> double;

        inline void add_signal_t(const BarCalData& signal)
        {
            // add_signal_t_sum(signal);
            add_signal_t_diff(signal);
        }
    };

} // namespace R3B::Neuland::Calibration
