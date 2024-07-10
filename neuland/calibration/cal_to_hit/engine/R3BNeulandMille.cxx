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
#include "R3BNeulandMille.h"
#include "R3BNeulandCalToHitPar.h"
#include "Utilities.h"
#include <R3BNeulandCommon.h>

namespace R3B::Neuland::Calibration
{
    constexpr auto SCALE_FACTOR = 10.F;
    constexpr auto DEFAULT_ERROR_NS = 1.F;

    void MilleHandler::add_signal_t_sum(const BarCalData& signal, Cal2HitPar* cal_to_hit_par_)
    {

        buffer_clear();

        const auto module_num = static_cast<int>(signal.module_num);
        const auto pos_z = ModuleNum2ZPos<float>(static_cast<int>(module_num));

        auto init_effective_c = cal_to_hit_par_->GetModuleParAt(module_num).effectiveSpeed.value;

        const auto& left_signal = signal.left.front();
        const auto& right_signal = signal.right.front();
        const auto t_sum = (left_signal.leading_time - left_signal.trigger_time) +
                           (right_signal.leading_time - right_signal.trigger_time) - average_t_sum_;

        input_data_buffer_.measurement =
            static_cast<float>(t_sum.value / SCALE_FACTOR / 2.F - BarLength / SCALE_FACTOR / init_effective_c);
        input_data_buffer_.sigma = static_cast<float>(t_sum.error / SCALE_FACTOR / 2.);
        // input_data_buffer_.sigma = static_cast<float>(DEFAULT_MEAS_ERROR);
        const auto local_derivs_t = std::array{ 0.F, 0.F, pos_z / SCALE_FACTOR, 0.F, 0.F, 1.F };
        std::copy(local_derivs_t.begin(), local_derivs_t.end(), std::back_inserter(input_data_buffer_.locals));
        input_data_buffer_.globals.emplace_back(
            to_global_label_id(module_num, GlobalLabel::tsync, max_number_of_modules_), 1.F);
        input_data_buffer_.globals.emplace_back(
            to_global_label_id(module_num, GlobalLabel::effective_c, max_number_of_modules_),
            -BarLength / SCALE_FACTOR / 2.F / init_effective_c / init_effective_c);

        write_to_buffer();
    }

    void MilleHandler::add_signal_t_diff(const BarCalData& signal, Cal2HitPar* cal_to_hit_par)
    {
        buffer_clear();
        const auto module_num = static_cast<int>(signal.module_num);
        const auto plane_id = ModuleID2PlaneID(static_cast<int>(module_num) - 1);
        const auto is_horizontal = IsPlaneIDHorizontal(plane_id);

        const auto pos_z = static_cast<float>(PlaneID2ZPos(plane_id));

        const auto& left_signal = signal.left.front();
        const auto& right_signal = signal.right.front();
        const auto t_diff = (right_signal.leading_time - right_signal.trigger_time) -
                            (left_signal.leading_time - left_signal.trigger_time);

        const auto measurement = calculate_position_along_bar(signal, cal_to_hit_par);
        // input_data_buffer_.measurement = 4 * static_cast<float>(t_diff.value / SCALE_FACTOR);
        input_data_buffer_.measurement = static_cast<float>(measurement) / SCALE_FACTOR;
        input_data_buffer_.sigma = static_cast<float>(DEFAULT_ERROR_NS / SCALE_FACTOR);
        const auto local_derivs = is_horizontal ? std::array{ pos_z / SCALE_FACTOR, 0.F, 1.F, 0.F }
                                                : std::array{ 0.F, pos_z / SCALE_FACTOR, 0.F, 1.F };
        // const auto local_derivs = is_horizontal ? std::array{ pos_z / SCALE_FACTOR, 0.F, 0.F, 1.F, 0.F, 0.F }
        //                                         : std::array{ 0.F, pos_z / SCALE_FACTOR, 0.F, 0.F, 1.F, 0.F };
        std::copy(local_derivs.begin(), local_derivs.end(), std::back_inserter(input_data_buffer_.locals));
        input_data_buffer_.globals.emplace_back(
            to_global_label_id(module_num, GlobalLabel::offset_effective_c, max_number_of_modules_), -0.5F);
        input_data_buffer_.globals.emplace_back(
            to_global_label_id(module_num, GlobalLabel::effective_c, max_number_of_modules_),
            static_cast<float>(t_diff.value / SCALE_FACTOR / 2.));

        write_to_buffer();
        R3BLOG(
            debug,
            fmt::format(
                "Writting Mille data to binary file with meas = {} and z = {}", input_data_buffer_.measurement, pos_z));
    }

    void MilleHandler::add_spacial_local_constraint(int plane_num, double displacement, double error)
    {
        auto bar_position = BarPositionFitter::BarPosition{};
        bar_position.is_horizontal = IsPlaneIDHorizontal(plane_num - 1);
        bar_position.pos_z = PlaneID2ZPos<float>(plane_num - 1);
        bar_position.displacement = displacement;
        bar_position.displacement_error = error;
        bar_positions_.add_one_bar_signal(bar_position);
    }

    void MilleHandler::write_local_constraint(const BarPositionFitter& bar_positions)
    {
        auto fill_action = [this](const auto& bar_position)
        {
            buffer_clear();
            const auto local_derivs =
                bar_position.is_horizontal
                    ? std::array{ 0.F, static_cast<float>(bar_position.pos_z) / SCALE_FACTOR, 0.F, 1.F }
                    : std::array{ static_cast<float>(bar_position.pos_z) / SCALE_FACTOR, 0.F, 1.F, 0.F };
            input_data_buffer_.measurement = static_cast<float>(bar_position.displacement / SCALE_FACTOR);
            input_data_buffer_.sigma = static_cast<float>(bar_position.displacement_error / SCALE_FACTOR);
            std::copy(local_derivs.begin(), local_derivs.end(), std::back_inserter(input_data_buffer_.locals));
            write_to_buffer();
        };
        bar_positions.fill_to_mille_writer(fill_action);
    }
} // namespace R3B::Neuland::Calibration
