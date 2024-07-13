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
#include "R3BNeulandCalData2.h"
#include <Mille.h>
#include <string>

namespace R3B::Neuland
{
    class Cal2HitPar;
    namespace Calibration
    {
        class MilleHandler
        {
          public:
            MilleHandler() = default;
            inline void init() { binary_data_writer_.set_buffer_size(MILLE_BUFFER_SIZE); }
            inline void close() { binary_data_writer_.close(); }

            inline void do_end_of_event()
            {
                if (not bar_positions_.is_skippable())
                {
                    write_local_constraint(bar_positions_);
                }
                binary_data_writer_.end();
                average_t_sum_ = 0.F;
                bar_positions_.clear();
            }

            inline void write_to_buffer() { binary_data_writer_.mille(input_data_buffer_); }
            void buffer_clear()
            {
                input_data_buffer_.locals.clear();
                input_data_buffer_.globals.clear();
                input_data_buffer_.measurement = 0.F;
                input_data_buffer_.sigma = 0.F;
            }

            inline auto fit_local() -> bool { return bar_positions_.fit_positions(); }
            inline void add_signal_t(const BarCalData& signal, Cal2HitPar* cal_to_hit_par)
            {
                // add_signal_t_sum(signal, cal_to_hit_par);
                add_signal_t_diff(signal, cal_to_hit_par);
            }

            void add_spacial_local_constraint(int plane_num, double displacement, double error);

            // Setters:
            void set_max_number_of_modules(int num_of_modules) { max_number_of_modules_ = num_of_modules; }
            void set_avarage_t_sum(float t_sum) { average_t_sum_ = t_sum; }
            void set_scale_factor(float scale_factor) { scale_factor_ = scale_factor; }

            // Getters:
            auto get_filename() const -> const auto& { return input_data_filename_; }
            auto get_bar_positions() const -> const auto& { return bar_positions_; }

          private:
            int max_number_of_modules_ = 0;
            float average_t_sum_ = 0.F;
            float scale_factor_ = 1.F;
            static constexpr auto MILLE_BUFFER_SIZE = std::size_t{ 100000 };
            std::string input_data_filename_ = "neuland_cosmic_mille.bin";
            Mille binary_data_writer_{ input_data_filename_ };
            MilleDataPoint input_data_buffer_;
            BarPositionFitter bar_positions_;

            void write_local_constraint(const BarPositionFitter& bar_positions);
            void add_signal_t_sum(const BarCalData& signal, Cal2HitPar* cal_to_hit_par_);
            void add_signal_t_diff(const BarCalData& signal, Cal2HitPar* cal_to_hit_par);
        };

    } // namespace Calibration
} // namespace R3B::Neuland
