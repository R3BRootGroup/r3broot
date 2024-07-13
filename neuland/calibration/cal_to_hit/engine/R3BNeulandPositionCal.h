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

#include <R3BNeulandCalData2.h>
#include <R3BValueError.h>

class TH2I;

namespace R3B
{
    class DataMonitor;
    namespace Neuland
    {
        class Cal2HitPar;
    }
} // namespace R3B

namespace R3B::Neuland::Calibration
{
    class PositionCalibrator
    {
      public:
        PositionCalibrator() = default;
        void set_max_number_of_modules(int num_of_modules) { max_number_of_modules_ = num_of_modules; }
        void fill_time_differences(const CalData& cal_data);
        void HistInit(DataMonitor& histograms);
        void calibrate(Cal2HitPar& cal_to_hit_par);

      private:
        int max_number_of_modules_ = 0;
        TH2I* hist_time_offsets_ = nullptr;
        auto calculate_time_offset_effective_speed(int module_num) -> std::pair<ValueErrorD, ValueErrorD>;
    };
} // namespace R3B::Neuland::Calibration
