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

#include "R3BNeulandCalData2.h"
#include "R3BNeulandMapToCalPar.h"
#include "R3BPaddleTamexMappedData2.h"
#include <R3BIOConnector.h>
#include <R3BNeulandBasePar.h>
#include <R3BNeulandCalData.h>
#include <R3BNeulandCalibrationTask.h>
#include <R3BNeulandCommon.h>

class R3BEventHeader;
namespace R3B::Neuland
{
    class Map2CalTask : public CalibrationTask
    {
      public:
        Map2CalTask();
        Map2CalTask(std::string_view name, int iVerbose);
        void SetPulserMode(bool pulser_mode = true) { is_pulse_mode_ = pulser_mode; }
        void SetNhitmin(unsigned int size) { signal_min_size_ = size; }
        void EnableWak(bool is_walk_enabled = true) { is_walk_enabled_ = is_walk_enabled; }

      private:
        bool is_pulse_mode_ = false;
        bool is_walk_enabled_ = true;
        bool is_trig_enabled_ = true;
        float coarse_time_frequency_ = 0.; // MHz
        unsigned int coarse_time_max_num_ = MAXCTValue;
        double max_coarse_time_ = R3B::Neuland::MaxCalTime;
        unsigned int total_pmt_nums_ = 0;
        unsigned int signal_min_size_ = 1;
        unsigned int plane_num_ = 0;

        // IO data and paramters:
        InputVectorConnector<PaddleTamexMappedData> mappedData_{ "NeulandMappedData" };
        InputMapConnector<unsigned int, PaddleTamexTrigMappedData> trigMappedData_{ "NeulandTrigMappedData" };
        OutputConnector<CalData> calData_{ "NeulandCalData" };

        CalibrationBasePar* base_par_ = InputPar<CalibrationBasePar>("NeulandCalibrationBasePar");
        Map2CalPar* calibrationPar_ = InputPar<Map2CalPar>("LandTCalPar");
        Map2CalPar* calibrationTrigPar_ = InputPar<Map2CalPar>("LandTrigTCalPar");

        void ExtraInit(FairRootManager* rootMan) override;
        void HistogramInit(DataMonitor& histograms) override;
        void BeginOfEvent() override { calData_.clear(); };
        void TriggeredExec() override;
        void FinishEvent() override;
        [[nodiscard]] auto CheckConditions() const -> bool override;

        void set_pmt_num();
        void set_ct_freq();
        void sample_parameters();
        void calibrate();
        void histogram_monitor(const BarCalData& cal, Side side);
        void fill_cal_data(BarCalData& cal, const MapBarSignal& signals);
        [[nodiscard]] auto doubleEdgeSignal_to_calSignal(const DoubleEdgeSignal& double_edge_signal,
                                                         R3B::Side side,
                                                         unsigned int module_num) const -> CalDataSignal;
        [[nodiscard]] auto mapBarSignal_to_calSignals(const MapBarSignal& map_bar_signals,
                                                      unsigned int module_num,
                                                      R3B::Side side) const -> std::vector<CalDataSignal>;
        [[nodiscard]] auto convert_to_real_time(R3BTCalPar2* calPar,
                                                SingleEdgeSignal signal,
                                                FTType ftType,
                                                unsigned int module_num) const -> ValueError<double>;
        [[nodiscard]] auto get_tot(DoubleEdgeSignal pmtSignal, unsigned int module_num, Side side) const
            -> ValueError<double>;
        [[nodiscard]] auto get_trigger_time(unsigned int module_num, Side side) const -> ValueError<double>;
        void overflow_correct(R3B::Neuland::CalDataSignal& calSignal) const;
    };
} // namespace R3B::Neuland

using R3BNeulandMapped2Cal2 = R3B::Neuland::Map2CalTask;
