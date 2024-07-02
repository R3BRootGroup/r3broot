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

#include "R3BNeulandCalToHitPar.h"
#include <R3BIOConnector.h>
#include <R3BNeulandBasePar.h>
#include <R3BNeulandCalData2.h>
#include <R3BNeulandCalibrationTask.h>
#include <R3BNeulandCosmicEngine.h>
#include <R3BNeulandMapToCalPar.h>

namespace R3B::Neuland
{
    enum class Cal2HitParMethod
    {
        LSQT,
        Millipede
    };

    class Cal2HitParTask : public CalibrationTask
    {
      public:
        explicit Cal2HitParTask(Cal2HitParMethod method = Cal2HitParMethod::LSQT,
                                std::string_view name = "NeulandCal2HitParTask",
                                int iVerbose = 1);
        void SetMinStat(int min) { engine_->SetMinStat(min); }
        void SetErrorScale(float scale) { engine_->SetErrorScale(scale); }
        auto GetCal2HitPar() -> auto* { return hit_par_; }
        void SetMaxModuleNum(int max_num);

      private:
        int max_module_num_ = MaxNumberOfBars;
        InputVectorConnector<BarCalData> cal_data_{ "NeulandCalData" };

        CalibrationBasePar* base_par_ = InputPar<CalibrationBasePar>("NeulandCalibrationBasePar");
        Map2CalPar* map_to_cal_par_ = InputPar<Map2CalPar>("LandTCalPar");
        Cal2HitPar* hit_par_ = OutputPar<Cal2HitPar>("NeulandHitPar");

        std::unique_ptr<Calibration::CosmicEngineInterface> engine_;

        // overriden functions:
        void HistogramInit(DataMonitor& histograms) override;
        void ExtraInit(FairRootManager* rootMan) override;
        void SetExtraPar(FairRuntimeDb* rtdb) override;
        void TriggeredExec() override;
        void EndOfTask() override;
        [[nodiscard]] auto CheckConditions() const -> bool override;

        // private non virtual functions:
    };

} // namespace R3B::Neuland
