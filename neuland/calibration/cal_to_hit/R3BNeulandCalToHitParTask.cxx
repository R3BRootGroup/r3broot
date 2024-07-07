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

#include "R3BNeulandCalToHitParTask.h"
#include <R3BLogger.h>
#include <R3BNeulandLSQREngineAdaptor.h>
#include <R3BNeulandMillepede.h>

namespace R3B::Neuland
{
    Cal2HitParTask::Cal2HitParTask(Cal2HitParMethod method, std::string_view name, int iVerbose)
        : CalibrationTask(name, iVerbose)
    {
        switch (method)
        {
            case Cal2HitParMethod::LSQT:
                R3BLOG(info, "Cal2HitPar method: LSQT.");
                engine_ = std::make_unique<Calibration::LSQREngineAdaptor>();
                break;
            case Cal2HitParMethod::Millipede:
                R3BLOG(info, "Cal2HitPar method: Millepede.");
                engine_ = std::make_unique<Calibration::MillepedeEngine>();
                break;
            default:
                throw R3B::logic_error("Unrecgonized calibration method!");
        }
    }

    void Cal2HitParTask::HistogramInit(DataMonitor& histograms) { engine_->HistInit(histograms); }

    void Cal2HitParTask::SetMaxModuleNum(int max_num)
    {
        max_module_num_ = max_num;
        engine_->SetMaxModuleNum(max_module_num_);
    }

    void Cal2HitParTask::ExtraInit(FairRootManager* /*rootMan*/)
    {
        cal_data_.init();

        const auto plane_num = base_par_->GetNumOfPlanes();
        if (plane_num == 0)
        {
            throw R3B::runtime_error("Plane number extracted from Map2CalPar is 0!");
        }
        engine_->SetTask(this);
        engine_->SetModuleSize(plane_num * BarsPerPlane);
        engine_->Init();
    }

    void Cal2HitParTask::SetExtraPar(FairRuntimeDb* rtdb) {}

    void Cal2HitParTask::TriggeredExec()
    {
        engine_->EventReset();
        engine_->AddSignals(cal_data_.get());
        const auto* eventHeader = GetEventHeader();
        engine_->EndOfEvent(eventHeader->GetEventno());
    }

    void Cal2HitParTask::EndOfTask()
    {
        if (hit_par_ == nullptr)
        {
            throw R3B::runtime_error("Hit parameter is nullptr during the calibration!");
        }

        R3BLOG(info, "Starting to do cosmic calibration calibration!");
        engine_->Calibrate(*hit_par_);
        engine_->EndOfTask();
    }

    auto Cal2HitParTask::CheckConditions() const -> bool
    {
        auto res = engine_->SignalFilter(cal_data_.get());
        return res;
    }
} // namespace R3B::Neuland
