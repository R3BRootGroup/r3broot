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

// NOTE: comments below are only meant for eduational purpose. DO NOT include them in your code!

// ROOT headers
#include <TClonesArray.h>
#include <TMath.h>
#include <TRandom.h>

// Fair headers
#include <FairLogger.h>
#include <FairRootManager.h>
#include <FairRunAna.h>
#include <FairRuntimeDb.h>

// R3B headers
#include "NewTask.h"
#include <R3BEventHeader.h>
#include <R3BLogger.h>

#include <fmt/core.h>

namespace R3B
{
    // ---- Standard Constructor ------------------------------------------
    NewTask::NewTask(const std::string& name, int iVerbose)
        : FairTask(name.c_str(), iVerbose)
    {
    }

    // ---- Default constructor -------------------------------------------
    // use constructor delegation
    NewTask::NewTask()
        : NewTask("NewTask", 1)
    {
    }

    // ----  Initialisation  ----------------------------------------------
    void NewTask::SetParContainers()
    {
        R3BLOG(debug2, "SetParContainers of NewTask");
        // Load all necessary parameter containers from the runtime data base
    }

    // ---- Init ----------------------------------------------------------
    auto NewTask::Init() -> InitStatus
    {
        R3BLOG(debug2, "NewTask::Init()");

        // Get a handle from the IO manager
        auto* ioman = FairRootManager::Instance();
        if (ioman == nullptr)
        {
            return kFATAL;
        }

        // ---------------- input initialisation ---------------
        // Case 1. Input std::vector:
        input_data_vector_.init();

        // Case 2. Input TClonesArray (NOT RECOMMENDED):
        input_data_TCA_ = dynamic_cast<TClonesArray*>(ioman->GetObject("DetectorCalData"));
        if (input_data_TCA_ == nullptr)
        {
            throw R3B::runtime_error("DetectorCalData is not found!");
        }

        // Case 3. Input a TObject:
        if (eventHeader_ = dynamic_cast<R3BEventHeader*>(ioman->GetObject("EventHeader.")); eventHeader_ == nullptr)
        {
            throw R3B::runtime_error("R3BEventHeader is nullptr!");
        }

        // ---------------- output initialisation ---------------
        // Case 1. Output std::vector
        output_data_vector_.init();

        // Case 2. Ouput TClonesArray (NOT RECOMMENDED)
        ioman->Register("DetectorHitData", "Detector", output_data_TCA_, kTRUE);

        // Case 3. Output a TObject
        output_data_ = std::make_unique<DetectorData>().release();
        ioman->Register(output_data_->GetName(), "Detector", output_data_, kTRUE);

        // Do whatever else is needed at the initilization stage
        // Create histograms to be filled
        // initialize variables

        return kSUCCESS;
    }

    // ---- Exec ----------------------------------------------------------
    void NewTask::Exec(Option_t* /*opt*/)
    {
        R3BLOG(debug, "Exec of NewTask");
        output_data_vector_.clear();
        output_data_TCA_->Clear();
        output_data_->reset();

        // --------- input data reading -------------
        // Case 1: read std::vector
        for (const auto& calData : input_data_vector_)
        {
            // do something with calData
            fmt::print("qdc value: {}\n", calData.left_time);
        }
        // Case 2: read TClonesArray:
        for (const auto* calData_ptr : TRangeDynCast<DetectorCalData>(input_data_TCA_))
        {
            // do something with calData
            fmt::print("qdc value: {}\n", (*calData_ptr).left_time);
        }
        // Case 3: read TObject:
        fmt::print("Event number: {}\n", eventHeader_->GetEventno());

        // --------- output data writing -------------
        // Case 1: write std::vector:
        auto& hit_data = output_data_vector_.get().emplace_back();
        hit_data.energy = 0.;
        hit_data.time = 1;
        hit_data.position = 0.;
        // Case 2: write TClonesArray:
        auto* hit_data_ptr =
            dynamic_cast<DetectorHitData*>(output_data_TCA_->ConstructedAt(output_data_TCA_->GetEntriesFast()));
        if (hit_data_ptr != nullptr)
        {
            (*hit_data_ptr).energy = 0.;
            (*hit_data_ptr).time = 1;
            (*hit_data_ptr).position = 0.;
        }
        // Case 3: write TObject:
        if (output_data_ != nullptr)
        {
            (*output_data_).energy = 0.;
        }
    }

    // ---- Finish --------------------------------------------------------
    void NewTask::Finish() { R3BLOG(debug2, "Finish of NewTask"); }
} // namespace R3B

ClassImp(R3B::NewTask);
