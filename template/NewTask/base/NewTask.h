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

#pragma once

// NOTE: comments below are only meant for eduational purpose. DO NOT include them in your code!

#include "DetectorData.h"
#include <FairTask.h>
#include <R3BIOConnector.h>
#include <R3BNeulandCalData.h>
#include <R3BNeulandHit.h>
#include <string>

class R3BEventHeader;
// namespace here is optional.
namespace R3B
{
    // If R3B namespace is not used, the task should be named with R3BNewTask
    class NewTask : public FairTask
    {
      public:
        // Default constructor
        NewTask();

        // Standard constructor
        explicit NewTask(const std::string& name, int iVerbose = 1);

        // Other speical functions. Either define all these 5 functions or none of them (rule of 5).
        // Defining none of them is preferred (rule of 0).
        // ~NewTask() override;
        // NewTask(const NewTask&) = delete;
        // NewTask(NewTask&&) = delete;
        // NewTask& operator=(const NewTask&) = delete;
        // NewTask& operator=(NewTask&&) = delete;

        // Method to setup online mode
        void SetOnline(bool is_online) { is_online_ = is_online; }

      private:
        // NOTE: all member variables should be default initiliazed right here

        // Store data for online
        // Naming convenction of a boolean variable should be started with is_ or has_
        bool is_online_ = false;

        // -----------------Input data from previous already existing data level---------------
        // Case 1: Input data from std::vector container
        R3B::InputVectorConnector<DetectorCalData> input_data_vector_{ "DetectorCalData" };
        // or
        // InputConnector<std::vector<R3BNeulandCalData>> input_data_{ "DetectorHitData" };
        //
        // Case 2: Input data from TClonesArray container
        TClonesArray* input_data_TCA_ = nullptr;
        //
        // Case 3: Input a TObject
        R3BEventHeader* eventHeader_ = nullptr;

        // -----------------Output array to new data level-----------------
        // Case 1: Ouput data container in std::vector (RECOMMENDED)
        R3B::OutputVectorConnector<DetectorHitData> output_data_vector_{ "DetectorHitData" };
        // or
        // OutputConnector<std::vector<R3BNeulandHit>> output_data_{ "DetectorHitData" };
        //
        // Case 2: Output data container in TClonesArray (NOT RECOMMENDED)
        TClonesArray* output_data_TCA_ = nullptr;
        //
        // Case 3: Output a TObject
        DetectorData* output_data_ = nullptr;

        // virtual functions should be private

        // Initiliazation of task at the beginning of a run
        auto Init() -> InitStatus override;

        // Executed for each event
        void Exec(Option_t* opt) override;

        // Load the parameter container from the runtime database
        void SetParContainers() override;

        // Finish task called at the end of the run
        void Finish() override;

      public:
        // Class definition
        ClassDefOverride(NewTask, 1); // NOLINT
    };
} // namespace R3B
