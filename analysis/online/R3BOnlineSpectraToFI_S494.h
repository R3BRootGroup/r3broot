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

#ifndef R3BONLINESPECTRATOFI_S494
#define R3BONLINESPECTRATOFI_S494
#define N_PLANE_MAX_TOFI 1
#define N_PADDLE_MAX_TOFI 24

#include "FairTask.h"
#include <array>
#include <fstream>
#include <iostream>
#include <sstream>
#include <vector>

#include "TClonesArray.h"
#include "TMath.h"
#include <cstdlib>

class TClonesArray;
class TH1F;
class TH2F;
class R3BEventHeader;

/**
 * This taks reads all detector data items and plots histograms
 * for online checks.
 */
class R3BOnlineSpectraToFI_S494 : public FairTask
{

  public:
    /**
     * Default constructor.
     * Creates an instance of the task with default parameters.
     */
    R3BOnlineSpectraToFI_S494();

    /**
     * Standard constructor.
     * Creates an instance of the task.
     * @param name a name of the task.
     * @param iVerbose a verbosity level.
     */
    R3BOnlineSpectraToFI_S494(const char* name, Int_t iVerbose = 1);

    /**
     * Destructor.
     * Frees the memory used by the object.
     */
    virtual ~R3BOnlineSpectraToFI_S494();

    /**
     * Method for task initialization.
     * This function is called by the framework before
     * the event loop.
     * @return Initialization status. kSUCCESS, kERROR or kFATAL.
     */
    virtual InitStatus Init();

    /**
     * Method for event loop implementation.
     * Is called by the framework every time a new event is read.
     * @param option an execution option.
     */
    virtual void Exec(Option_t* option);

    /**
     * A method for finish of processing of an event.
     * Is called by the framework for each event after executing
     * the tasks.
     */
    virtual void FinishEvent();

    /**
     * Method for finish of the task execution.
     * Is called by the framework after processing the event loop.
     */
    virtual void FinishTask();

    /**
     * Methods for setting position offset and effective velocity of light
     */

    /**
     * Method for setting the trigger value.
     * @param trigger 1 - physics, 2 - offspill, -1 - all events.
     */
    inline void SetTrigger(Int_t trigger) { fTrigger = trigger; }
    inline void SetTpat(Int_t tpat1, Int_t tpat2)
    {
        fTpat1 = tpat1;
        fTpat2 = tpat2;
    }

    /**
     * Methods for setting number of planes and paddles
     */
    inline void SetNofModules(Int_t planes, Int_t ppp)
    {
        fNofPlanes = planes;
        fPaddlesPerPlane = ppp;
    }

    void Reset_TOFI_Histo();

  private:
    std::vector<TClonesArray*> fMappedItems;
    std::vector<TClonesArray*> fCalItems;
    std::vector<TClonesArray*> fHitItems;

    TClonesArray* fCalTriggerItems; /**< Array with trigger Cal items - input data. */

    enum DetectorInstances
    {
        DET_TOFI,
        DET_MAX
    };

    const char* fDetectorNames[DET_MAX + 1] = { "Tofi", NULL };

    // check for trigger should be done globablly (somewhere else)
    R3BEventHeader* header; /**< Event header. */
    Int_t fTrigger;         /**< Trigger value. */
    Int_t fTpat1, fTpat2;
    Double_t fClockFreq; /**< Clock cycle in [ns]. */
    UInt_t fNofPlanes;
    UInt_t fPaddlesPerPlane; /**< Number of paddles per plane. */

    unsigned long long time_spill_start = 0, time_spill_end = 0, time = 0, time_prev = -1, time_prev1 = -1,
                       time_start = 0, time_mem = 0;
    unsigned long long time_previous_event = 0;

    unsigned long fNEvents = 0, fNEvents_start = 0; /**< Event counter. */

    TH1F* fh_num_bars;
    TH1F* fh_tofi_channels[N_PLANE_MAX_TOFI];
    TH2F* fh_tofi_multihit[N_PLANE_MAX_TOFI];
    TH2F* fh_tofi_TotPm[N_PLANE_MAX_TOFI];
    TH2F* fh_tofi_timePm[N_PLANE_MAX_TOFI];
    TH2F* fh_tofi_dt;
    TH2F* fh_tofi_time;
    TH1F* fh_TimePreviousEvent;
    TH2F* fh_tofi_multihit_coinc[N_PLANE_MAX_TOFI];
    TH2F* fh_tofi_TotPm_coinc[N_PLANE_MAX_TOFI];
    TH2F* fh_test;
    TH2F* fh_test1;
    TH2F* fh_test2;
    TH2F* fh_test3;
    TH2F* fh_test4;
    TH2F* fh_num_side;
    TH2F* fh_dt_hits_ToT_top;
    TH2F* fh_dt_hits_ToT_bot;
    TH2F* fh_tofi_Tot_hit[N_PLANE_MAX_TOFI];
    TH2F* fh_tofi_time_hit[N_PLANE_MAX_TOFI];
    TH1F* fh_tofi_multihit_hit[N_PLANE_MAX_TOFI];
    TH1F* fh_tofi_bars[N_PLANE_MAX_TOFI];
    TH2F* fh_tofi_dt_hit[N_PLANE_MAX_TOFI];
    TH1F* fh_tofi_Z[N_PLANE_MAX_TOFI];

  public:
    ClassDef(R3BOnlineSpectraToFI_S494, 2)
};

#endif
