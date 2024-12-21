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

// ------------------------------------------------------------
// -----               R3BOnlineSpectraSfib               -----
// -----            Created 13-04-2021 by M.Heil          -----
// -----               Fill online histograms             -----
// ------------------------------------------------------------

#ifndef R3BONLINESPECTRASFIB
#define R3BONLINESPECTRASFIB

#include "FairTask.h"
#include <array>
#include <fstream>
#include <iostream>
#include <sstream>
#include <vector>

#include "TClonesArray.h"
#include "TF1.h"
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
class R3BOnlineSpectraSfib : public FairTask
{

  public:
    /**
     * Default constructor.
     * Creates an instance of the task with default parameters.
     */
    R3BOnlineSpectraSfib();

    /**
     * Standard constructor.
     * Creates an instance of the task.
     * @param name a name of the task.
     * @param iVerbose a verbosity level.
     */
    R3BOnlineSpectraSfib(const char* name, Int_t iVerbose = 1);

    /**
     * Destructor.
     * Frees the memory used by the object.
     */
    virtual ~R3BOnlineSpectraSfib();

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
     * Method for setting the trigger value.
     * @param trigger 1 - physics, 2 - offspill, -1 - all events.
     */
    inline void SetTrigger(Int_t trigger) { fTrigger = trigger; }
    inline void SetTpat(Int_t tpat) { fTpat = tpat; }

    void Reset_SFIB_Histo();

  private:
    std::vector<TClonesArray*> fMappedItems;
    std::vector<TClonesArray*> fCalItems;
    std::vector<TClonesArray*> fHitItems;

    enum DetectorInstances
    {
        DET_SFIB,
        DET_MAX
    };
    const char* fDetectorNames[DET_MAX + 1] = { "Sfib", NULL };

    // check for trigger should be done globablly (somewhere else)
    R3BEventHeader* header; /**< Event header. */
    Int_t fTrigger;         /**< Trigger value. */
    Int_t fTpat;
    UInt_t fNofPlanes;
    UInt_t fPaddlesPerPlane; /**< Number of paddles per plane. */

    unsigned long fNEvents = 0, fNEvents_start = 0; /**< Event counter. */

    TF1* xyline;

    TH2F* fh_sfib_Tot_top;
    TH1F* fh_sfib_Tot_top1d;
    TH2F* fh_sfib_Tot_top_g;
    TH1F* fh_sfib_Tot_top1d_g;

    TH2F* fh_sfib_Tot_bot;
    TH1F* fh_sfib_Tot_bot1d;
    TH2F* fh_sfib_Tot_bot_g;
    TH1F* fh_sfib_Tot_bot1d_g;

    TH1F* fh_sfib_channels;
    TH1F* fh_sfib_channels_top;
    TH1F* fh_sfib_channels_bot;
    TH2F* fh_sfib_channels_topvsbot;

    TH2F* fh_sfib_tot_tvb;

    TH1F* fh_sfib_channels_g;
    // TH2F* fh_sfib_tofd;
    // TH2F* fh_sfib_pspx;
    TH1F* fh_sfib_fmult;
    TH1F* fh_sfib_cmult;
    TH1F* fh_sfib_fmult_g;
    TH1F* fh_sfib_cmult_g;
    TH1F* fh_sfib_tot_ch_top1;
    TH1F* fh_sfib_tot_ch_top2;
    TH1F* fh_sfib_tot_ch_bot1;
    TH1F* fh_sfib_tot_ch_bot2;
    TH1F* fh_sfib_tot_max_top1;
    TH1F* fh_sfib_tot_max_top2;
    TH1F* fh_sfib_tot_max_bot1;
    TH1F* fh_sfib_tot_max_bot2;
    TH2F* fh_sfib_tot_ch_max_top1;
    TH2F* fh_sfib_tot_ch_max_top2;
    TH2F* fh_sfib_tot_ch_max_bot1;
    TH2F* fh_sfib_tot_ch_max_bot2;
    TH2F* fh_sfib_tot_ch_topbot1;
    TH1F* fh_sfib_multi;
    TH1F* fh_sfib_clus;
    TH2F* fh_sfib_multi_clus;
    TH2F* fh_sfib_multi_ch_topbot1;
    TH1F* fh_sfib_tot_top_multi1;
    TH1F* fh_sfib_tot_top_multi2;

  public:
    ClassDef(R3BOnlineSpectraSfib, 1)
};

#endif
