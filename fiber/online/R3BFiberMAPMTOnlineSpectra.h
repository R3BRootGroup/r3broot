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

// -------------------------------------------------------------
// -----      R3BFiberMAPMTOnlineSpectra source file       -----
// -----    Created 25/04/22 by J.L. Rodriguez-Sanchez     -----
// -------------------------------------------------------------

#pragma once

#include <FairTask.h>

#include <TMath.h>
#include <array>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <list>
#include <sstream>
#include <vector>

class TClonesArray;
class TH1F;
class TH2F;
class R3BEventHeader;
class R3BFiberMappingPar;
class R3BFiberMAPMTCalData;

class R3BFiberMAPMTOnlineSpectra : public FairTask
{
  public:
    struct Channel
    {
        std::list<R3BFiberMAPMTCalData const*> lead_list;
    };
    /**
     * Default constructor.
     * Creates an instance of the task with default parameters.
     */
    R3BFiberMAPMTOnlineSpectra();

    /**
     * Standard constructor.
     * Creates an instance of the task.
     * @param name a name of the task.
     * @param iVerbose a verbosity level.
     */
    R3BFiberMAPMTOnlineSpectra(const TString, Int_t iVerbose = 1);

    /**
     * Destructor.
     * Frees the memory used by the object.
     */
    virtual ~R3BFiberMAPMTOnlineSpectra();

    /**
     * Method for task initialization.
     * This function is called by the framework before
     * the event loop.
     * @return Initialization status. kSUCCESS, kERROR or kFATAL.
     */
    virtual InitStatus Init();

    virtual InitStatus ReInit();

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

    virtual void SetParContainers();

    /**
     * Method for setting the trigger value.
     * @param trigger 1 - physics, 2 - offspill, -1 - all events.
     */
    void SetTrigger(Int_t trigger) { fTrigger = trigger; }
    void SetTpat(Int_t tpat1, Int_t tpat2)
    {
        fTpat1 = tpat1;
        fTpat2 = tpat2;
    }
    void SetClock(Double_t CF, Double_t CP)
    {
        fClockFreq = CF;
        fClockPeriods = CP;
    }

    void Reset_Histo();

  private:
    TString fName;
    UInt_t fNbfibersplot;
    UInt_t fNbfibers;
    R3BFiberMappingPar* fMapPar;
    TClonesArray* fMappedItems;
    TClonesArray* fCalItems;
    TClonesArray* fHitItems;
    TClonesArray* fCalTriggerItems;

    // check for trigger should be done globablly (somewhere else)
    R3BEventHeader* header; /**< Event header. */
    Int_t fTrigger;         /**< Trigger value. */
    Int_t fTpat1, fTpat2;
    Double_t fClockFreq; /**< Clock cycle in [ns]. */
    Double_t fClockPeriods;
    std::vector<Channel> fChannelArray[2];
    unsigned const* fTriggerMap[2];
    unsigned long fNEvents; /**< Event counter. */

    TH1F* fh_channels_Fib;
    TH1F* fh_channels_single_Fib;
    TH1F* fh_fibers_Fib;
    TH1F* fh_mult_Fib;
    TH2F* fh_Fib_pos;
    TH2F* fh_time_Fib;
    TH2F* fh_multihit_m_Fib;
    TH2F* fh_multihit_s_Fib;
    TH2F* fh_ToT_Fib;
    TH2F* fh_Fib_vs_Events;
    TH2F* fh_ToTup_vs_ToTdown;
    TH2F* fh_chan_corell;
    TH2F* fh_raw_tot_up;
    TH2F* fh_raw_tot_down;
    TH2F* fh_chan_dt_cal;

  public:
    ClassDef(R3BFiberMAPMTOnlineSpectra, 1)
};
