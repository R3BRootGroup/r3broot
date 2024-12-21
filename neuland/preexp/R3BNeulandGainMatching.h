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

#ifndef R3BNEULANDGAINMATCHING_H
#define R3BNEULANDGAINMATCHING_H

#include "FairTask.h"
#include "R3BChannelAccessEPICS.h"
#include "R3BEventHeader.h"
#include "TH1.h"
#include <TClonesArray.h>

class R3BNeulandGainMatching : public FairTask
{

  public:
    /**
     * Default constructor.
     * Creates an instance of the task with default parameters.
     */
    R3BNeulandGainMatching();

    /**
     * Standard constructor.
     * Creates an instance of the task.
     * @param name a name of the task.
     * @param iVerbose a verbosity level.
     */
    R3BNeulandGainMatching(const char* name, Int_t iVerbose = 1);

    /**
     * Destructor.
     * Frees the memory used by the object.
     */
    virtual ~R3BNeulandGainMatching();

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
     * Method for setting the update rate.
     * @param rate an update rate value (events).
     */
    inline void SetUpdateRate(Int_t rate) { fUpdateRate = rate; }

    /**
     * Method for selecting events with certain trigger value.
     * @param trigger 1 - onspill, 2 - offspill, -1 - all events.
     */
    inline void SetTrigger(Int_t trigger) { fTrigger = trigger; }

    /**
     * Method for setting number of modules in NeuLAND setup.
     * @param nPMTs a number of photomultipliers.
     */
    inline void SetNofModules(Int_t firstPlane, Int_t nPlanes, Int_t nPaddles)
    {
        fFirstPlane = firstPlane;
        fNofPlanes = nPlanes;
        fNofBarsPerPlane = nPaddles;
        fNofPMTs = nPlanes * nPaddles * 2;
    }

    /**
     * Method for setting needed statistics for data analyis.
     */
    inline void SetNeededStat(Int_t nevents) { fNEventsNeeded = nevents; }

  private:
    UInt_t fNofPlanes;       /**< Number of planes. */
    UInt_t fNofBarsPerPlane; /**< Number of bars per plane. */
    UInt_t fNofPMTs;

    UInt_t fFirstPlane;

    Int_t fUpdateRate; /**< An update rate. */
    Int_t fTrigger;    /**< Trigger value. */
    Int_t fNEvents;    /**< Event counter. */

    Int_t fNEventsNeeded;

    Bool_t finished;

    R3BChannelAccessMasterEPICS epics;

    struct
    {
        R3BChannelAccessGroup* group;
        R3BChannelAccess* vmon;
        R3BChannelAccess* vtarget;
    } ca[60][50][2];

    TH1F* hCosmicPeak[60][50][2];

    Int_t iteration[60][50][2];
    Double_t esum[60][50][2];
    Double_t ealt[60][50][2];

    Double_t hv[60][50][2];

    TClonesArray* fPmt;     /**< Array with time items - input data. */
    R3BEventHeader* header; /**< Event header - input data. */

    Int_t peakmethod;
    Double_t getcosmicpeak(TH1*);
    Double_t searchcosmicpeak(TH1*, Double_t);
    Double_t searchcosmicpeaknb(TH1*, Double_t);

  public:
    ClassDef(R3BNeulandGainMatching, 1)
};

#endif
