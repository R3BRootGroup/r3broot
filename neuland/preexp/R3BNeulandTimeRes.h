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

#ifndef R3BNEULANDTIMERES_H
#define R3BNEULANDTIMERES_H

#include "FairTask.h"
#include "R3BEventHeader.h"
#include "TClonesArray.h"
#include "TF1.h"
#include "TH1.h"
#include "TH2.h"

class R3BNeulandTimeRes : public FairTask
{

  public:
    /**
     * Default constructor.
     * Creates an instance of the task with default parameters.
     */
    R3BNeulandTimeRes();

    /**
     * Standard constructor.
     * Creates an instance of the task.
     * @param name a name of the task.
     * @param iVerbose a verbosity level.
     */
    R3BNeulandTimeRes(const char* name, Int_t iVerbose = 1);

    /**
     * Destructor.
     * Frees the memory used by the object.
     */
    virtual ~R3BNeulandTimeRes();

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
    inline void SetNofModules(Int_t nPlanes, Int_t nPaddles)
    {
        fNofPlanes = nPlanes;
        fNofBarsPerPlane = nPaddles;
    }

    /**
     * Method for setting needed statistics for data analyis.
     */
    inline void SetNeededStat(Int_t nevents) { fNEventsNeeded = nevents; }

  private:
    UInt_t fNofPlanes;       /**< Number of planes. */
    UInt_t fNofBarsPerPlane; /**< Number of bars per plane. */

    Int_t fUpdateRate; /**< An update rate. */
    Int_t fTrigger;    /**< Trigger value. */
    Int_t fNEvents;    /**< Event counter. */

    Int_t fNEventsNeeded;

    TF1* fitfunc[60][50];

    Int_t bars_done;
    Int_t bar_done[60][50];
    Bool_t finished;

    TH1F* hTimeRes[60][50];
    TH2F* hTimeResQ[60][50];

    TClonesArray* fPmt;     /**< Array with time items - input data. */
    R3BEventHeader* header; /**< Event header - input data. */

    Double_t wlk(Double_t x);

  public:
    ClassDef(R3BNeulandTimeRes, 1)
};

#endif
