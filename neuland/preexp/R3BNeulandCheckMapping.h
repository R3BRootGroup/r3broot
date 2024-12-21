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

#ifndef R3BNEULANDCHECKMAPPING_H
#define R3BNEULANDCHECKMAPPING_H

#include "FairTask.h"
#include "R3BChannelAccessEPICS.h"
#include "R3BEventHeader.h"
#include "TH1.h"
#include "TTimeStamp.h"
#include <TClonesArray.h>

class R3BNeulandCheckMapping : public FairTask
{

  public:
    /**
     * Default constructor.
     * Creates an instance of the task with default parameters.
     */
    R3BNeulandCheckMapping();

    /**
     * Standard constructor.
     * Creates an instance of the task.
     * @param name a name of the task.
     * @param iVerbose a verbosity level.
     */
    R3BNeulandCheckMapping(const char* name, Int_t iVerbose = 1);

    /**
     * Destructor.
     * Frees the memory used by the object.
     */
    virtual ~R3BNeulandCheckMapping();

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

    inline void SetNofModules(Int_t nPlanes, Int_t nBars)
    {
        fNofPlanes = nPlanes;
        fPaddlesPerPlane = nBars;
    }

  private:
    UInt_t fNofPlanes;       /**< Number of planes. */
    UInt_t fPaddlesPerPlane; /**< Number of bars per plane. */

    Int_t fTrigger; /**< Trigger value. */

    UInt_t iPlnHvOn;
    UInt_t iBarHvOn;
    UInt_t iPmtHvOn;

    Bool_t finished;

    TTimeStamp timestamp0, timestamp1;
    R3BChannelAccessMasterEPICS epics;

    time_t timer0; // start time
    time_t timer1; // current

    struct
    {
        R3BChannelAccessGroup* group;
        R3BChannelAccess* vmon;
        R3BChannelAccess* vtarget;
    } ca[60][50][2];

    UInt_t cntOk[60][50][2];

    TClonesArray* fMapped;  /**< Array with mapped data - input data. */
    R3BEventHeader* header; /**< Event header - input data. */

    TH1F* h_countsok;

  public:
    ClassDef(R3BNeulandCheckMapping, 1)
};

#endif
