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
// -----              R3BPtofCal2Hit                      -----
// -----          Created Mar 1st 2017 by V. Wagner       -----
// -----       Convert cal data to detector hit data      -----
// ------------------------------------------------------------

#ifndef R3BPTOFCAL2HIT_H
#define R3BPTOFCAL2HIT_H

#include "FairTask.h"
#include <TArrayD.h>

#define PtofPaddlesPerPlane 8

class TClonesArray;
class R3BPtofHitPar;

class R3BPtofCal2Hit : public FairTask
{

  public:
    /**
     * Default constructor.
     * Creates an instance of the task with default parameters.
     */
    R3BPtofCal2Hit();

    /**
     * Standard constructor.
     * Creates an instance of the task.
     * @param name a name of the task.
     * @param iVerbose a verbosity level.
     */
    R3BPtofCal2Hit(const char* name, Int_t iVerbose = 1);

    /**
     * Destructor.
     * Frees the memory used by the object.
     */
    virtual ~R3BPtofCal2Hit();

    /**
     * Method for task initialization.
     * This function is called by the framework before
     * the event loop.
     * @return Initialization status. kSUCCESS, kERROR or kFATAL.
     */
    virtual InitStatus Init();

    /**
     * Method for initialization of the parameter containers.
     * Called by the framework prior to Init() method.
     */
    virtual void SetParContainers();

    /**
     * Method for re-initialization of parameter containers
     * in case the Run ID has changed.
     */
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

  private:
    void SetParameter();

    TClonesArray* fCalItems; /**< Array with CAL items - input data. */
    TClonesArray* fHitItems; /**< Array with HIT items - output data. */

    R3BPtofHitPar* fHitPar; /**< HIT parameter container. */

    Int_t fNEvents; /**< Event counter. */

    TArrayD fTOffset1; /**< an array with timeoffset parameters of modules */
    TArrayD fTOffset2; /**< an array with timeoffset parameters of modules */
    TArrayD fZScale;   /**< an array with zscale parameters of all modules */
    TArrayD fVEff;     /**< an array with veff parameters of all modules */

  public:
    ClassDef(R3BPtofCal2Hit, 1)
};

#endif
