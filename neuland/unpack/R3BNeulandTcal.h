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
// -----                  R3BNeuLandTcal                  -----
// -----            Created 27-01-2015 by M.Heil          -----
// ------------------------------------------------------------

#ifndef R3BNEULANDTCAL
#define R3BNEULANDTCAL

#include "FairTask.h"

class TClonesArray;
class R3BTCalModulePar;
class R3BTCalPar;
class R3BEventHeader;

/**
 * An analysis task to apply TCAL calibration for NeuLAND.
 * This class reads NeuLAND mapped items with TDC values and
 * produces time items with time in [ns]. It requires TCAL
 * calibration parameters, which are produced in a separate
 * analysis run containing R3BNeulandTcalFill task.
 */
class R3BNeulandTcal : public FairTask
{

  public:
    /**
     * Default constructor.
     * Creates an instance of the task with default parameters.
     */
    R3BNeulandTcal();

    /**
     * Standard constructor.
     * Creates an instance of the task.
     * @param name a name of the task.
     * @param iVerbose a verbosity level.
     */
    R3BNeulandTcal(const char* name, Int_t iVerbose = 1);

    /**
     * Destructor.
     * Frees the memory used by the object.
     */
    virtual ~R3BNeulandTcal();

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

    /**
     * Method for setting the trigger value.
     * @param trigger 1 - onspill, 2 - offspill, -1 - all events.
     */
    inline void SetTrigger(Int_t trigger) { fTrigger = trigger; }

    /**
     * Method for setting the number of NeuLAND modules.
     * @param nPMTs a number of photomultipliers.
     * @param n17 a number of channels with stop signal (17-th channel).
     */
    inline void SetNofModules(Int_t nPlanes, Int_t nBars, Int_t n17)
    {
        fNofPlanes = nPlanes;
        fNofBars = nBars;
        fNof17 = n17;
    }

  private:
    Int_t fNEvents;           /**< Event counter. */
    R3BEventHeader* header;   /**< Event header. */
    TClonesArray* fMappedHit; /**< Array with raw items - input data. */
    TClonesArray* fPmt;       /**< Array with time items - output data. */
    Int_t fNPmt;              /**< Number of produced time items per event. */
    R3BTCalPar* fTcalPar;     /**< TCAL parameter container. */
    Int_t fTrigger;           /**< Trigger value. */
    Int_t fNofPlanes;         /**< Number of planes. */
    Int_t fNofBars;           /**< Number of bars per plane. */
    Int_t fNof17;             /**< Number of channels with stop signal. */
    //    std::map<Int_t, Bool_t> fMap17Seen;         /**< Map with flag of observed stop signal. */
    //    std::map<Int_t, Double_t> fMapStopTime;     /**< Map with value of stop time. */
    //    std::map<Int_t, Int_t> fMapStopClock;       /**< Map with value of stop clock. */
    Double_t fClockFreq; /**< Clock cycle in [ns]. */

  public:
    ClassDef(R3BNeulandTcal, 1)
};

#endif
