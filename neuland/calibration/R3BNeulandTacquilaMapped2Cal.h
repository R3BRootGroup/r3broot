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

#ifndef R3BNEULANDTACQUILAMAPPED2CAL_H
#define R3BNEULANDTACQUILAMAPPED2CAL_H

#include "FairTask.h"

class TClonesArray;
class TH1F;
class R3BTCalModulePar;
class R3BTCalPar;
class R3BEventHeader;
class R3BNeulandQCalPar;

/**
 * An analysis task to apply TCAL calibration for NeuLAND.
 * This class reads NeuLAND raw items with TDC values and
 * produces time items with time in [ns]. It requires TCAL
 * calibration parameters, which are produced in a separate
 * analysis run containing R3BNeulandTacquilaMapped2CalPar task.
 * @author D. Kresan
 * @since September 7, 2015
 */
class R3BNeulandTacquilaMapped2Cal : public FairTask
{

  public:
    /**
     * Default constructor.
     * Creates an instance of the task with default parameters.
     */
    R3BNeulandTacquilaMapped2Cal();

    /**
     * Standard constructor.
     * Creates an instance of the task.
     * @param name a name of the task.
     * @param iVerbose a verbosity level.
     */
    R3BNeulandTacquilaMapped2Cal(const char* name, Int_t iVerbose = 1);

    /**
     * Destructor.
     * Frees the memory used by the object.
     */
    virtual ~R3BNeulandTacquilaMapped2Cal();

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
     */
    inline void SetNofModules(Int_t nPMTs) { fNofPMTs = nPMTs; }

    /**
     * Method to set running mode for pulser data analysis.
     * @param mode a boolean flag - if TRUE events with all PMT's fired will be taken.
     */
    inline void SetPulserMode(Bool_t mode = kTRUE) { fPulserMode = mode; }

    /**
     * Method to enable / disable walk corrections.
     * @param walk a boolean flag - if TRUE, walk corrections will be applied.
     */
    inline void EnableWalk(Bool_t walk = kTRUE) { fWalkEnabled = walk; }

  private:
    void SetParameter();

    Int_t fNEvents;                          /**< Event counter. */
    Bool_t fPulserMode;                      /**< Running with pulser data. */
    Bool_t fWalkEnabled;                     /**< Enable / Disable walk correction. */
    R3BEventHeader* header;                  /**< Event header. */
    TClonesArray* fRawHit;                   /**< Array with raw items - input data. */
    TClonesArray* fPmt;                      /**< Array with time items - output data. */
    Int_t fNPmt;                             /**< Number of produced time items per event. */
    R3BTCalPar* fTcalPar;                    /**< TCAL parameter container. */
    R3BNeulandQCalPar* fQCalPar;             /**< QCAL parameter container. */
    Int_t fTrigger;                          /**< Trigger value. */
    Int_t fNofPMTs;                          /**< Number of photomultipliers. */
    std::map<Int_t, Bool_t> fMap17Seen;      /**< Map with flag of observed stop signal. */
    std::map<Int_t, Double_t> fMapStopTime;  /**< Map with value of stop time. */
    std::map<Int_t, Int_t> fMapStopClock;    /**< Map with value of stop clock. */
    std::map<Int_t, Double_t> fMapQdcOffset; /**< Map with value of qdc offset. */
    Double_t fClockFreq;                     /**< Clock cycle in [ns]. */
    TH1F* fh_pulser_5_2;                     /**< Resolution of one PMT. */
    TH1F* fh_pulser_105_2;                   /**< Resolution of one PMT. */

    void MakeCal();

  public:
    ClassDef(R3BNeulandTacquilaMapped2Cal, 1)
};

#endif
