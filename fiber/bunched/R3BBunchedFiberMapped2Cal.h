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
// -----              R3BBunchedFiberMapped2Cal           -----
// -----            Created Jan 13th 2018 by M.Heil       -----
// ----- Convert mapped data to time calibrated data      -----
// ------------------------------------------------------------

#ifndef R3BBUNCHEDFIBERMAPPED2CAL
#define R3BBUNCHEDFIBERMAPPED2CAL 1

#include "FairTask.h"
#include "Rtypes.h"
#include <R3BTCalEngine.h>
#include <string.h>

class TClonesArray;
class R3BTCalPar;

/**
 * An analysis task to apply TCAL calibration.
 * This class reads mapped items with TDC values and
 * produces time items with time in [ns]. It requires TCAL
 * calibration parameters, which are produced in a separate
 * analysis run containing R3BBunchedFiberMapped2CalPar task.
 */
class R3BBunchedFiberMapped2Cal : public FairTask
{
  public:
    /**
     * Standard constructor.
     * Creates an instance of the task.
     * @param name a name of the task.
     * @param iVerbose a verbosity level.
     * @param a_variant CTDC firmware variant, see R3BTCalEngine.
     * @param a_skip_spmt Don't process SPMT side for pure MAPMT tests.
     */
    R3BBunchedFiberMapped2Cal(const char*, Int_t, enum R3BTCalEngine::CTDCVariant, Bool_t = false);

    enum Electronics
    {
        CTDC,
        TAMEX
    };
    /**
     * Standard constructor.
     * Creates an instance of the task.
     * @param name a name of the task.
     * @param iVerbose a verbosity level.
     * @param a_variant CTDC firmware variant, see R3BTCalEngine.
     * @param a_skip_spmt Don't process SPMT side for pure MAPMT tests.
     */
    R3BBunchedFiberMapped2Cal(const char*, Int_t, enum Electronics, enum R3BTCalEngine::CTDCVariant, Bool_t = false);

    /**
     * Destructor.
     * Frees the memory used by the object.
     */
    virtual ~R3BBunchedFiberMapped2Cal();

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

    /** Accessor to select online mode **/
    void SetOnline(Bool_t option) { fOnline = option; }

  private:
    // Selector for online data storage
    Bool_t fOnline;
    enum Electronics fSPMTElectronics;
    enum R3BTCalEngine::CTDCVariant fCTDCVariant;
    TString fName;
    UInt_t fSideNum[2];
    Bool_t fSkipSPMT;
    R3BTCalPar* fMAPMTTCalPar;
    R3BTCalPar* fMAPMTTrigTCalPar;
    R3BTCalPar* fSPMTTCalPar;
    TClonesArray* fMappedItems;
    TClonesArray* fCalItems;
    TClonesArray* fCalTriggerItems;
    Double_t fClockFreq;
    Double_t fTamexFreq;
    Int_t fnEvents;

  public:
    ClassDef(R3BBunchedFiberMapped2Cal, 4)
};

#endif
