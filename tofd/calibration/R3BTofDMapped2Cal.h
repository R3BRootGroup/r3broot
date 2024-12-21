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

#pragma once

#include <vector>

#include "FairTask.h"

class TClonesArray;
class R3BTofDMappingPar;
class R3BTCalPar;
class R3BTofdMappedData;
class R3BTofdCalData;
class R3BEventHeader;

class R3BTofDMapped2Cal : public FairTask
{
  public:
    /**
     * Default constructor.
     * Creates an instance of the task with default parameters.
     */
    R3BTofDMapped2Cal();

    /**
     * Standard constructor.
     * Creates an instance of the task.
     * @param name a name of the task.
     * @param iVerbose a verbosity level.
     */
    R3BTofDMapped2Cal(const char*, Int_t = 1);

    /**
     * Destructor.
     * Frees the memory used by the object.
     */
    virtual ~R3BTofDMapped2Cal();

    /**
     * Method for task initialization.
     * This function is called by the framework before
     * the event loop.
     * @return Initialization status. kSUCCESS, kERROR or kFATAL.
     */
    InitStatus Init();

    /**
     * Method for initialization of the parameter containers.
     * Called by the framework prior to Init() method.
     */
    void SetParContainers();

    /**
     * Method for re-initialization of parameter containers
     * in case the Run ID has changed.
     */
    InitStatus ReInit();

    /**
     * Method for event loop implementation.
     * Is called by the framework every time a new event is read.
     * @param option an execution option.
     */
    void Exec(Option_t*);

    /**
     * A method for finish of processing of an event.
     * Is called by the framework for each event after executing
     * the tasks.
     */
    void FinishEvent();

    /**
     * Method for setting the trigger value.
     * @param trigger 1 - onspill, 2 - offspill, -1 - all events.
     */
    inline void SetTrigger(Int_t trigger) { fTrigger = trigger; }

    /**
     * Method for setting the number of TofD modules.
     */
    void SetNofModules(Int_t, Int_t);

    // Method to setup online mode
    void SetOnline(Bool_t option) { fOnline = option; }

  private:
    void SetParameter();
    size_t GetCalLookupIndex(R3BTofdMappedData const&) const;

    R3BTofDMappingPar* fMapPar;

    TClonesArray* fMappedItems;        /**< Array with mapped items - input data. */
    TClonesArray* fMappedTriggerItems; /**< Array with mapped items - input data. */
    TClonesArray* fCalItems;           /**< Array with cal items - output data. */
    TClonesArray* fCalTriggerItems;    /**< Array with cal trigger items - output data. */

    R3BTCalPar* fTcalPar; /**< TCAL parameter container. */
    UInt_t fNofTcalPars;  /**< Number of modules in parameter file. */

    UInt_t fNofPlanes;
    UInt_t fPaddlesPerPlane; /**< Number of paddles per plane. */
    Double_t fClockFreq;     /**< Clock cycle in [ns]. */
    R3BEventHeader* header;  /**< Event header. */
    Int_t fTrigger;          /**< Trigger value. */
    Bool_t fOnline;          // Don't store data for online

    // Fast lookup for matching mapped data.
    std::vector<std::vector<R3BTofdCalData*>> fCalLookup;

    R3BTofdCalData* AddTCalData(UInt_t detid, UInt_t barid, UInt_t sideid, Double_t lead_time, Double_t trail_time);

    R3BTofdCalData* AddTriggerTCalData(UInt_t detid, UInt_t barid, Double_t lead_time);

  public:
    ClassDef(R3BTofDMapped2Cal, 1)
};
