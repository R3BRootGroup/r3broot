/******************************************************************************
 *   Copyright (C) 2025 GSI Helmholtzzentrum für Schwerionenforschung GmbH    *
 *   Copyright (C) 2025 Members of R3B Collaboration                          *
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
// -----             R3BActafOnlineSpectra                -----
// -----    Created 27/02/25 by J.L. Rodriguez-Sanchez    -----
// -----          Fill Actaf online histograms            -----
// ------------------------------------------------------------

#pragma once

#include <FairTask.h>
#include <TCanvas.h>
#include <vector>

class TClonesArray;
class TH1F;
class TH2F;
class R3BEventHeader;
// class R3BActafMappingPar;

class R3BActafOnlineSpectra : public FairTask
{
  public:
    /**
     * Default constructor.
     * Creates an instance of the task with default parameters.
     */
    R3BActafOnlineSpectra();

    /**
     * Standard constructor.
     * Creates an instance of the task.
     * @param name a name of the task.
     * @param iVerbose a verbosity level.
     */
    R3BActafOnlineSpectra(const TString& name, Int_t iVerbose = 1);

    /**
     * Destructor.
     * Frees the memory used by the object.
     */
    virtual ~R3BActafOnlineSpectra() = default;

    /**
     * Method for task initialization.
     * This function is called by the framework before
     * the event loop.
     * @return Initialization status. kSUCCESS, kERROR or kFATAL.
     */
    InitStatus Init() override;

    /** Virtual method ReInit **/
    InitStatus ReInit() override;

    /**
     * Method for event loop implementation.
     * Is called by the framework every time a new event is read.
     * @param option an execution option.
     */
    void Exec(Option_t* /*option*/) override;

    void SetParContainers() override;

    /**
     * A method for finish of processing of an event.
     * Is called by the framework for each event after executing
     * the tasks.
     */
    void FinishEvent() override;

    void FinishTask() override;

    /**
     * Methods to clean histograms.
     */
    virtual void Reset_Histo();

    /**
     * Method for setting the trigger value.
     */
    inline void SetTrigger(int trigger) { fTrigger = trigger; }
    inline void SetTpat(int tpat1, int tpat2)
    {
        fTpat1 = tpat1;
        fTpat2 = tpat2;
    }

  private:
    void SetParameter();

    TClonesArray* fMappedItems = nullptr;
    TClonesArray* fCalItems = nullptr;
    TClonesArray* fHitItems = nullptr;

    R3BEventHeader* header = nullptr; /**< Event header. */
    // R3BActafMappingPar* fMap_Par = nullptr; /**< Parameter container. >*/
    int fTrigger = -1; /**< Trigger value. */
    int fTpat1 = 0, fTpat2 = 0;
    unsigned long fNEvents = 0;
    int fPads = 64;
    int fChn = 16;
    int fFadcs = 4;

    std::vector<TH2F*> fh2_RawTraces;

    // TH2F* fh2_xxx;

  public:
    ClassDefOverride(R3BActafOnlineSpectra, 1)
};
