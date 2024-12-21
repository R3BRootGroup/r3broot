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
// -----             R3BCalifaJulichOnlineSpectra                 -----
// -----    Created 16/07/21  by J.L. Rodriguez-Sanchez   -----
// -----          Fill CalifaJulich online histograms             -----
// ------------------------------------------------------------

#ifndef R3BCalifaJulichOnlineSpectra_H
#define R3BCalifaJulichOnlineSpectra_H

#include "FairTask.h"

#include "TCanvas.h"
#include "TMath.h"
#include <Rtypes.h>
#include <array>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <sstream>
#include <vector>

class TClonesArray;
class R3BEventHeader;
class TH1F;
class TH2F;

/**
 * This taks reads CalifaJulich data and plots online histograms
 */
class R3BCalifaJulichOnlineSpectra : public FairTask
{

  public:
    /**
     * Default constructor.
     * Creates an instance of the task with default parameters.
     */
    R3BCalifaJulichOnlineSpectra();

    /**
     * Standard constructor.
     * Creates an instance of the task.
     * @param name a name of the task.
     * @param iVerbose a verbosity level.
     */
    R3BCalifaJulichOnlineSpectra(const TString& name, Int_t iVerbose = 1);

    /**
     * Destructor.
     * Frees the memory used by the object.
     */
    virtual ~R3BCalifaJulichOnlineSpectra();

    /**
     * Method for task initialization.
     * This function is called by the framework before
     * the event loop.
     * @return Initialization status. kSUCCESS, kERROR or kFATAL.
     */
    virtual InitStatus Init() override;

    /**
     * Method for event loop implementation.
     * Is called by the framework every time a new event is read.
     * @param option an execution option.
     */
    virtual void Exec(Option_t* option) override;

    /**
     * A method for finish of processing of an event.
     * Is called by the framework for each event after executing
     * the tasks.
     */
    virtual void FinishEvent() override;

    /**
     * Method for finish of the task execution.
     * Is called by the framework after processing the event loop.
     */
    virtual void FinishTask() override;

    /** Virtual method SetParContainers **/
    // virtual void SetParContainers();

    /**
     * Method to reset histograms
     */
    void Reset_CalifaJulich_Histo();

  private:
    TClonesArray* fMappedItemsCalifa; // Array with mapped items.
    // TClonesArray* fCalItems;    // Array with cal items.
    // TClonesArray* fHitItems;    // Array with hit items.
    TClonesArray* fMappedItemsSi; // Array with mapped items.

    R3BEventHeader* header; // Event header.
    Int_t fTrigger;         // Trigger value.
    Int_t fNEvents;         // Event counter.
    Int_t fNbDet;           // Number of AMS detectors.

    // Histograms for map data
    std::vector<TH2F*> fh2_EnergyVsStrip;

  public:
    ClassDefOverride(R3BCalifaJulichOnlineSpectra, 0)
};

#endif /* R3BCalifaJulichOnlineSpectra_H */
