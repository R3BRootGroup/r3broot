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
// -----         R3BMwpcCorrelationOnlineSpectra          -----
// -----    Created 02/12/19 by J.L. Rodriguez-Sanchez    -----
// -----     Fill MWPC correlation online histograms      -----
// ------------------------------------------------------------

#pragma once

#include "FairTask.h"
#include "TCanvas.h"
#include "TH1.h"
#include "TH2F.h"
#include "TMath.h"
#include <array>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <sstream>

class TClonesArray;
class R3BEventHeader;

class R3BMwpcCorrelationOnlineSpectra : public FairTask
{
  public:
    /**
     * Default constructor.
     * Creates an instance of the task with default parameters.
     */
    R3BMwpcCorrelationOnlineSpectra();

    /**
     * Standard constructor.
     * Creates an instance of the task.
     * @param name a name of the task.
     * @param iVerbose a verbosity level.
     * @param namedet1 a name of the detector 1.
     * @param namedet2 a name of the detector 2.
     */
    R3BMwpcCorrelationOnlineSpectra(const TString& name,
                                    Int_t iVerbose = 1,
                                    const TString& namedet1 = "Mwpc1",
                                    const TString& namedet2 = "Mwpc2");

    /**
     * Destructor.
     * Frees the memory used by the object.
     */
    virtual ~R3BMwpcCorrelationOnlineSpectra() = default;

    /**
     * Method for task initialization.
     * This function is called by the framework before
     * the event loop.
     * @return Initialization status. kSUCCESS, kERROR or kFATAL.
     */
    InitStatus Init() override;

    /**
     * Method for event loop implementation.
     * Is called by the framework every time a new event is read.
     * @param option an execution option.
     */
    void Exec(Option_t* /*option*/) override;

    /**
     * A method for finish of processing of an event.
     * Is called by the framework for each event after executing
     * the tasks.
     */
    void FinishEvent() override;

    /**
     * Method for finish of the task execution.
     * Is called by the framework after processing the event loop.
     */
    void FinishTask() override;

    /**
     * Methods to clean histograms.
     */
    virtual void Reset_Histo();

  private:
    TClonesArray* fCalItemsMwpc1 = nullptr; /**< Array with cal items. */
    TClonesArray* fHitItemsMwpc1 = nullptr; /**< Array with hit items. */
    TClonesArray* fCalItemsMwpc2 = nullptr; /**< Array with cal items. */
    TClonesArray* fHitItemsMwpc2 = nullptr; /**< Array with hit items. */

    // check for trigger should be done globablly (somewhere else)
    int fNEvents = 0;
    TString fNameDet1, fNameDet2;

    // Canvas
    TCanvas *cMWPCCal2D, *cMWPCHit2D;

    // Histograms for cal and hit data
    TH2F* fh2_mwpc_calx;
    TH2F* fh2_mwpc_caly;
    TH2F* fh2_mwpc_hitx;
    TH2F* fh2_mwpc_hity;

  public:
    ClassDefOverride(R3BMwpcCorrelationOnlineSpectra, 1)
};
