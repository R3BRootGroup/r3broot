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

class R3BMwpcVsMusliCorrelationOnlineSpectra : public FairTask
{
  public:
    /**
     * Default constructor.
     * Creates an instance of the task with default parameters.
     */
    R3BMwpcVsMusliCorrelationOnlineSpectra();

    /**
     * Standard constructor.
     * Creates an instance of the task.
     * @param name a name of the task.
     * @param iVerbose a verbosity level.
     * @param namedet1 a name of the detector 1.
     * @param namedet2 a name of the detector 2.
     */
    R3BMwpcVsMusliCorrelationOnlineSpectra(const TString& name, Int_t iVerbose = 1, const TString& namedet = "Mwpc2");

    /**
     * Destructor.
     * Frees the memory used by the object.
     */
    virtual ~R3BMwpcVsMusliCorrelationOnlineSpectra() = default;

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
    TClonesArray* fCalItemsMwpc = nullptr;  /**< Array with Mwpc cal items. */
    TClonesArray* fMapItemsMusli = nullptr; /**< Array with Musli cal items. */
    TClonesArray* fCalItemsMusli = nullptr; /**< Array with Musli cal items. */

    // check for trigger should be done globablly (somewhere else)
    int fNEvents = 0;
    TString fNameDet;

    // Canvas
    TCanvas* cQvsPad;
    TCanvas* cPosX;

    // Histograms for cal and hit data
    TH2F* fh2_map_QmaxVsPadXdown;
    TH2F* fh2_map_QmaxVsPadXup;
    TH2F* fh2_map_QmaxVsPadY;
    TH2F* fh2_map_padxup;
    TH2F* fh2_map_padxdown;

  public:
    ClassDefOverride(R3BMwpcVsMusliCorrelationOnlineSpectra, 1)
};
