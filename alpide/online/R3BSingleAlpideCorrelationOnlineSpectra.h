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
// -----      R3BSingleAlpideCorrelationOnlineSpectra     -----
// -----    Created 28/10/22 by J.L. Rodriguez-Sanchez    -----
// ----- Fill Single ALPIDE correlation online histograms -----
// ------------------------------------------------------------

#pragma once

#include <FairTask.h>
#include <TCanvas.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TMath.h>

class TClonesArray;
class R3BEventHeader;

class R3BSingleAlpideCorrelationOnlineSpectra : public FairTask
{
  public:
    /**
     * Default constructor.
     * Creates an instance of the task with default parameters.
     */
    R3BSingleAlpideCorrelationOnlineSpectra();

    /**
     * Standard constructor.
     * Creates an instance of the task.
     * @param name a name of the task.
     * @param iVerbose a verbosity level.
     * @param namedet1 a name of the detector 1.
     * @param namedet2 a name of the detector 2.
     */
    R3BSingleAlpideCorrelationOnlineSpectra(const TString& name,
                                            Int_t iVerbose = 1,
                                            const TString& namedet1 = "1",
                                            const int id1 = 1,
                                            const TString& namedet2 = "2",
                                            const int id2 = 2);

    /**
     * Destructor.
     * Frees the memory used by the object.
     */
    virtual ~R3BSingleAlpideCorrelationOnlineSpectra();

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

    /**
     * Methods to clean histograms.
     */
    virtual void Reset_Histo();

  private:
    TClonesArray* fHitItems; /**< Array with hit items. */

    R3BEventHeader* header; /**< Event header.      */
    Int_t fNEvents;         /**< Event counter.     */
    Int_t fId1, fId2;
    TString fNameDet1, fNameDet2;

    // Canvas
    TCanvas* cHit2D;

    // Histograms for hit data
    TH1F* fh1_hitxdif;
    TH1F* fh1_hitydif;
    TH2F* fh2_hitx;
    TH2F* fh2_hity;

  public:
    ClassDef(R3BSingleAlpideCorrelationOnlineSpectra, 1)
};
