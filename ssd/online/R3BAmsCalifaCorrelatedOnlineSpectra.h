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

// --------------------------------------------------------------
// -----       R3BAmsCalifaCorrelatedOnlineSpectra             --
// -----    Created 17/02/19  by J.L. Rodriguez-Sanchez        --
// ----- Fill AMS and CALIFA correlations in online histograms --
// --------------------------------------------------------------

#ifndef R3BAmsCalifaCorrelatedOnlineSpectra_H
#define R3BAmsCalifaCorrelatedOnlineSpectra_H

#include "FairTask.h"
#include "TCanvas.h"
#include "TH2F.h"
#include "TMath.h"
#include <array>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <sstream>

class TClonesArray;
class R3BEventHeader;

/**
 *  This taks reads hit data from AMS and califa detectors and plots
 *  online histograms
 */
class R3BAmsCalifaCorrelatedOnlineSpectra : public FairTask
{

  public:
    /**
     * Default constructor.
     * Creates an instance of the task with default parameters.
     */
    R3BAmsCalifaCorrelatedOnlineSpectra();

    /**
     * Standard constructor.
     * Creates an instance of the task.
     * @param name a name of the task.
     * @param iVerbose a verbosity level.
     */
    R3BAmsCalifaCorrelatedOnlineSpectra(const TString& name, Int_t iVerbose = 1);

    /**
     * Destructor.
     * Frees the memory used by the object.
     */
    virtual ~R3BAmsCalifaCorrelatedOnlineSpectra();

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
     * Method for setting the number of AMS detectors
     */
    void SetNumDet(Int_t NbDet) { fNbDet = NbDet; }

    void Reset_AMS_CALIFA_Histo();

  private:
    TClonesArray* fHitItemsAms;    /**< Array with AMS hit items. */
    TClonesArray* fHitItemsCalifa; /**< Array with CALIFA hit items. */
    TClonesArray* fCalItemsLos;

    Double_t fClockFreq; /**< Clock cycle in [ns]. */

    // check for trigger should be done globablly (somewhere else)
    R3BEventHeader* header; /**< Event header. */
    Int_t fTrigger;         /**< Trigger value. */
    Int_t fNEvents;         /**< Event counter. */

    TCanvas *cMap1, *cMap2;
    TCanvas *cAMSangles, *cHitone, *cHittwo, *cHitPhis;
    TCanvas *cMap1e, *cMap2e;
    // TCanvas* cCal;
    // TCanvas* cHit[4];

    TH2F *fh_Califa_coinc_petal1, *fh_Califa_coinc_petal2;
    TH2F *fh_Califa_coinc_petal1e, *fh_Califa_coinc_petal2e;
    TH2F* fh_los_pos_MCFD;
    TH2F* fh_los_pos_TAMEX;
    TH1F *fh_Ams_openangle1, *fh_Ams_openangle2;
    TH2F *fh_ams_theta_phi, *fh_ams_thetatop_thetabottom, *fh_ams_thetaright_thetaleft, *fh_ams_phis;
    // TH2F* fh_Ams_energy_allCalStrips[8];
    // TH2F* fh_Ams_hit_Pos[4];
    // TH2F* fh_Ams_hit_E[4];
    // TH1F* fh_Ams_hit_Mul[4];

    // TString fAmsFile;        	      /**< Config file name. */
    Int_t fNbDet; /**< Number of AMS detectors. */

  public:
    ClassDef(R3BAmsCalifaCorrelatedOnlineSpectra, 1)
};

#endif
