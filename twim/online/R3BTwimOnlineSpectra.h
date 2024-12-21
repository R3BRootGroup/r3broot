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
// -----               R3BTwimOnlineSpectra               -----
// -----    Created 06/10/19 by J.L. Rodriguez-Sanchez    -----
// -----           Fill TWIM online histograms            -----
// ------------------------------------------------------------

#pragma once

#include <FairTask.h>
#include <R3BLogger.h>
#include <TCanvas.h>
#include <TH1.h>
#include <TH2F.h>
#include <TMath.h>
#include <array>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <sstream>

class TClonesArray;
class R3BEventHeader;

/**
 * This taks reads TWIM data and plots online histograms
 */
class R3BTwimOnlineSpectra : public FairTask
{
  public:
    /**
     * Default constructor.
     * Creates an instance of the task with default parameters.
     */
    R3BTwimOnlineSpectra();

    /**
     * Standard constructor.
     * Creates an instance of the task.
     * @param name a name of the task.
     * @param iVerbose a verbosity level.
     */
    R3BTwimOnlineSpectra(const TString& name, Int_t iVerbose = 1);

    /**
     * Destructor.
     * Frees the memory used by the object.
     */
    virtual ~R3BTwimOnlineSpectra() = default;

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

    void SetExpId(Int_t exp)
    {
        R3BLOG(info, "fExpId is set locally. Original:" << fExpId << ", New value:" << exp);
        R3BLOG(info, "Using R3BEventHeader::SetExpId() is recommended instead.");
        fExpId = exp;
    }

  private:
    void s444_s467();
    void s455();

    TClonesArray* fMappedItemsTwim = nullptr; /**< Array with mapped items. */
    TClonesArray* fCalItemsTwim = nullptr;    /**< Array with cal items. */
    TClonesArray* fHitItemsTwim = nullptr;    /**< Array with hit items. */
    TClonesArray* fHitItemsMwpc3 = nullptr;   /**< Array with mw3-hit items. */
    TClonesArray* fHitItemsTofW = nullptr;    /**< Array with tofw-hit items. */

    // check for trigger should be done globablly (somewhere else)
    int fNEvents = 0;
    int fExpId = 0;
    int fNbSections = 1;
    int fNbAnodes = 16;
    int fNbTref = 1;
    int fNbTrig = 1;

    // Canvas
    TCanvas** cTwimMap_E;
    TCanvas** cTwimMap_T;
    TCanvas** cTwimMult;
    TCanvas** cTwimMap_EvsDT;
    TCanvas* cTwimMap_EsumvsDT;
    TCanvas** cTwimMap_DT;
    TCanvas** cTwimCal_Pos;
    TCanvas *cTwimMap_ESum, *cTwimMap_ESum1, *cTwimMap_ESum2, **cTwim_DTvsDT;
    TCanvas* cTwimMap_ESum_vs_diffDT;
    TCanvas** cTwimMap_DeltaTrefTrig;
    TCanvas* cTwimTheta_vs_mwpc3x;
    TCanvas* cTwimZ_vs_mwpc3x;
    TCanvas* cTwimZs[2];
    TCanvas* cTwimZsum;
    TCanvas* cTwimZsvsTofL[4];
    TCanvas* cTwimZsvsTofR[4];
    TCanvas* cTwimZsvsTof;

    // Histograms for Mapped data
    TH1F** fh1_twimmap_E;
    TH1F** fh1_twimmap_T;
    TH1F** fh1_twimmult;
    TH1F** fh1_twimmap_DT;
    TH2F** fh2_twim_DTvsDT;
    TH1F** fh1_Twimmap_mult;
    TH1F* fh1_twim_ESum[3];
    TH2F* fh2_twim_ESum;
    TH2F** fh2_twim_ESum_vs_diffDT;
    TH2F** fh2_twim_EneRawVsDriftTime;
    TH2F** fh2_twim_EneRawSumVsDriftTime;
    TH1F** fh1_twimmap_DeltaTrefTrig;

    // Histograms for Cal data
    TH1F** fh1_Twimcal_Pos;

    // Histograms for Hit data
    // s444 and s467
    TH1F* fh1_Twimhit_z;
    TH1F* fh1_Twimhit_theta;
    TH2F* fh2_Twimhit_zvstheta;
    TH2F* fh2_TwimTheta_vs_mwpc3x;
    TH2F* fh2_TwimZ_vs_mwpc3x;
    // s455
    TH2F* fh2_Twimhit_ZrZl[2];
    TH1F* fh1_Twimhit_Zr[2];
    TH1F* fh1_Twimhit_Zl[2];
    TH2F* fh2_Twimhit_ZlvsTof[2];
    TH2F* fh2_Twimhit_ZrvsTof[2];
    TH2F* fh2_Twimhit_ZlvsTofl[24];
    TH2F* fh2_Twimhit_ZrvsTofr[24];
    TH1F* fh1_twim_ZSum[2];

  public:
    ClassDefOverride(R3BTwimOnlineSpectra, 1)
};
