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
// -----               R3BTttxOnlineSpectra               -----
// -----    Created 23/01/24 by M.Whitehead               -----
// -----           Fill TTT10 online histograms           -----
// ------------------------------------------------------------

#pragma once

#include "FairTask.h"
#include "R3BLogger.h"
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

/**
 * This taks reads TTT10 data and plots online histograms
 */
class R3BTttxOnlineSpectra : public FairTask
{

  public:
    /**
     * Default constructor.
     * Creates an instance of the task with default parameters.
     */
    R3BTttxOnlineSpectra();

    /**
     * Standard constructor.
     * Creates an instance of the task.
     * @param name a name of the task.
     * @param iVerbose a verbosity level.
     */
    R3BTttxOnlineSpectra(const TString& name, Int_t iVerbose = 1);

    /**
     * Destructor.
     * Frees the memory used by the object.
     */
    virtual ~R3BTttxOnlineSpectra();

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

    void SetExpId(Int_t exp)
    {
        R3BLOG(info, "fExpId is set locally. Original:" << fExpId << ", New value:" << exp);
        R3BLOG(info, "Using R3BEventHeader::SetExpId() is recommended instead.");
        fExpId = exp;
    }

    // Number of entries before some histograms reset
    void SetNumEntries(Int_t num) { fNbEntries = num; }

  private:
    TClonesArray* fMappedItemsTttx; /**< Array with mapped items. */
    TClonesArray* fCalItemsTttx;    /**< Array with cal items. */
    TClonesArray* fHitItemsTttx;    /**< Array with hit items. */

    // check for trigger should be done globablly (somewhere else)
    R3BEventHeader* header; /**< Event header.      */
    Int_t fNEvents;         /**< Event counter.     */
    Int_t fExpId;
    Int_t fNbDets;
    Int_t fNbStrips;
    Int_t fNbTref;
    Int_t fNbTrig;
    Int_t fNbEntries;

    // Canvas
    TCanvas* cMap;
    TCanvas* cCal;
    TCanvas* cHit;

    // Histograms for Mapped data
    std::vector<TH2F*> fh2_EnergyVsStrip_all;
    std::vector<TH2F*> fh2_EnergyVsStrip_recent;
    std::vector<TH2F*> fh2_TimeVsStrip;
    TH2F* fh2_E2VsE1;
    TH2F* fh2_E2StripVsE1Strip;
    std::vector<TH1F*> fh1_Strip;
    std::vector<TH2F*> fh2_MultVsStrip;

    // Histograms for Cal data
    std::vector<TH2F*> fh2_CalEnergyVsStrip;
    std::vector<TH2F*> fh2_CalTimeVsStrip;
    TH2F* fh2_CalMultCorr;
    TH2F* fh2_CalECorr;
    TH2F* fh2_CalTCorr;

    // Histograms for Hit data
    std::vector<TH2F*> fh2_HitEnergyVsPos;
    std::vector<TH2F*> fh2_HitTimeVsPos;
    TH2F* fh2_HitMultCorr;
    TH2F* fh2_HitECorr;
    TH2F* fh2_HitZCorr;
    TH2F* fh2_HitPosCorr;
    TH2F* fh2_Zcharge_tpat;

    // Temporal value containers
    std::vector<std::vector<Double_t>> energy;
    std::vector<std::vector<Double_t>> time;
    std::vector<std::vector<Int_t>> mult;
    std::vector<Double_t> highest_e;
    std::vector<UInt_t> highest_e_strip;

    std::vector<Double_t> cal_mult;

    std::vector<Double_t> hit_mult;

  public:
    ClassDef(R3BTttxOnlineSpectra, 1)
};
