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
// -----                R3BMusicOnlineSpectra             -----
// -----    Created 06/10/19  by J.L. Rodriguez-Sanchez   -----
// -----           Fill MusIC online histograms           -----
// ------------------------------------------------------------

#ifndef R3BMusicOnlineSpectra_H
#define R3BMusicOnlineSpectra_H

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

#define NbAnodesMus 8

class TClonesArray;
class R3BEventHeader;

/**
 * This taks reads MusIC data and plots online histograms
 */
class R3BMusicOnlineSpectra : public FairTask
{

  public:
    /**
     * Default constructor.
     * Creates an instance of the task with default parameters.
     */
    R3BMusicOnlineSpectra();

    /**
     * Standard constructor.
     * Creates an instance of the task.
     * @param name a name of the task.
     * @param iVerbose a verbosity level.
     */
    R3BMusicOnlineSpectra(const TString& name, Int_t iVerbose = 1);

    /**
     * Destructor.
     * Frees the memory used by the object.
     */
    virtual ~R3BMusicOnlineSpectra();

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

    void SetExpId(Int_t exp) { fExpId = exp; }

  private:
    TClonesArray* fMappedItemsMus; /**< Array with mapped items. */
    TClonesArray* fCalItemsMus;    /**< Array with cal items. */
    TClonesArray* fHitItemsMus;    /**< Array with hit items. */

    // check for trigger should be done globablly (somewhere else)
    R3BEventHeader* header; /**< Event header.      */
    Int_t fNEvents;         /**< Event counter.     */
    Int_t fExpId;
    Int_t fERange;
    Int_t fZRange;
    Double_t fE[NbAnodesMus], fT[NbAnodesMus + 2];
    Int_t multhit[NbAnodesMus + 2];

    // Canvas
    TCanvas* cMusMap_E;
    TCanvas* cMusMap_EvsDT;
    TCanvas* cMusMap_DT;
    TCanvas* cMusCal_Pos;
    TCanvas *cMusMap_ESum, *cMus_EsumvsDT, *cMusMap_ESum1, *cMusMap_ESum2, *cMusMap_DTvsDT;
    TCanvas* cMusMap_ESum_vs_diffDT;

    // Histograms for Mapped data
    TH1F* fh1_Musmap_mult;
    TH2F* fh2_Musmap_multhit;
    TH1F* fh1_Mus_treftrigger;
    TH1F* fh1_Musmap_E[NbAnodesMus];
    TH1F* fh1_Musmap_DT[NbAnodesMus];
    TH1F* fh1_Mus_ESum[3];
    TH2F* fh2_Mus_ESum;
    TH2F* fh2_Mus_ESum_vs_DT;
    TH2F* fh2_Mus_EneRawVsDriftTime[NbAnodesMus];
    TH2F* fh2_Mus_DTvsDT[NbAnodesMus - 1];
    TH2F* fh2_Mus_ESum_vs_diffDT;

    // Histograms for Cal data
    TH1F* fh1_Muscal_Pos[NbAnodesMus];

    // Histograms for Hit data
    TH1F* fh1_Mushit_z;
    TH1F* fh1_Mushit_theta;
    TH2F* fh2_Mushit_zvstheta;

  public:
    ClassDef(R3BMusicOnlineSpectra, 1)
};

#endif
