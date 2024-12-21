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

#ifndef R3BMusliVsMwpcOnlineSpectra_H
#define R3BMusliVsMwpcOnlineSpectra_H

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
class R3BTGeoPar;

/**
 * This taks reads MusIC data and plots online histograms
 */
class R3BMusliVsMwpcOnlineSpectra : public FairTask
{

  public:
    /**
     * Default constructor.
     * Creates an instance of the task with default parameters.
     */
    R3BMusliVsMwpcOnlineSpectra();

    /**
     * Standard constructor.
     * Creates an instance of the task.
     * @param name a name of the task.
     * @param iVerbose a verbosity level.
     */
    R3BMusliVsMwpcOnlineSpectra(const TString& name,
                                Int_t iVerbose = 1,
                                const TString& namedet1 = "Mwpc1",
                                const TString& namedet2 = "Mwpc0");

    /**
     * Destructor.
     * Frees the memory used by the object.
     */
    virtual ~R3BMusliVsMwpcOnlineSpectra();

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

    virtual void SetParContainers();

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
    TClonesArray* fMapMusli;    /**< Array with mapped level items. */
    TClonesArray* fCalMusli;    /**< Array with cal level items. */
    TClonesArray* fHitMusli;    /**< Array with hit level items. */
    TClonesArray* fHitMwpcDet1; /**< Array with hit items of mwpc */
    TClonesArray* fHitMwpcDet2; /**< Array with hit items of mwpc */
    R3BTGeoPar* fMw1Geo_Par;
    R3BTGeoPar* fMw2Geo_Par;

    // check for trigger should be done globablly (somewhere else)
    R3BEventHeader* header; /**< Event header.      */
    Int_t fNEvents;         /**< Event counter.     */
    TString fNameDet1;
    TString fNameDet2;

    // Canvas - MapMusli and HitMwpc
    TCanvas* cMusliMap_DTvsX1;
    TCanvas* cMusliMap_DTvsX2;
    TH2D* fh2_mapDTvsX1[15];
    TH2D* fh2_mapDTvsX2[15];

    // Canvas - CalMusli and HitMwpc
    TCanvas* cMusliCal_DTvsX1;
    TCanvas* cMusliCal_DTvsX2;

    // TH2 - CalMusli and HitMwpc
    TH2D* fh2_calDTvsX1[15];
    TH2D* fh2_calDTvsX2[15];

    // Canvas - HitMusli and HitMwpc
    TCanvas* cMusliHit_EvsX1;
    TCanvas* cMusliHit_ZvsX1;
    TCanvas* cMusliHit_XvsX1;
    TCanvas* cMusliHit_TvsX1;

    TCanvas* cMusliHit_EvsX2;
    TCanvas* cMusliHit_ZvsX2;
    TCanvas* cMusliHit_XvsX2;
    TCanvas* cMusliHit_TvsX2;

    // TH2 - HitMusli and HitMwpc
    TH2D* fh2_hitEvsX1[4];
    TH2D* fh2_hitZvsX1[4];
    TH2D* fh2_hitXvsX1[4];
    TH2D* fh2_hitTvsX1;

    TH2D* fh2_hitEvsX2[4];
    TH2D* fh2_hitZvsX2[4];
    TH2D* fh2_hitXvsX2[4];
    TH2D* fh2_hitTvsX2;

  public:
    ClassDef(R3BMusliVsMwpcOnlineSpectra, 1)
};

#endif
