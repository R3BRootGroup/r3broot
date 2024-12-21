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

// ---------------------------------------------------------------
// -----         R3BMwpc2Mapped2CalPar source file           -----
// -----    Created 10/10/19 by J.L. Rodriguez-Sanchez       -----
// ---------------------------------------------------------------

#include "R3BMwpc2Mapped2CalPar.h"
#include "R3BEventHeader.h"
#include "R3BMwpc2CalPar.h"
#include "R3BMwpcMappedData.h"

#include "FairLogger.h"
#include "FairRootManager.h"
#include "FairRunAna.h"
#include "FairRuntimeDb.h"
#include "TGeoManager.h"

#include "TCanvas.h"
#include "TClonesArray.h"
#include "TF1.h"
#include "TGeoMatrix.h"
#include "TGraph.h"
#include "TH1F.h"
#include "TMath.h"
#include "TObjArray.h"
#include "TRandom.h"
#include "TVector3.h"

#include <iostream>
#include <stdlib.h>

using namespace std;

// R3BMwpc2Mapped2CalPar: Default Constructor --------------------------
R3BMwpc2Mapped2CalPar::R3BMwpc2Mapped2CalPar()
    : R3BMwpc2Mapped2CalPar("R3BMwpc2Mapped2CalPar", 1)
{
}

// R3BMwpc2Mapped2CalPar: Standard Constructor --------------------------
R3BMwpc2Mapped2CalPar::R3BMwpc2Mapped2CalPar(const char* name, Int_t iVerbose)
    : FairTask(name, iVerbose)
    , fNumPadX(128)
    , fNumPadY(40)
    , fNumParams(2)
    , fMinStadistics(100)
    , fMapHistos_left(0)
    , fMapHistos_right(270000)
    , fMapHistos_bins(27000)
    , fPad_Par(NULL)
    , fMwpcMappedDataCA(NULL)
{
}

// R3BMwpc2Mapped2CalPar: Destructor ----------------------------------------
R3BMwpc2Mapped2CalPar::~R3BMwpc2Mapped2CalPar() { LOG(info) << "R3BMwpc2Mapped2CalPar: Delete instance"; }

// -----   Public method Init   --------------------------------------------
InitStatus R3BMwpc2Mapped2CalPar::Init()
{
    LOG(info) << "R3BMwpc2Mapped2CalPar: Init";

    char name[100];

    fh_Map_q_pad = new TH1F*[fNumPadX + fNumPadY];
    for (Int_t i = 0; i < fNumPadX + fNumPadY; i++)
    {
        if (i < fNumPadX)
            sprintf(name, "fh_Map_q_padx_%i", i + 1);
        else
            sprintf(name, "fh_Map_q_pady_%i", i + 1 - fNumPadX);
        fh_Map_q_pad[i] = new TH1F(name, name, fMapHistos_bins, fMapHistos_left, fMapHistos_right);
    }

    FairRootManager* rootManager = FairRootManager::Instance();
    if (!rootManager)
    {
        return kFATAL;
    }

    fMwpcMappedDataCA = dynamic_cast<TClonesArray*>(rootManager->GetObject("Mwpc2MappedData"));
    if (!fMwpcMappedDataCA)
    {
        return kFATAL;
    }

    FairRuntimeDb* rtdb = FairRuntimeDb::instance();
    if (!rtdb)
    {
        return kFATAL;
    }

    fPad_Par = dynamic_cast<R3BMwpc2CalPar*>(rtdb->getContainer("mwpc2CalPar"));
    if (!fPad_Par)
    {
        LOG(error) << "R3BMwpc2Mapped2CalPar::Init() Couldn't get handle on mwpc2CalPar container";
        return kFATAL;
    }
    return kSUCCESS;
}

// -----   Public method ReInit   --------------------------------------------
InitStatus R3BMwpc2Mapped2CalPar::ReInit() { return kSUCCESS; }

// -----   Public method Exec   --------------------------------------------
void R3BMwpc2Mapped2CalPar::Exec(Option_t* opt)
{

    Int_t nHits = fMwpcMappedDataCA->GetEntriesFast();
    if (nHits == 0)
        return;

    R3BMwpcMappedData* MapHit;
    UChar_t planeid, padid;

    for (Int_t i = 0; i < nHits; i++)
    {
        MapHit = dynamic_cast<R3BMwpcMappedData*>(fMwpcMappedDataCA->At(i));
        planeid = MapHit->GetPlane();
        padid = MapHit->GetPad() - 1;

        // Fill Histos
        if (planeid == 1) // plane X down
            fh_Map_q_pad[padid]->Fill(MapHit->GetQ());
        else if (planeid == 2) // plane X up
            fh_Map_q_pad[padid + fNumPadX / 2]->Fill(MapHit->GetQ());
        else if (planeid == 3) // plane Y
            fh_Map_q_pad[fNumPadX + padid]->Fill(MapHit->GetQ());
        else
            LOG(error) << "Plane " << planeid << " does not exist in MWPC2";
    }
}

// ---- Public method Finish   --------------------------------------------------
void R3BMwpc2Mapped2CalPar::FinishTask()
{

    SearchPedestals();
    // fPad_Par->printParams();
}

//---- Search Pedestals   -------------------------------------------------------
void R3BMwpc2Mapped2CalPar::SearchPedestals()
{

    LOG(info) << "R3BMwpc2Mapped2CalPar: Search pedestals";

    Int_t numPars = fNumParams;

    fPad_Par->SetNumPadsX(fNumPadX);
    fPad_Par->SetNumPadsY(fNumPadY);
    fPad_Par->SetNumParametersFit(fNumParams);
    fPad_Par->GetPadCalParams()->Set(fNumParams * (fNumPadX + fNumPadY));

    Int_t nbpad = 0;
    for (Int_t i = 0; i < fNumPadX + fNumPadY; i++)
    {
        nbpad = i * fNumParams;
        if (fh_Map_q_pad[i]->GetEntries() > fMinStadistics)
        {
            Int_t tmp = fh_Map_q_pad[i]->GetMaximumBin();
            TF1* f1 = new TF1("f1", "gaus", fMapHistos_left, fMapHistos_right);
            fh_Map_q_pad[i]->Fit("f1", "QON", "", tmp - 80, tmp + 80);
            fPad_Par->SetPadCalParams(f1->GetParameter(1), nbpad);
            fPad_Par->SetPadCalParams(f1->GetParameter(2), nbpad + 1);
        }
        else
        {
            fPad_Par->SetPadCalParams(-1, nbpad); // dead pad
            fPad_Par->SetPadCalParams(0, nbpad + 1);
            if (i < fNumPadX / 2) // plane X down
                LOG(warn) << "Histogram NO Fitted in mwpc2, plane 1 and pad (1-64): " << i + 1;
            else if (i < fNumPadX) // plane X up
                LOG(warn) << "Histogram NO Fitted in mwpc2, plane 2 and pad (1-64): " << i + 1 - fNumPadX / 2;
            else // plane y
                LOG(warn) << "Histogram NO Fitted in mwpc2, plane 3 and pad (1-40): " << i + 1 - fNumPadX;
        }
    }
    fPad_Par->setChanged();
    return;
}

ClassImp(R3BMwpc2Mapped2CalPar)
