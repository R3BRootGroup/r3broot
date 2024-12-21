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

// ----------------------------------------------------------------------
// -----                   R3BMwpc0Mapped2CalPar                    -----
// -----             Created 14/10/19 by G. García Jiménez          -----
// ----------------------------------------------------------------------

#include "R3BMwpc3Mapped2CalPar.h"
#include "R3BEventHeader.h"
#include "R3BMwpc3CalPar.h"
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

/* ---- R3BMwpc3Mapped2CalPar: Default Constructor ---- */
R3BMwpc3Mapped2CalPar::R3BMwpc3Mapped2CalPar()
    : R3BMwpc3Mapped2CalPar("R3BMwpc3Mapped2CalPar", 1)
{
}

/* ---- R3BMwpc3Mapped2CalPar: Standard Constructor ---- */
R3BMwpc3Mapped2CalPar::R3BMwpc3Mapped2CalPar(const char* name, Int_t iVerbose)
    : FairTask(name, iVerbose)
    , fNumPadX(288)
    , fNumPadY(120)
    , fNumParams(2)
    , fMinStadistics(100)
    , fMapHistos_left(0)
    , fMapHistos_right(270000)
    , fMapHistos_bins(27000)
    , fPad_Par(NULL)
    , fMwpcMappedDataCA(NULL)
{
}

/* ---- R3BMwpc3Mapped2CalPar : Destructor ---- */
R3BMwpc3Mapped2CalPar::~R3BMwpc3Mapped2CalPar() { LOG(info) << "R3BMwpc3Mapped2CalPar: Delete Instance"; }

/* ---- Public Method Init ---- */
InitStatus R3BMwpc3Mapped2CalPar::Init()
{

    LOG(info) << "R3BMwpc3Mapped2CalPar: Init";

    char name[100];

    fh_Map_q_pad = new TH1F*[fNumPadX + fNumPadY]; // Array of Histograms

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

    fMwpcMappedDataCA = dynamic_cast<TClonesArray*>(rootManager->GetObject("Mwpc3MappedData"));
    if (!fMwpcMappedDataCA)
    {
        return kFATAL;
    }

    FairRuntimeDb* rtdb = FairRuntimeDb::instance();
    if (!rtdb)
    {
        return kFATAL;
    }

    fPad_Par = dynamic_cast<R3BMwpc3CalPar*>(rtdb->getContainer("mwpc3CalPar"));
    if (!fPad_Par)
    {
        LOG(error) << "R3BMwpc3Mapped2CalPar::Init() Couldn't get handle on mwpc3CalPar container";
        return kFATAL;
    }
    return kSUCCESS;
}

/* -----   Public method ReInit   ----*/
InitStatus R3BMwpc3Mapped2CalPar::ReInit() { return kSUCCESS; }

/* ---- Public Method Exec ---- */
void R3BMwpc3Mapped2CalPar::Exec(Option_t* opt)
{

    Int_t nHits = fMwpcMappedDataCA->GetEntries();
    if (!nHits)
        return;

    R3BMwpcMappedData* MapHit;
    UChar_t planeid, padid;

    for (Int_t i = 0; i < nHits; i++)
    {
        MapHit = dynamic_cast<R3BMwpcMappedData*>(fMwpcMappedDataCA->At(i));
        planeid = MapHit->GetPlane();
        padid = MapHit->GetPad() - 1;

        // Fill Histos
        if (planeid == 1) // plane X
            fh_Map_q_pad[padid]->Fill(MapHit->GetQ());
        else if (planeid == 3) // plane Y
            fh_Map_q_pad[fNumPadX + padid]->Fill(MapHit->GetQ());
        else
            LOG(error) << "Plane " << planeid << " does not exist in MWPC3";
    }
}

/* ---- Public method Finish   ---- */
void R3BMwpc3Mapped2CalPar::FinishTask()
{

    SearchPedestals();
    // fPad_Par->printParams();
}

/* ---- Searching Pedestals ----*/

void R3BMwpc3Mapped2CalPar::SearchPedestals()
{

    LOG(info) << "R3BMwpc3Mapped2CalPar: Search pedestals";
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
            if (i < fNumPadX) // plane X
                LOG(warn) << "Histogram NO Fitted in mwpc3, plane 1 and pad " << i + 1;
            else // plane y
                LOG(warn) << "Histogram NO Fitted in mwpc3, plane 3 and pad " << i + 1 - fNumPadX;
        }
    }

    fPad_Par->setChanged();
    return;
}
ClassImp(R3BMwpc3Mapped2CalPar);
