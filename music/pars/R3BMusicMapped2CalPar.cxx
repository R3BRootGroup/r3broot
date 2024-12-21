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

// ---------------------------------------------------------------------
// -----           R3BMusicMapped2CalPar source file               -----
// -----      Created 29/01/20 by J.L. Rodriguez-Sanchez           -----
// ---------------------------------------------------------------------

// ROOT headers
#include "TClonesArray.h"
#include "TF1.h"
#include "TGraph.h"
#include "TMath.h"
#include "TRandom.h"
#include "TVector3.h"

// Fair headers
#include "FairLogger.h"
#include "FairRootManager.h"
#include "FairRunAna.h"
#include "FairRuntimeDb.h"

// R3B headers
#include "R3BCoarseTimeStitch.h"
#include "R3BMusicCalPar.h"
#include "R3BMusicMapped2CalPar.h"
#include "R3BMusicMappedData.h"
#include "R3BMwpcHitData.h"

// R3BMusicMapped2CalPar: Default Constructor --------------------------
R3BMusicMapped2CalPar::R3BMusicMapped2CalPar()
    : R3BMusicMapped2CalPar("R3BMusicMapped2CalPar", 1)
{
}

// R3BMusicMapped2CalParPar: Standard Constructor --------------------------
R3BMusicMapped2CalPar::R3BMusicMapped2CalPar(const TString& name,
                                             Int_t iVerbose,
                                             const TString& namedeta,
                                             const TString& namedetb)
    : FairTask(name, iVerbose)
    , fNumAnodes(MAX_NB_MUSICANODE)   // 8 anodes
    , fNumAnodesRef(MAX_NB_MUSICTREF) // 1 anode for TREF + 1 for trigger
    , fMaxMult(MAX_MULT_MUSIC_CAL)
    , fMinStadistics(1000)
    , fLimit_left(10000)
    , fLimit_right(24000)
    , fNumParams(3)
    , fNumPosParams(2)
    , fMaxSigma(200)
    , CalParams(NULL)
    , PosParams(NULL)
    , fCal_Par(NULL)
    , fNameDetA(namedeta)
    , fPosDetA(0.)
    , fNameDetB(namedetb)
    , fPosDetB(0.)
    , fPosMusic(0.)
    , fMusicMappedDataCA(NULL)
    , fHitItemsDetA(NULL)
    , fHitItemsDetB(NULL)
{
}

// Virtual R3BMusicMapped2CalPar: Destructor
R3BMusicMapped2CalPar::~R3BMusicMapped2CalPar() { LOG(info) << "R3BMusicMapped2CalPar: Delete instance"; }

// -----   Public method Init   --------------------------------------------
InitStatus R3BMusicMapped2CalPar::Init()
{
    LOG(info) << "R3BMusicMapped2CalPar: Init";

    // INPUT DATA
    FairRootManager* rootManager = FairRootManager::Instance();
    if (!rootManager)
    {
        return kFATAL;
    }

    fMusicMappedDataCA = dynamic_cast<TClonesArray*>(rootManager->GetObject("MusicMappedData"));
    if (!fMusicMappedDataCA)
    {
        LOG(error) << "R3BMusicMapped2CalPar: MusicMappedData not found";
        return kFATAL;
    }

    // get access to hit data of detA and detB

    if (fNameDetA == "Los")
    {
        fHitItemsDetA = dynamic_cast<TClonesArray*>(rootManager->GetObject(fNameDetA + "Hit"));
        if (!fHitItemsDetA)
        {
            LOG(error) << "R3BMusicMapped2CalPar: " + fNameDetA + "Hit not found";
            return kFATAL;
        }
    }
    else
    {
        fHitItemsDetA = dynamic_cast<TClonesArray*>(rootManager->GetObject(fNameDetA + "HitData"));
        if (!fHitItemsDetA)
        {
            LOG(error) << "R3BMusicMapped2CalPar: " + fNameDetA + "HitData not found";
            return kFATAL;
        }
    }

    fHitItemsDetB = dynamic_cast<TClonesArray*>(rootManager->GetObject(fNameDetB + "HitData"));
    if (!fHitItemsDetB)
    {
        LOG(error) << "R3BMusicMapped2CalPar: " + fNameDetB + "HitData not found";
        return kFATAL;
    }

    FairRuntimeDb* rtdb = FairRuntimeDb::instance();
    if (!rtdb)
    {
        return kFATAL;
    }

    fCal_Par = dynamic_cast<R3BMusicCalPar*>(rtdb->getContainer("musicCalPar"));
    if (!fCal_Par)
    {
        LOG(error) << "R3BMusicMapped2CalPar:: Couldn't get handle on musicCalPar container";
        return kFATAL;
    }

    // Definition of a time stich object to correlate VFTX times
    fTimeStitch = new R3BCoarseTimeStitch();

    // Define TGraph for fits
    char Name1[255];
    fg_anode = new TGraph*[fNumAnodes];
    for (Int_t i = 0; i < fNumAnodes; i++)
    {
        fg_anode[i] = new TGraph();
        sprintf(Name1, "fg1_Anode_%d", i + 1);
        fg_anode[i]->SetName(Name1);
        fg_anode[i]->SetTitle(Name1);
        fg_anode[i]->SetFillColor(1);
        fg_anode[i]->SetLineColor(0);
        fg_anode[i]->SetMarkerColor(4);
        fg_anode[i]->SetMarkerStyle(20);
        fg_anode[i]->SetMarkerSize(1.2);
    }

    return kSUCCESS;
}

// -----   Public method ReInit   ----------------------------------------------
InitStatus R3BMusicMapped2CalPar::ReInit() { return kSUCCESS; }

// -----   Public method Execution   --------------------------------------------
void R3BMusicMapped2CalPar::Exec(Option_t*)
{
    // Reading the Input -- Mapped Data --
    Int_t nHits = fMusicMappedDataCA->GetEntriesFast();
    Int_t nHitsA = fHitItemsDetA->GetEntriesFast();
    Int_t nHitsB = fHitItemsDetB->GetEntriesFast();
    if (nHits < 3 || nHitsA != 1 || nHitsB != 1)
        return;

    TVector3 PosDetA(0., 0., fPosDetA);
    TVector3 PosDetB(0., 0., fPosDetB);
    R3BMwpcHitData** hitMwAData = new R3BMwpcHitData*[nHitsA];
    for (Int_t i = 0; i < nHitsA; i++)
    {
        hitMwAData[i] = dynamic_cast<R3BMwpcHitData*>(fHitItemsDetA->At(i));
        PosDetA.SetX(hitMwAData[i]->GetX());
        // LOG(info) <<hitMwAData[i]->GetX();
    }
    R3BMwpcHitData** hitMwBData = new R3BMwpcHitData*[nHitsB];
    for (Int_t i = 0; i < nHitsB; i++)
    {
        hitMwBData[i] = dynamic_cast<R3BMwpcHitData*>(fHitItemsDetB->At(i));
        PosDetB.SetX(hitMwBData[i]->GetX());
        // LOG(info) <<hitMwBData[i]->GetX();
    }

    R3BMusicMappedData** mappedData = new R3BMusicMappedData*[nHits];
    UShort_t anodeId = 0;

    for (Int_t i = 0; i < (fNumAnodes + fNumAnodesRef); i++)
    {
        mulanode[i] = 0;
        for (Int_t j = 0; j < fMaxMult; j++)
        {
            energy[j][i] = 0.;
            dtime[j][i] = 0.;
        }
    }

    for (Int_t i = 0; i < nHits; i++)
    {
        mappedData[i] = dynamic_cast<R3BMusicMappedData*>(fMusicMappedDataCA->At(i));
        anodeId = mappedData[i]->GetAnodeID();

        if (anodeId < fNumAnodes)
        {
            energy[mulanode[anodeId]][anodeId] = mappedData[i]->GetEnergy();
            dtime[mulanode[anodeId]][anodeId] = mappedData[i]->GetTime();
            mulanode[anodeId]++;
        }
        else if (anodeId >= fNumAnodes)
        { // Ref. Time and trigger
            dtime[mulanode[anodeId]][anodeId] = mappedData[i]->GetTime();
            mulanode[anodeId]++;
        }
    }

    // Fill data only if there are trigger and TREF signals
    if (mulanode[fNumAnodes] == 1 && mulanode[fNumAnodes + 1] == 1)
    {
        TF1* fa = new TF1("fa", "pol1", fPosDetA, fPosDetB);
        fa->SetParameter(0, PosDetA.X());
        fa->SetParameter(1, (PosDetB - PosDetA).X() / (fPosDetB - fPosDetA));
        for (Int_t i = 0; i < fNumAnodes; i++)
        {
            for (Int_t j = 0; j < mulanode[fNumAnodes]; j++)
                for (Int_t k = 0; k < mulanode[i]; k++)
                {
                    if (energy[k][i] > 0.)
                    { // Anode is 50mm, first anode is at 175mm with respect to the center of music detector
                        fg_anode[i]->SetPoint(fg_anode[i]->GetN() + 1,
                                              fTimeStitch->GetTime(dtime[k][i] - dtime[j][fNumAnodes], "vftx", "vftx"),
                                              fa->Eval(fPosMusic - 175.0 + i * 50.0));
                    }
                }
        }
    }
    if (mappedData)
        delete[] mappedData;
    if (hitMwAData)
        delete[] hitMwAData;
    if (hitMwBData)
        delete[] hitMwBData;
    return;
}

void R3BMusicMapped2CalPar::FinishTask()
{
    fCal_Par->SetNumAnodes(fNumAnodes);
    fCal_Par->SetNumParamsEFit(fNumParams);
    fCal_Par->SetNumParamsPosFit(fNumPosParams);
    fCal_Par->GetAnodeCalParams()->Set(fNumParams * fNumAnodes);
    fCal_Par->GetPosParams()->Set(fNumPosParams * fNumAnodes);

    TF1* fit = new TF1("fit", "pol1", fLimit_left, fLimit_right);
    fit->SetLineColor(2);
    for (Int_t i = 0; i < fNumAnodes; i++)
    {
        if (fg_anode[i]->GetN() > fMinStadistics)
        {
            fCal_Par->SetInUse(1, i + 1);
            fg_anode[i]->Fit("fit", "QR0");
            Double_t par[fNumPosParams];
            fit->GetParameters(&par[0]);
            fCal_Par->SetPosParams(par[0], i * fNumPosParams);
            fCal_Par->SetPosParams(par[1], i * fNumPosParams + 1);
        }
        else
            fCal_Par->SetAnodeCalParams(-1.0, i * fNumParams + 1);
        // fg_anode[i]->Draw("p");
        // fit->Draw("same");
        fg_anode[i]->Write();
    }
    fCal_Par->setChanged();
}

ClassImp(R3BMusicMapped2CalPar);
