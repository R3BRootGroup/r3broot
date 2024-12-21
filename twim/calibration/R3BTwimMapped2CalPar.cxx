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
// -----            R3BTwimMapped2CalPar source file               -----
// -----      Created 29/01/20 by J.L. Rodriguez-Sanchez           -----
// ---------------------------------------------------------------------

// ROOT headers
#include "TClonesArray.h"
#include "TF1.h"
#include "TGraph.h"
#include "TMath.h"
#include "TRandom.h"
#include "TVector3.h"

// FAIR headers
#include "FairLogger.h"
#include "FairRootManager.h"
#include "FairRunAna.h"
#include "FairRuntimeDb.h"

// Twim headers
#include "R3BEventHeader.h"
#include "R3BLogger.h"
#include "R3BMwpcHitData.h"
#include "R3BTwimCalPar.h"
#include "R3BTwimMapped2CalPar.h"
#include "R3BTwimMappedData.h"

// R3BTwimMapped2CalPar: Default Constructor --------------------------
R3BTwimMapped2CalPar::R3BTwimMapped2CalPar()
    : R3BTwimMapped2CalPar("R3BTwimMapped2CalPar", 1, "Mwpc1", "Mwpc2")
{
}

// R3BTwimMapped2CalParPar: Standard Constructor --------------------------
R3BTwimMapped2CalPar::R3BTwimMapped2CalPar(const TString& name,
                                           Int_t iVerbose,
                                           const TString& namedeta,
                                           const TString& namedetb)
    : FairTask(name, iVerbose)
    , fNumSec(4)
    , fNumAnodes(16)   // 16 anodes
    , fNumAnodesRef(2) // 2 anode for TREF
    , fMaxMult(20)
    , fMinStadistics(1000)
    , fLimit_left(0)
    , fLimit_right(24000)
    , fNumParams(3)
    , fNumPosParams(2)
    , fMaxSigma(200)
    , CalParams(NULL)
    , PosParams(NULL)
    , fCal_Par(NULL)
    , fNameDetA(namedeta)
    , fPosMwpcA(0.)
    , fNameDetB(namedetb)
    , fPosMwpcB(0.)
    , fPosTwim(0.)
    , fTwimMappedDataCA(NULL)
    , fHitItemsMwpcA(NULL)
    , fHitItemsMwpcB(NULL)
    , fExpId(0)
{
}

// Virtual R3BTwimMapped2CalPar: Destructor
R3BTwimMapped2CalPar::~R3BTwimMapped2CalPar() { R3BLOG(info, "Delete instance"); }

// -----   Public method Init   --------------------------------------------
InitStatus R3BTwimMapped2CalPar::Init()
{
    R3BLOG(info, "Init");

    // INPUT DATA
    FairRootManager* rootManager = FairRootManager::Instance();
    if (!rootManager)
    {
        return kFATAL;
    }

    header = dynamic_cast<R3BEventHeader*>(rootManager->GetObject("EventHeader."));
    if (!header)
        header = dynamic_cast<R3BEventHeader*>(rootManager->GetObject("R3BEventHeader"));
    if (fExpId == 0) // Obtain global ExpId if it's not set locally.
    {
        fExpId = header->GetExpId();
        R3BLOG(info, "fExpId :" << fExpId);
    }

    fTwimMappedDataCA = dynamic_cast<TClonesArray*>(rootManager->GetObject("TwimMappedData"));
    if (!fTwimMappedDataCA)
    {
        R3BLOG(error, "TwimMappedData not found");
        return kFATAL;
    }

    // get access to hit data of mwpcs
    fHitItemsMwpcA = dynamic_cast<TClonesArray*>(rootManager->GetObject(fNameDetA + "HitData"));
    if (!fHitItemsMwpcA)
    {
        R3BLOG(error, fNameDetA << "HitData not found");
        return kFATAL;
    }

    fHitItemsMwpcB = dynamic_cast<TClonesArray*>(rootManager->GetObject(fNameDetB + "HitData"));
    if (!fHitItemsMwpcB)
    {
        R3BLOG(error, fNameDetB << "HitData not found");
        return kFATAL;
    }

    FairRuntimeDb* rtdb = FairRuntimeDb::instance();
    if (!rtdb)
    {
        return kFATAL;
    }

    fCal_Par = dynamic_cast<R3BTwimCalPar*>(rtdb->getContainer("twimCalPar"));
    if (!fCal_Par)
    {
        R3BLOG(error, "Couldn't get handle on twimCalPar container");
        return kFATAL;
    }

    if (fExpId == 444 || fExpId == 467)
    {
        fNumSec = 1;
        fNumAnodes = 16;   // 16 anodes
        fNumAnodesRef = 2; // 2 anode for TREF
        fMaxMult = 10;
    }
    else
    {
        fNumSec = 4;
        fNumAnodes = 16;   // 16 anodes
        fNumAnodesRef = 1; // 1 anode for TREF
        fMaxMult = 20;
    }

    // Define TGraph for fits
    char Name1[255];
    fg_anode = new TGraph*[fNumSec * fNumAnodes];
    for (Int_t s = 0; s < fNumSec; s++)
        for (Int_t i = 0; i < fNumAnodes; i++)
        {
            fg_anode[s * fNumAnodes + i] = new TGraph();
            sprintf(Name1, "fg_sec%d_Anode_%d", s + 1, i + 1);
            fg_anode[s * fNumAnodes + i]->SetName(Name1);
            fg_anode[s * fNumAnodes + i]->SetTitle(Name1);
            fg_anode[s * fNumAnodes + i]->SetFillColor(1);
            fg_anode[s * fNumAnodes + i]->SetLineColor(0);
            fg_anode[s * fNumAnodes + i]->SetMarkerColor(4);
            fg_anode[s * fNumAnodes + i]->SetMarkerStyle(20);
            fg_anode[s * fNumAnodes + i]->SetMarkerSize(1.2);
        }

    return kSUCCESS;
}

// -----   Public method ReInit   ----------------------------------------------
InitStatus R3BTwimMapped2CalPar::ReInit() { return kSUCCESS; }

// -----   Public method Execution   --------------------------------------------
void R3BTwimMapped2CalPar::Exec(Option_t* option)
{
    // Reading the Input -- Mapped Data --
    Int_t nHits = fTwimMappedDataCA->GetEntriesFast();
    Int_t nHitsA = fHitItemsMwpcA->GetEntriesFast();
    Int_t nHitsB = fHitItemsMwpcB->GetEntriesFast();
    if (nHits < 3 || nHitsA != 1 || nHitsB != 1)
        return;

    TVector3 PosMwpcA(0., 0., fPosMwpcA);
    TVector3 PosMwpcB(0., 0., fPosMwpcB);
    R3BMwpcHitData** hitMwAData = new R3BMwpcHitData*[nHitsA];
    for (Int_t i = 0; i < nHitsA; i++)
    {
        hitMwAData[i] = dynamic_cast<R3BMwpcHitData*>(fHitItemsMwpcA->At(i));
        PosMwpcA.SetX(hitMwAData[i]->GetX());
        // R3BLOG(info,hitMwAData[i]->GetX());
    }
    R3BMwpcHitData** hitMwBData = new R3BMwpcHitData*[nHitsB];
    for (Int_t i = 0; i < nHitsB; i++)
    {
        hitMwBData[i] = dynamic_cast<R3BMwpcHitData*>(fHitItemsMwpcB->At(i));
        PosMwpcB.SetX(hitMwBData[i]->GetX());
        // R3BLOG(info,hitMwBData[i]->GetX());
    }

    R3BTwimMappedData** mappedData = new R3BTwimMappedData*[nHits];
    Int_t secId = 0;
    Int_t anodeId = 0;

    for (Int_t s = 0; s < fNumSec; s++)
        for (Int_t i = 0; i < (fNumAnodes + fNumAnodesRef); i++)
        {
            mulanode[s][i] = 0;
            for (Int_t j = 0; j < fMaxMult; j++)
            {
                fE[s][j][i] = 0.;
                fDT[s][j][i] = 0.;
            }
        }

    for (Int_t i = 0; i < nHits; i++)
    {
        mappedData[i] = dynamic_cast<R3BTwimMappedData*>(fTwimMappedDataCA->At(i));
        secId = mappedData[i]->GetSecID() - 1;
        anodeId = mappedData[i]->GetAnodeID() - 1;

        if (anodeId < fNumAnodes)
        {

            fE[secId][mulanode[secId][anodeId]][anodeId] = mappedData[i]->GetEnergy();
            fDT[secId][mulanode[secId][anodeId]][anodeId] = mappedData[i]->GetTime();
            mulanode[secId][anodeId]++;
        }
        else if (anodeId >= fNumAnodes)
        {                                                                             // Ref. Time
            fDT[secId][mulanode[secId][anodeId]][anodeId] = mappedData[i]->GetTime(); // Ref. Time
            mulanode[secId][anodeId]++;
        }
    }

    // Fill data only if there is TREF signal
    for (Int_t s = 0; s < fNumSec; s++)
        if (TMath::Abs(PosMwpcA.X() + (PosMwpcB - PosMwpcA).X() / (fPosMwpcB - fPosMwpcA) * fPosMwpcA) < 100.)
            if (mulanode[s][fNumAnodes] == 1 && mulanode[s][fNumAnodes + 1] == 1)
            {
                TF1* fa = new TF1("fa", "pol1", fPosMwpcA, fPosMwpcB);
                fa->SetParameter(0, PosMwpcA.X() - (PosMwpcB - PosMwpcA).X() / (fPosMwpcB - fPosMwpcA) * fPosMwpcA);
                fa->SetParameter(1, (PosMwpcB - PosMwpcA).X() / (fPosMwpcB - fPosMwpcA));
                for (Int_t i = 0; i < fNumAnodes; i++)
                {
                    // for (Int_t j = 0; j < mulanode[s][fNumAnodes]; j++)
                    for (Int_t j = 0; j < 1; j++) // Mult 1
                        // for (Int_t k = 0; k < mulanode[s][i]; k++)
                        for (Int_t k = 0; k < 1; k++) // Mult 1
                        {
                            if (fE[s][k][i] > 0. && (fExpId == 444 || fExpId == 467))
                            { // Anode is 25mm, first anode is at -187.5mm with respect to the center of twim detector
                                if (i < fNumAnodes / 2)
                                    fg_anode[s * fNumAnodes + i]->SetPoint(fg_anode[s * fNumAnodes + i]->GetN() + 1,
                                                                           fDT[s][k][i] - fDT[s][j][fNumAnodes],
                                                                           fa->Eval(fPosTwim - 187.5 + i * 25.0));
                                else
                                    fg_anode[s * fNumAnodes + i]->SetPoint(fg_anode[s * fNumAnodes + i]->GetN() + 1,
                                                                           fDT[s][k][i] - fDT[s][j][fNumAnodes + 1],
                                                                           fa->Eval(fPosTwim - 187.5 + i * 25.0));
                            }

                            if (fE[s][k][i] > 0. && fExpId == 455)
                            { // Anode is 25mm, first anode is at -187.5mm with respect to the center of twim detector

                                fg_anode[s * fNumAnodes + i]->SetPoint(fg_anode[s * fNumAnodes + i]->GetN() + 1,
                                                                       fDT[s][k][i] - fDT[s][j][fNumAnodes],
                                                                       fa->Eval(fPosTwim - 187.5 + i * 25.0));
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

void R3BTwimMapped2CalPar::FinishTask()
{
    fCal_Par->SetNumSec(fNumSec);
    fCal_Par->SetNumAnodes(fNumAnodes);
    fCal_Par->SetNumParamsEFit(fNumParams);
    fCal_Par->SetNumParamsPosFit(fNumPosParams);

    TF1* fit = new TF1("fit", "pol1", fLimit_left, fLimit_right);
    fit->SetLineColor(2);

    for (Int_t s = 0; s < fNumSec; s++)
    {
        fCal_Par->GetAnodeCalEParams(s)->Set(fNumParams * fNumAnodes);
        fCal_Par->GetAnodeCalPosParams(s)->Set(fNumPosParams * fNumAnodes);

        for (Int_t i = 0; i < fNumAnodes; i++)
        {
            if (fg_anode[s * fNumAnodes + i]->GetN() > fMinStadistics)
            {
                fCal_Par->SetInUse(1, s + 1, i + 1);
                fg_anode[s * fNumAnodes + i]->Fit("fit", "QR0");
                Double_t par[fNumPosParams];
                fit->GetParameters(&par[0]);
                fCal_Par->SetPosParams(par[0], s, i, 1);
                fCal_Par->SetPosParams(par[1], s, i, 2);
            }
            else
            {
                fCal_Par->SetAnodeCalParams(0.0, s, i, 1);
                fCal_Par->SetAnodeCalParams(-1.0, s, i, 2);
            }

            fg_anode[s * fNumAnodes + i]->Write();
        }
    }

    fCal_Par->setChanged();
}

ClassImp(R3BTwimMapped2CalPar);
