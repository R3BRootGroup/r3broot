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

#include "FairLogger.h"
#include "FairRootManager.h"
#include "FairRunOnline.h"
#include "FairRuntimeDb.h"
#include "TCanvas.h"
#include "TClonesArray.h"
#include "TFolder.h"
#include "TH1.h"
#include "TH2.h"
#include "THttpServer.h"
#include "TMath.h"

#include "R3BEventHeader.h"
#include "R3BFrsSciPosCalData.h"
#include "R3BFrsSciTofCalData.h"
#include "R3BMusliCalData.h"
#include "R3BMusliHitData.h"
#include "R3BMusliMappedData.h"
#include "R3BOnlineSpectraFrsSciVsMusli.h"

R3BOnlineSpectraFrsSciVsMusli::R3BOnlineSpectraFrsSciVsMusli()
    : R3BOnlineSpectraFrsSciVsMusli("FrsSciOnlineSpectra", 1)
{
}

R3BOnlineSpectraFrsSciVsMusli::R3BOnlineSpectraFrsSciVsMusli(const char* name, Int_t iVerbose)
    : FairTask(name, iVerbose)
    , fFrsSci_PosCal(NULL)
    , fFrsSci_TofCal(NULL)
    , fMusli_Map(NULL)
    , fMusli_Cal(NULL)
    , fMusli_Hit(NULL)
    , fNEvents(0)
    , fIdS2(1)
    , fIdCaveC(2)
{
}

R3BOnlineSpectraFrsSciVsMusli::~R3BOnlineSpectraFrsSciVsMusli()
{
    LOG(debug) << "R3BOnlineSpectraFrsSciVsMusli::Destructor";
    if (fFrsSci_PosCal)
        delete fFrsSci_PosCal;
    if (fFrsSci_TofCal)
        delete fFrsSci_TofCal;
    if (fMusli_Map)
        delete fMusli_Map;
    if (fMusli_Cal)
        delete fMusli_Cal;
    if (fMusli_Hit)
        delete fMusli_Hit;
}

InitStatus R3BOnlineSpectraFrsSciVsMusli::Init()
{
    LOG(info) << "R3BOnlineSpectraFrsSciVsMusli::Init()";

    char Name1[255];

    // try to get a handle on the EventHeader. EventHeader may not be
    // present though and hence may be null. Take care when using.

    FairRootManager* mgr = FairRootManager::Instance();
    if (NULL == mgr)
        LOG(fatal) << "R3BOnlineSpectraFrsSciVsMusli::Init FairRootManager not found";

    FairRunOnline* run = FairRunOnline::Instance();
    run->GetHttpServer()->Register("", this);

    // --- ------------------ --- //
    // --- GET ACCESS TO DATA --- //
    // --- ------------------ --- //

    // === get access to frssci pos cal data ===//
    fFrsSci_PosCal = dynamic_cast<TClonesArray*>(mgr->GetObject("FrsSciPosCalData"));
    if (!fFrsSci_PosCal)
    {
        LOG(info) << "R3BOnlineSpectraFrsSciVsMusli::Init() :: FrsSciPosCal not found";
    }
    // === get access to frssci tof cal data ===//
    fFrsSci_TofCal = dynamic_cast<TClonesArray*>(mgr->GetObject("FrsSciTofCalData"));
    if (!fFrsSci_TofCal)
    {
        LOG(info) << "R3BOnlineSpectraFrsSciVsMusli::Init() :: FrsSciTofCal not found";
    }

    // === get access to musli cal data ===//
    fMusli_Map = dynamic_cast<TClonesArray*>(mgr->GetObject("MusliMappedData"));
    if (!fMusli_Map)
    {
        LOG(info) << "R3BOnlineSpectraFrsSciVsMusli::Init() :: MusliMap not found";
    }

    // === get access to musli cal data ===//
    fMusli_Cal = dynamic_cast<TClonesArray*>(mgr->GetObject("MusliCalData"));
    if (!fMusli_Cal)
    {
        LOG(info) << "R3BOnlineSpectraFrsSciVsMusli::Init() :: MusliCal not found";
    }
    // === get access to musli hit data ===//
    fMusli_Hit = dynamic_cast<TClonesArray*>(mgr->GetObject("MusliHitData"));
    if (!fMusli_Hit)
    {
        LOG(info) << "R3BOnlineSpectraFrsSciVsMusli::Init() :: MusliHit not found";
    }

    // --- ---------------------------- --- //
    // --- DECLARATION TCANVAS + HISTOS --- //
    // --- ---------------------------- --- //

    if (fFrsSci_PosCal && fMusli_Cal && fIdCaveC > 0)
    {
        sprintf(Name1, "SciVsMusli_PosRaw");
        cMap_Pos = new TCanvas(Name1, Name1, 10, 10, 800, 700);
        cMap_Pos->Divide(2, 4);
        fh2_Map_PosRawVsDT = new TH2D*[8];
        sprintf(Name1, "SciVsMusli_PosCal");
        cCal_Pos = new TCanvas(Name1, Name1, 10, 10, 800, 700);
        cCal_Pos->Divide(2, 4);
        fh2_Cal_PosCalVsDT = new TH2D*[8];

        for (UShort_t i = 0; i < 8; i++)
        {
            sprintf(Name1, "PosRaw_vs_DTraw%i", i + 1);
            fh2_Map_PosRawVsDT[i] = new TH2D(Name1, Name1, 2000, -500, 35000, 1000, -5, 5);
            fh2_Map_PosRawVsDT[i]->GetXaxis()->SetTitle(Form("Drift time pair %i", i + 1));
            fh2_Map_PosRawVsDT[i]->GetYaxis()->SetTitle(Form("PosRaw at Cave C FrsSci%i", fIdCaveC));
            fh2_Map_PosRawVsDT[i]->GetXaxis()->CenterTitle(true);
            fh2_Map_PosRawVsDT[i]->GetYaxis()->CenterTitle(true);
            fh2_Map_PosRawVsDT[i]->GetXaxis()->SetLabelSize(0.05);
            fh2_Map_PosRawVsDT[i]->GetXaxis()->SetTitleSize(0.05);
            fh2_Map_PosRawVsDT[i]->GetYaxis()->SetLabelSize(0.05);
            fh2_Map_PosRawVsDT[i]->GetYaxis()->SetTitleSize(0.05);
            cMap_Pos->cd(i + 1);
            fh2_Map_PosRawVsDT[i]->Draw("col");

            sprintf(Name1, "PosCal_vs_DTcal%i", i + 1);
            fh2_Cal_PosCalVsDT[i] = new TH2D(Name1, Name1, 1000, -200, 200, 500, -100, 100);
            fh2_Cal_PosCalVsDT[i]->GetXaxis()->SetTitle(Form("Drift time pair %i", i + 1));
            fh2_Cal_PosCalVsDT[i]->GetYaxis()->SetTitle(Form("PosCal at Cave C FrsSci%i", fIdCaveC));
            fh2_Cal_PosCalVsDT[i]->GetXaxis()->CenterTitle(true);
            fh2_Cal_PosCalVsDT[i]->GetYaxis()->CenterTitle(true);
            fh2_Cal_PosCalVsDT[i]->GetXaxis()->SetLabelSize(0.05);
            fh2_Cal_PosCalVsDT[i]->GetXaxis()->SetTitleSize(0.05);
            fh2_Cal_PosCalVsDT[i]->GetYaxis()->SetLabelSize(0.05);
            fh2_Cal_PosCalVsDT[i]->GetYaxis()->SetTitleSize(0.05);
            cCal_Pos->cd(i + 1);
            fh2_Cal_PosCalVsDT[i]->Draw("col");
        }

        if (fFrsSci_TofCal && fIdS2 > 0)
        {
            sprintf(Name1, "SciVsMusliCal_EvsAoQ");
            cCal_EvsAoQ = new TCanvas(Name1, Name1, 10, 10, 800, 700);
            cCal_EvsAoQ->Divide(5, 3);
            fh2_Cal_EvsAoQ = new TH2D*[15];

            for (UShort_t i = 0; i < 15; i++)
            {
                sprintf(Name1, "E_vs_AoQ_sig%i", i + 1);
                fh2_Cal_EvsAoQ[i] = new TH2D(Name1, Name1, 3000, 0.5, 3.5, 1500, 0, 60000);
                fh2_Cal_EvsAoQ[i]->GetXaxis()->SetTitle("A/Q, from S2 to cave C");
                fh2_Cal_EvsAoQ[i]->GetYaxis()->SetTitle(Form("E@cal for signal%i", i + 1));
                fh2_Cal_EvsAoQ[i]->GetXaxis()->CenterTitle(true);
                fh2_Cal_EvsAoQ[i]->GetYaxis()->CenterTitle(true);
                fh2_Cal_EvsAoQ[i]->GetXaxis()->SetLabelSize(0.05);
                fh2_Cal_EvsAoQ[i]->GetXaxis()->SetTitleSize(0.05);
                fh2_Cal_EvsAoQ[i]->GetYaxis()->SetLabelSize(0.05);
                fh2_Cal_EvsAoQ[i]->GetYaxis()->SetTitleSize(0.05);
                cCal_EvsAoQ->cd(i + 1);
                fh2_Cal_EvsAoQ[i]->Draw("col");
            }
            if (fMusli_Hit && fIdCaveC > fIdS2)
            {
                sprintf(Name1, "SciVsMusliHit_EvsTofRaw");
                cHit_EvsTofRaw = new TCanvas(Name1, Name1, 10, 10, 800, 700);
                cHit_EvsTofRaw->Divide(2, 2);
                fh2_Hit_EvsTofRaw = new TH2D*[4];
                sprintf(Name1, "SciVsMusliHit_EvsBeta");
                cHit_EvsBeta = new TCanvas(Name1, Name1, 10, 10, 800, 700);
                cHit_EvsBeta->Divide(2, 2);
                fh2_Hit_EvsBeta = new TH2D*[4];
                sprintf(Name1, "SciVsMusliHit_EvsAoQ");
                cHit_EvsAoQ = new TCanvas(Name1, Name1, 10, 10, 800, 700);
                cHit_EvsAoQ->Divide(2, 2);
                fh2_Hit_EvsAoQ = new TH2D*[4];
                sprintf(Name1, "SciVsMusliHit_ZvsAoQ");
                cHit_ZvsAoQ = new TCanvas(Name1, Name1, 10, 10, 800, 700);
                cHit_ZvsAoQ->Divide(2, 2);
                fh2_Hit_ZvsAoQ = new TH2D*[4];
                for (UShort_t i = 0; i < 4; i++)
                {
                    sprintf(Name1, "Eave_vs_TofRaw_typ%i", i + 1);
                    fh2_Hit_EvsTofRaw[i] = new TH2D(Name1, Name1, 1200, 1380, 1500, 2000, 0, 40000);
                    fh2_Hit_EvsTofRaw[i]->GetXaxis()->SetTitle("TofRaw, from S2 to cave C");
                    fh2_Hit_EvsTofRaw[i]->GetYaxis()->SetTitle(
                        Form("<E>@hit for type%i (data per %i anodes)", i + 1, (int)pow(2, i + 1)));
                    fh2_Hit_EvsTofRaw[i]->GetXaxis()->CenterTitle(true);
                    fh2_Hit_EvsTofRaw[i]->GetYaxis()->CenterTitle(true);
                    fh2_Hit_EvsTofRaw[i]->GetXaxis()->SetLabelSize(0.05);
                    fh2_Hit_EvsTofRaw[i]->GetXaxis()->SetTitleSize(0.05);
                    fh2_Hit_EvsTofRaw[i]->GetYaxis()->SetLabelSize(0.05);
                    fh2_Hit_EvsTofRaw[i]->GetYaxis()->SetTitleSize(0.05);
                    cHit_EvsTofRaw->cd(i + 1);
                    fh2_Hit_EvsTofRaw[i]->Draw("col");

                    sprintf(Name1, "Eave_vs_Beta_typ%i", i + 1);
                    fh2_Hit_EvsBeta[i] = new TH2D(Name1, Name1, 1200, 0.69, 0.81, 2000, 0, 40000);
                    fh2_Hit_EvsBeta[i]->GetXaxis()->SetTitle("Beta, from S2 to cave C");
                    fh2_Hit_EvsBeta[i]->GetYaxis()->SetTitle(
                        Form("<E>@hit for type%i (data per %i anodes)", i + 1, (int)pow(2, i + 1)));
                    fh2_Hit_EvsBeta[i]->GetXaxis()->CenterTitle(true);
                    fh2_Hit_EvsBeta[i]->GetYaxis()->CenterTitle(true);
                    fh2_Hit_EvsBeta[i]->GetXaxis()->SetLabelSize(0.05);
                    fh2_Hit_EvsBeta[i]->GetXaxis()->SetTitleSize(0.05);
                    fh2_Hit_EvsBeta[i]->GetYaxis()->SetLabelSize(0.05);
                    fh2_Hit_EvsBeta[i]->GetYaxis()->SetTitleSize(0.05);
                    cHit_EvsBeta->cd(i + 1);
                    fh2_Hit_EvsBeta[i]->Draw("col");

                    sprintf(Name1, "Eave_vs_AoQ_typ%i", i + 1);
                    fh2_Hit_EvsAoQ[i] = new TH2D(Name1, Name1, 3000, 0.5, 3.5, 2000, 0, 40000);
                    fh2_Hit_EvsAoQ[i]->GetXaxis()->SetTitle("A/Q, from S2 to cave C");
                    fh2_Hit_EvsAoQ[i]->GetYaxis()->SetTitle(
                        Form("<E>@hit for type%i (data per %i anodes)", i + 1, (int)pow(2, i + 1)));
                    fh2_Hit_EvsAoQ[i]->GetXaxis()->CenterTitle(true);
                    fh2_Hit_EvsAoQ[i]->GetYaxis()->CenterTitle(true);
                    fh2_Hit_EvsAoQ[i]->GetXaxis()->SetLabelSize(0.05);
                    fh2_Hit_EvsAoQ[i]->GetXaxis()->SetTitleSize(0.05);
                    fh2_Hit_EvsAoQ[i]->GetYaxis()->SetLabelSize(0.05);
                    fh2_Hit_EvsAoQ[i]->GetYaxis()->SetTitleSize(0.05);
                    cHit_EvsAoQ->cd(i + 1);
                    fh2_Hit_EvsAoQ[i]->Draw("col");

                    sprintf(Name1, "Z_vs_AoQ_typ%i", i + 1);
                    fh2_Hit_ZvsAoQ[i] = new TH2D(Name1, Name1, 3000, 0.5, 3.5, 1000, -0.5, 9.5);
                    fh2_Hit_ZvsAoQ[i]->GetXaxis()->SetTitle("A/Q, from S2 to cave C");
                    fh2_Hit_ZvsAoQ[i]->GetYaxis()->SetTitle(
                        Form("Z@hit for type%i (data per %i anodes)", i + 1, (int)pow(2, i + 1)));
                    fh2_Hit_ZvsAoQ[i]->GetXaxis()->CenterTitle(true);
                    fh2_Hit_ZvsAoQ[i]->GetYaxis()->CenterTitle(true);
                    fh2_Hit_ZvsAoQ[i]->GetXaxis()->SetLabelSize(0.05);
                    fh2_Hit_ZvsAoQ[i]->GetXaxis()->SetTitleSize(0.05);
                    fh2_Hit_ZvsAoQ[i]->GetYaxis()->SetLabelSize(0.05);
                    fh2_Hit_ZvsAoQ[i]->GetYaxis()->SetTitleSize(0.05);
                    cHit_ZvsAoQ->cd(i + 1);
                    fh2_Hit_ZvsAoQ[i]->Draw("col");
                }
            }
        }
    }

    // --- --------------- --- //
    // --- MAIN FOLDER-Sci --- //
    // --- --------------- --- //
    TFolder* mainfol = new TFolder("FrsSciVsMusli", "FrsSci vs Musli info");

    if (fFrsSci_PosCal && fMusli_Map && fMusli_Cal && fIdCaveC > 0)
    {
        mainfol->Add(cMap_Pos);
        mainfol->Add(cCal_Pos);
        if (fFrsSci_TofCal && fIdS2 > 0)
        {
            mainfol->Add(cCal_EvsAoQ);
            if (fMusli_Hit && fIdCaveC > fIdS2)
            {
                mainfol->Add(cHit_EvsTofRaw);
                mainfol->Add(cHit_EvsBeta);
                mainfol->Add(cHit_EvsAoQ);
                mainfol->Add(cHit_ZvsAoQ);
            }
        }
    }
    run->AddObject(mainfol);

    // Register command to reset histograms
    run->GetHttpServer()->RegisterCommand("Reset_FRSSCIvsMUSLI_HIST", Form("/Objects/%s/->Reset_Histo()", GetName()));

    return kSUCCESS;
}

void R3BOnlineSpectraFrsSciVsMusli::Reset_Histo()
{

    LOG(info) << "R3BOnlineSpectraFrsSciVsMusli::Reset_Histo";

    if (fFrsSci_PosCal && fMusli_Cal && fIdCaveC > 0)
    {
        for (UShort_t i = 0; i < 8; i++)
        {
            fh2_Map_PosRawVsDT[i]->Reset();
            fh2_Cal_PosCalVsDT[i]->Reset();
        }
        if (fFrsSci_TofCal && fIdS2 > 0)
        {
            for (UShort_t i = 0; i < 15; i++)
            {
                fh2_Cal_EvsAoQ[i]->Reset();
            }

            if (fMusli_Hit && fIdCaveC > fIdS2)
            {
                for (UShort_t i = 0; i < 4; i++)
                {
                    fh2_Hit_EvsTofRaw[i]->Reset();
                    fh2_Hit_EvsBeta[i]->Reset();
                    fh2_Hit_EvsAoQ[i]->Reset();
                    fh2_Hit_ZvsAoQ[i]->Reset();
                }
            }
        }
    }
}

void R3BOnlineSpectraFrsSciVsMusli::Exec(Option_t* option)
{

    FairRootManager* mgr = FairRootManager::Instance();
    if (NULL == mgr)
    {
        LOG(fatal) << "R3BOnlineSpectraFrsSciVsMusli::Exec FairRootManager not found";
    }
    UInt_t nHits;

    UInt_t iSig;
    Double_t iRawPosCaveC = -1000.;
    Double_t iCalPosCaveC = -1000.;
    Double_t iAoQ = -1.;
    Double_t iRawTof = -1.;
    Double_t iCalTof = -1.;
    Double_t iBeta = -1.;
    UInt_t iTyp;

    if (fFrsSci_PosCal && fMusli_Cal && fIdCaveC > 0)
    {
        // --- read frssci poscal data and extract raw and cal pos at cave C --- //
        if (fFrsSci_PosCal->GetEntriesFast() > 0)
        {
            nHits = fFrsSci_PosCal->GetEntriesFast();
            for (UInt_t ihit = 0; ihit < nHits; ihit++)
            {
                R3BFrsSciPosCalData* hitposcal = dynamic_cast<R3BFrsSciPosCalData*>(fFrsSci_PosCal->At(ihit));
                if (!hitposcal)
                    continue;
                if (hitposcal->GetDetector() == fIdCaveC)
                {
                    iRawPosCaveC = hitposcal->GetRawPosNs();
                    iCalPosCaveC = hitposcal->GetCalPosMm();
                }
            }
        }

        // --- read musli mapped data and fill histo at MusliMapped level--- //
        if (fMusli_Cal->GetEntriesFast() > 0)
        {
            nHits = fMusli_Map->GetEntriesFast();
            UInt_t Tref = 0;
            for (UInt_t ihit = 0; ihit < nHits; ihit++)
            {
                R3BMusliMappedData* hitmap = dynamic_cast<R3BMusliMappedData*>(fMusli_Map->At(ihit));
                if (!hitmap)
                    continue;
                iSig = hitmap->GetSignal();
                if (iSig != 17)
                    continue;
                else
                    Tref = hitmap->GetTime();
            }
            if (Tref > 0)
            {
                for (UInt_t ihit = 0; ihit < nHits; ihit++)
                {
                    R3BMusliMappedData* hitmap = dynamic_cast<R3BMusliMappedData*>(fMusli_Map->At(ihit));
                    if (!hitmap)
                        continue;
                    iSig = hitmap->GetSignal();
                    if (iSig > 8)
                        continue;
                    if (iRawPosCaveC != -1000.)
                        fh2_Map_PosRawVsDT[iSig - 1]->Fill(hitmap->GetTime() - Tref, iRawPosCaveC);
                }
            }
        }

        if (fFrsSci_TofCal)
        {
            // --- read tofcal data and extract AoQ for S2-caveC tof--- //
            if (fFrsSci_TofCal->GetEntriesFast() > 0)
            {
                nHits = fFrsSci_TofCal->GetEntriesFast();
                for (UInt_t ihit = 0; ihit < nHits; ihit++)
                {
                    R3BFrsSciTofCalData* hittofcal = dynamic_cast<R3BFrsSciTofCalData*>(fFrsSci_TofCal->At(ihit));
                    if (!hittofcal)
                        continue;
                    if (hittofcal->GetDetIdSta() == fIdS2 && hittofcal->GetDetIdSto() == fIdCaveC)
                    {
                        iRawTof = hittofcal->GetRawTofNs();
                        iCalTof = hittofcal->GetCalTofNs();
                        iBeta = hittofcal->GetBeta();
                        iAoQ = hittofcal->GetAoQ();
                    }
                }
            }
        }

        // --- read musli cal data and fill histo at MusliCal level--- //
        if (fMusli_Cal->GetEntriesFast() > 0)
        {
            nHits = fMusli_Cal->GetEntriesFast();
            for (UInt_t ihit = 0; ihit < nHits; ihit++)
            {
                R3BMusliCalData* hitcal = dynamic_cast<R3BMusliCalData*>(fMusli_Cal->At(ihit));
                if (!hitcal)
                    continue;
                iSig = hitcal->GetSignal();
                if (iSig > 15)
                    continue;
                if (iAoQ > 0)
                    fh2_Cal_EvsAoQ[iSig - 1]->Fill(iAoQ, hitcal->GetEnergy());
                if (iSig > 8)
                    continue;
                if (iCalPosCaveC != -1000.)
                    fh2_Cal_PosCalVsDT[iSig - 1]->Fill(hitcal->GetDT(), iCalPosCaveC);
            }
        }

        // --- read musli hit data and fill histo at MusliHit level--- //
        if (fMusli_Hit)
        {
            if (fMusli_Hit->GetEntriesFast() > 0)
            {
                nHits = fMusli_Hit->GetEntriesFast();
                for (UInt_t ihit = 0; ihit < nHits; ihit++)
                {
                    R3BMusliHitData* hithit = dynamic_cast<R3BMusliHitData*>(fMusli_Hit->At(ihit));
                    if (!hithit)
                        continue;
                    iTyp = hithit->GetType() - 1;
                    if (iAoQ > 0)
                    {
                        fh2_Hit_EvsTofRaw[iTyp]->Fill(iRawTof, hithit->GetEave());
                        fh2_Hit_EvsBeta[iTyp]->Fill(iRawTof, hithit->GetEave());
                        fh2_Hit_EvsAoQ[iTyp]->Fill(iAoQ, hithit->GetEave());
                        fh2_Hit_ZvsAoQ[iTyp]->Fill(iAoQ, hithit->GetZcharge());
                    }
                }
            }
        }
        fNEvents++;
    } // end of if fFrsSci_PosCal && fMusli_Cal && fIdCaveC>0
}

// -----   Public method Finish   -----------------------------------------------
void R3BOnlineSpectraFrsSciVsMusli::FinishEvent()
{
    if (fFrsSci_PosCal)
    {
        fFrsSci_PosCal->Clear();
    }
    if (fFrsSci_TofCal)
    {
        fFrsSci_TofCal->Clear();
    }
    if (fMusli_Map)
    {
        fMusli_Map->Clear();
    }
    if (fMusli_Cal)
    {
        fMusli_Cal->Clear();
    }
    if (fMusli_Hit)
    {
        fMusli_Hit->Clear();
    }
}

void R3BOnlineSpectraFrsSciVsMusli::FinishTask()
{

    if (fFrsSci_PosCal && fMusli_Cal && fIdCaveC > 0)
    {
        for (UShort_t i = 0; i < 8; i++)
        {
            fh2_Map_PosRawVsDT[i]->Write();
            fh2_Cal_PosCalVsDT[i]->Write();
        }

        if (fFrsSci_TofCal && fIdS2 > 0)
        {
            for (UShort_t i = 0; i < 15; i++)
            {
                fh2_Cal_EvsAoQ[i]->Write();
            }
            if (fMusli_Hit)
            {
                for (UShort_t i = 0; i < 4; i++)
                {
                    fh2_Hit_EvsTofRaw[i]->Write();
                    fh2_Hit_EvsBeta[i]->Write();
                    fh2_Hit_EvsAoQ[i]->Write();
                    fh2_Hit_ZvsAoQ[i]->Write();
                }
            }
        }
    }
}
ClassImp(R3BOnlineSpectraFrsSciVsMusli);
