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
// -----                  R3BPspxOnlineSpectra            -----
// -----          Created Feb 23th 2019 by J.L. Rodriguez -----
// -----              Modified Dec 2019 by M. Holl        -----
// ------------------------------------------------------------

/*
 * This task should fill histograms with detector variables which allow
 * to test the detectors online
 *
 */

#include "R3BPspxOnlineSpectra.h"
#include "R3BPspxCalData.h"
#include "R3BPspxHitData.h"
#include "R3BPspxMappedData.h"
#include "R3BPspxPrecalData.h"

#include "R3BEventHeader.h"

#include "FairLogger.h"
#include "FairRootManager.h"
#include "FairRunAna.h"
#include "FairRunOnline.h"
#include "FairRuntimeDb.h"
#include "TCanvas.h"
#include "TFolder.h"
#include "TGaxis.h"
#include "TH1F.h"
#include "TH2F.h"

#include "THttpServer.h"

#include "TClonesArray.h"
#include "TMath.h"
#include "TVector3.h"
#include <array>
#include <cstdlib>
#include <ctime>
#include <fstream>
#include <iostream>
#include <sstream>

using namespace std;

R3BPspxOnlineSpectra::R3BPspxOnlineSpectra()
    : FairTask("PspOnlineSpectra", 1)
    , fMappedItemsPspx()
    , fPrecalItemsPspx()
    , fCalItemsPspx()
    , fHitItemsPspx()
    , fTrigger(-1)
    , fNEvents(0)
{
}

R3BPspxOnlineSpectra::R3BPspxOnlineSpectra(const char* name, Int_t iVerbose)
    : FairTask(name, iVerbose)
    , fMappedItemsPspx()
    , fPrecalItemsPspx()
    , fCalItemsPspx()
    , fHitItemsPspx()
    , fTrigger(-1)
    , fNEvents(0)
{
}

R3BPspxOnlineSpectra::~R3BPspxOnlineSpectra() {}

InitStatus R3BPspxOnlineSpectra::Init()
{
    LOG(info) << "R3BPspxOnlineSpectra::Init ";
    // try to get a handle on the EventHeader. EventHeader may not be
    // present though and hence may be null. Take care when using.

    FairRootManager* mgr = FairRootManager::Instance();
    if (NULL == mgr)
        LOG(fatal) << "FairRootManager not found";
    header = dynamic_cast<R3BEventHeader*>(mgr->GetObject("R3BEventHeader"));

    FairRunOnline* run = FairRunOnline::Instance();
    run->GetHttpServer()->Register("", this);

    // create histograms of all detectors
    //
    // PSPX detector
    // get access to data
    const char c_xy[2] = { 'x', 'y' }; // orientation of detector face
    // LOG(info) << "TEST ---- " << PSPX;

    for (Int_t d = 0; d < PSPX; d++)
    {
        for (Int_t f = 0; f < 2; f++)
        {
            fMappedItemsPspx.push_back(
                dynamic_cast<TClonesArray*>(mgr->GetObject(Form("Pspx%d_%cMapped", d + 1, c_xy[f]))));
            fPrecalItemsPspx.push_back(
                dynamic_cast<TClonesArray*>(mgr->GetObject(Form("Pspx%d_%cPrecal", d + 1, c_xy[f]))));
            fCalItemsPspx.push_back(dynamic_cast<TClonesArray*>(mgr->GetObject(Form("Pspx%d_%cCal", d + 1, c_xy[f]))));
            fHitItemsPspx.push_back(dynamic_cast<TClonesArray*>(mgr->GetObject(Form("Pspx%d_%cHit", d + 1, c_xy[f]))));
        }

        if (fMappedItemsPspx[0] == NULL)
        {
            printf("R3BPspxOnlineSpectra: Couldn't get handle on PSPX mapped items\n");
            return kFATAL;
        }
    }
    Int_t mappedSize = fMappedItemsPspx.size();
    Int_t precalSize = fPrecalItemsPspx.size();
    Int_t calSize = fCalItemsPspx.size();
    // LOG(info) << "TEST 2 ---- " << fMappedItemsPspx.size();
    Int_t hitSize = fHitItemsPspx.size();

    Int_t Emax = 6000;

    std::string xy[2] = { "x", "y" };
    // std::string ep[2] = {"Energy","Position"};
    std::string histoName, histoTitle;
    fh_pspx_multiplicity.resize(mappedSize);
    fh_pspx_strip_1.resize(mappedSize);
    fh_pspx_energy_strip_1.resize(mappedSize);
    fh_pspx_strip_2.resize(mappedSize);
    fh_pspx_energy_strip_2.resize(mappedSize);

    fh_pspx_cal_strip_frontback.resize(calSize);
    fh_pspx_cal_pos_frontback.resize(calSize);
    fh_pspx_cal_energy_frontback.resize(calSize);

    fh_pspx_cal_energyfront.resize(calSize);
    fh_pspx_cal_energyback.resize(calSize);

    fh_pspx_cal_strip.resize(precalSize);

    // Create histograms
    // Mapped level

    for (UInt_t i = 0; i < mappedSize; i++)
    {
        histoName = "pspx_" + std::to_string((i / 2) + 1) + "_" + xy[i % 2] + "_multiplicity";
        histoTitle = "Pspx " + std::to_string((i / 2) + 1) + ": " + xy[i % 2] + " Multiplicity;Multiplicity;Counts";
        fh_pspx_multiplicity[i] = new TH1F(histoName.data(), histoTitle.data(), 10, 0, 10);

        histoName = "pspx_" + std::to_string((i / 2) + 1) + "_" + xy[i % 2] + "_strip_1";
        histoTitle = "Pspx " + std::to_string((i / 2) + 1) + ": " + xy[i % 2] + " Strip Side 1;Strip;Counts";
        fh_pspx_strip_1[i] = new TH1F(histoName.data(), histoTitle.data(), 32, 1, 33);

        histoName = "pspx_" + std::to_string((i / 2) + 1) + "_" + xy[i % 2] + "_energy_strip_1";
        histoTitle = "Pspx " + std::to_string((i / 2) + 1) + ": " + xy[i % 2] + " Energy vs Strip Side 1;Strip;Energy";
        fh_pspx_energy_strip_1[i] = new TH2F(histoName.data(), histoTitle.data(), 32, 1, 33, 1000, 0, Emax);

        histoName = "pspx_" + std::to_string((i / 2) + 1) + "_" + xy[i % 2] + "_strip_2";
        histoTitle = "Pspx " + std::to_string((i / 2) + 1) + ": " + xy[i % 2] + " Strip Side 2;Strip;Counts";
        fh_pspx_strip_2[i] = new TH1F(histoName.data(), histoTitle.data(), 32, 1, 33);

        histoName = "pspx_" + std::to_string((i / 2) + 1) + "_" + xy[i % 2] + "_energy_strip_2";
        histoTitle = "Pspx " + std::to_string((i / 2) + 1) + ": " + xy[i % 2] + " Energy vs Strip Side 2;Strip;Energy";
        fh_pspx_energy_strip_2[i] = new TH2F(histoName.data(), histoTitle.data(), 32, 1, 33, 1000, 0, Emax);
    }

    // Cal level
    for (UInt_t i = 0; i < calSize / 2; i++)
    {
        histoName = "pspx_" + std::to_string(i + 1) + "_cal_strip_frontback";
        histoTitle = "Pspx " + std::to_string(i + 1) + ": " + " Cal Strip;XStrip;YStrip";
        fh_pspx_cal_strip_frontback[i] = new TH2F(histoName.data(), histoTitle.data(), 32, 1, 33, 32, 1, 33);
        histoName = "pspx_" + std::to_string(i + 1) + "_cal_pos_frontback";
        histoTitle = "Pspx " + std::to_string(i + 1) + ": " + " Cal Pos;XPos;YPos";
        fh_pspx_cal_pos_frontback[i] = new TH2F(histoName.data(), histoTitle.data(), 200, -1, 1, 200, -1, 1);
        histoName = "pspx_" + std::to_string(i + 1) + "_cal_energy_frontback";
        histoTitle = "Pspx " + std::to_string(i + 1) + ": " + " Cal Energy;XEnergy;YEnergy";
        fh_pspx_cal_energy_frontback[i] = new TH2F(histoName.data(), histoTitle.data(), 1e3, 0, Emax, 1e3, 0, Emax);

        histoName = "pspx_" + std::to_string(i + 1) + "_cal_energyfront";
        histoTitle = "Pspx " + std::to_string(i + 1) + ": " + " Cal X Energy;XEnergy";
        fh_pspx_cal_energyfront[i] = new TH1F(histoName.data(), histoTitle.data(), 1e3, 0, Emax);
        histoName = "pspx_" + std::to_string(i + 1) + "_cal_energyback";
        histoTitle = "Pspx " + std::to_string(i + 1) + ": " + " Cal Y Energy;YEnergy";
        fh_pspx_cal_energyback[i] = new TH1F(histoName.data(), histoTitle.data(), 1e3, 0, Emax);
    }
    for (UInt_t i = 0; i < precalSize; i++)
    {
        histoName = "pspx_" + std::to_string((i / 2) + 1) + "_cal_strip_" + xy[i % 2];
        histoTitle = "Pspx " + std::to_string((i / 2) + 1) + ": " + xy[i % 2] + " Cal Strip;XStrip;";
        fh_pspx_cal_strip[i] = new TH1F(histoName.data(), histoTitle.data(), 32, 1, 33);
    }

    // Create canvasses
    TFolder* mainfolPspx = new TFolder("PSPX", "PSPX info");
    TCanvas* cPspx_strips = new TCanvas("Pspx_strips", "Pspx Strips", 10, 10, 1100, 1000);
    cPspx_strips->Divide(4, PSPX);
    for (UInt_t i = 0; i < mappedSize; i++)
    {
        cPspx_strips->cd(2 * i + 1);
        fh_pspx_strip_1[i]->Draw();
        cPspx_strips->cd(2 * i + 2);
        fh_pspx_strip_2[i]->Draw();
    }
    mainfolPspx->Add(cPspx_strips);

    TCanvas* cPspx_multiplicity = new TCanvas("Pspx_multiplicity", "Pspx Multiplicity", 10, 10, 1100, 1000);
    cPspx_multiplicity->Divide(mappedSize, 1);
    for (UInt_t i = 0; i < mappedSize; i++)
    {
        cPspx_multiplicity->cd(i + 1);
        fh_pspx_multiplicity[i]->Draw();
    }
    mainfolPspx->Add(cPspx_multiplicity);

    TCanvas* cPspx_energy_strips = new TCanvas("Pspx_energy_vs_strip", "Pspx Energy vs Strips", 10, 10, 1100, 1000);
    cPspx_energy_strips->Divide(4, PSPX);
    for (UInt_t i = 0; i < mappedSize; i++)
    {
        cPspx_energy_strips->cd(2 * i + 1);
        fh_pspx_energy_strip_1[i]->Draw("colz");
        cPspx_energy_strips->cd(2 * i + 2);
        fh_pspx_energy_strip_2[i]->Draw("colz");
    }
    mainfolPspx->Add(cPspx_energy_strips);

    TCanvas* cPspx_cal = new TCanvas("Pspx_cal", "Pspx Cal Level front vs back", 10, 10, 1100, 1000);
    cPspx_cal->Divide(3, PSPX);
    for (UInt_t i = 0; i < calSize / 2; i++)
    {
        cPspx_cal->cd(i + 1); // i*2
        fh_pspx_cal_strip_frontback[i]->Draw("colz");
        cPspx_cal->cd(i + 1 + PSPX); // i*2
        fh_pspx_cal_pos_frontback[i]->Draw("colz");
        cPspx_cal->cd(i + 1 + 2 * PSPX); // i*2
        fh_pspx_cal_energy_frontback[i]->Draw("colz");
    }
    mainfolPspx->Add(cPspx_cal);

    TCanvas* cPspx_precal_strips = new TCanvas("Pspx_cal_strips", "Pspx Cal Strip x", 10, 10, 1100, 1000);
    cPspx_precal_strips->Divide(2, PSPX);
    for (UInt_t i = 0; i < 2 * PSPX; i++)
    {
        cPspx_precal_strips->cd(i + 1);
        fh_pspx_cal_strip[i]->Draw();
    }
    mainfolPspx->Add(cPspx_precal_strips);

    TCanvas* cPspx_energy = new TCanvas("Pspx_energy", "Pspx_energy", 10, 10, 1100, 1000);
    cPspx_energy->Divide(PSPX, 2);
    for (UInt_t i = 0; i < PSPX; i++)
    {
        cPspx_energy->cd(i + 1);
        fh_pspx_cal_energyfront[i]->Draw();
        cPspx_energy->cd(i + 1 + PSPX);
        fh_pspx_cal_energyback[i]->Draw();
    }
    mainfolPspx->Add(cPspx_energy);

    run->AddObject(mainfolPspx);

    run->GetHttpServer()->RegisterCommand("Reset_PSPX", Form("/Objects/%s/->Reset_PSPX_Histo()", GetName()));

    // -------------------------------------------------------------------------
    // LOG(info) << PSPX ;
    LOG(info) << "END of INIT ";
    return kSUCCESS;
}

void R3BPspxOnlineSpectra::Reset_PSPX_Histo()
{
    for (UInt_t i = 0; i < 2 * PSPX; i++)
    {
        fh_pspx_multiplicity[i]->Reset();
        fh_pspx_strip_1[i]->Reset();
        fh_pspx_strip_2[i]->Reset();
        fh_pspx_energy_strip_1[i]->Reset();
        fh_pspx_energy_strip_2[i]->Reset();

        fh_pspx_cal_strip[i]->Reset();
    }

    for (UInt_t i = 0; i < 2 * PSPX; i++)
    {
        fh_pspx_cal_strip_frontback[i]->Reset();
        fh_pspx_cal_pos_frontback[i]->Reset();
        fh_pspx_cal_energy_frontback[i]->Reset();
        fh_pspx_cal_energyfront[i]->Reset();
        fh_pspx_cal_energyback[i]->Reset();
    }
}

void R3BPspxOnlineSpectra::Exec(Option_t* option)
{
    FairRootManager* mgr = FairRootManager::Instance();
    if (NULL == mgr)
        LOG(fatal) << "R3BPspxOnlineSpectra::Exec FairRootManager not found";
    // check for requested trigger (Todo: should be done globablly / somewhere else)
    if ((fTrigger >= 0) && (header) && (header->GetTrigger() != fTrigger))
        return;

    //----------------------------------------------------------------------
    // PSPX
    //----------------------------------------------------------------------
    for (UInt_t d = 0; d < 2 * PSPX; d++)
    {
        Int_t nHits = fMappedItemsPspx[d]->GetEntriesFast();
        fh_pspx_multiplicity[d]->Fill(nHits);
        for (Int_t ihit = 0; ihit < nHits; ihit++)
        {
            R3BPspxMappedData* mappedData = dynamic_cast<R3BPspxMappedData*>(fMappedItemsPspx[d]->At(ihit));

            fh_pspx_strip_1[d]->Fill(mappedData->GetStrip1());
            fh_pspx_strip_2[d]->Fill(mappedData->GetStrip2());
            fh_pspx_energy_strip_1[d]->Fill(mappedData->GetStrip1(), mappedData->GetEnergy1());
            fh_pspx_energy_strip_2[d]->Fill(mappedData->GetStrip2(), mappedData->GetEnergy2());
        }
    }
    for (UInt_t d = 0; d < PSPX; d++)
    {
        Int_t nHits1 = fCalItemsPspx[2 * d]->GetEntriesFast();
        Int_t nHits2 = fCalItemsPspx[2 * d + 1]->GetEntriesFast();
        Int_t nHits = (nHits1 < nHits2) ? nHits1 : nHits2;

        for (Int_t ihit = 0; ihit < nHits; ihit++)
        {
            R3BPspxCalData* calData1 = dynamic_cast<R3BPspxCalData*>(fCalItemsPspx[2 * d]->At(ihit));
            R3BPspxCalData* calData2 = dynamic_cast<R3BPspxCalData*>(fCalItemsPspx[2 * d + 1]->At(ihit));
            fh_pspx_cal_strip_frontback[d]->Fill(calData1->GetStrip(), calData2->GetStrip());
            fh_pspx_cal_pos_frontback[d]->Fill(calData1->GetPos(), calData2->GetPos());
            fh_pspx_cal_energy_frontback[d]->Fill(calData1->GetEnergy(), calData2->GetEnergy());
            fh_pspx_cal_energyfront[d]->Fill(calData1->GetEnergy());
            fh_pspx_cal_energyback[d]->Fill(calData2->GetEnergy());
        }
    }
    for (UInt_t d = 0; d < 2 * PSPX; d++)
    {
        Int_t nHits_cal = fPrecalItemsPspx[d]->GetEntriesFast();
        for (Int_t ihit = 0; ihit < nHits_cal; ihit++)
        {
            R3BPspxPrecalData* precalData = dynamic_cast<R3BPspxPrecalData*>(fPrecalItemsPspx[d]->At(ihit));
            fh_pspx_cal_strip[d]->Fill(precalData->GetStrip());
        }
    }
    fNEvents += 1;
}

void R3BPspxOnlineSpectra::FinishEvent()
{
    fMappedItemsPspx.clear();
    fPrecalItemsPspx.clear();
    fCalItemsPspx.clear();
    fHitItemsPspx.clear();
}

void R3BPspxOnlineSpectra::FinishTask()
{

    LOG(info) << "Finish MappedPspx";

    // for (UInt_t i = 0; i < fMappedItemsPspx.size(); i++)
    for (UInt_t i = 0; i < 2 * PSPX; i++)
    {
        fh_pspx_multiplicity[i]->Write();
        fh_pspx_strip_1[i]->Write();
        fh_pspx_strip_2[i]->Write();
        fh_pspx_energy_strip_1[i]->Write();
        fh_pspx_energy_strip_2[i]->Write();

        fh_pspx_cal_strip[i]->Write();
    }
    for (UInt_t i = 0; i < PSPX; i++)
    {
        fh_pspx_cal_strip_frontback[i]->Write();
        fh_pspx_cal_pos_frontback[i]->Write();
        fh_pspx_cal_energy_frontback[i]->Write();
        fh_pspx_cal_energyfront[i]->Write();
        fh_pspx_cal_energyback[i]->Write();
    }
}

ClassImp(R3BPspxOnlineSpectra)
