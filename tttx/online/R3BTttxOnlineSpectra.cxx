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
// -----         Created 23/01/24 by M. Whitehead         -----
// -----           Fill TTT10 online histograms           -----
// ------------------------------------------------------------

/*
 * This task should fill histograms with TTTX online data
 */

#include "R3BTttxOnlineSpectra.h"
#include "R3BEventHeader.h"
#include "R3BLogger.h"
#include "R3BShared.h"
#include "R3BTttxCalData.h"
#include "R3BTttxHitData.h"
#include "R3BTttxMappedData.h"

#include "FairLogger.h"
#include "FairRootManager.h"
#include "FairRunAna.h"
#include "FairRunOnline.h"
#include "FairRuntimeDb.h"
#include "TCanvas.h"
#include "TClonesArray.h"
#include "TFolder.h"
#include "TH1F.h"
#include "TH2F.h"
#include "THttpServer.h"
#include "TMath.h"
#include "TRandom.h"
#include "TVector3.h"
#include <fstream>
#include <iostream>
#include <sstream>

R3BTttxOnlineSpectra::R3BTttxOnlineSpectra()
    : R3BTttxOnlineSpectra("TttxOnlineSpectra", 1)
{
}

R3BTttxOnlineSpectra::R3BTttxOnlineSpectra(const TString& name, Int_t iVerbose)
    : FairTask(name, iVerbose)
    , fMappedItemsTttx(NULL)
    , fNEvents(0)
    , fExpId(0)
    , fNbDets(2)
    , fNbStrips(32)
    , fNbTref(1)
    , fNbTrig(1)
    , fNbEntries(50000)
{
}

R3BTttxOnlineSpectra::~R3BTttxOnlineSpectra() { R3BLOG(debug1, ""); }

InitStatus R3BTttxOnlineSpectra::Init()
{
    R3BLOG(info, "");

    // try to get a handle on the EventHeader. EventHeader may not be
    // present though and hence may be null. Take care when using.

    FairRootManager* rootManager = FairRootManager::Instance();
    R3BLOG_IF(fatal, NULL == rootManager, "FairRootManager not found");

    FairRunOnline* run = FairRunOnline::Instance();
    run->GetHttpServer()->Register("", this);

    header = dynamic_cast<R3BEventHeader*>(rootManager->GetObject("EventHeader."));
    if (fExpId == 0) // Obtain global ExpId if it's not set locally.
    {
        fExpId = header->GetExpId();
        R3BLOG(info, "fExpId: " << fExpId);
    }

    // get access to mapped data of TTT10
    fMappedItemsTttx = dynamic_cast<TClonesArray*>(rootManager->GetObject("tttxMappedData"));
    if (!fMappedItemsTttx)
    {
        R3BLOG(fatal, "tttxMappedData not found");
        return kFATAL;
    }

    // get access to cal data of the TTTX
    fCalItemsTttx = dynamic_cast<TClonesArray*>(rootManager->GetObject("tttxCalData"));

    // get access to hit data of the TTTX
    fHitItemsTttx = dynamic_cast<TClonesArray*>(rootManager->GetObject("tttxHitData"));

    // Create histograms for detectors
    char Name1[255];
    char Name2[255];

    // Channel range for pileup/overflow
    Double_t binsE_pileup = 1000;
    Double_t minE_pileup = 0;
    Double_t maxE_pileup = 750000;

    // ADC Channel range
    Double_t binsE = 1000;
    Double_t minE = 0;
    Double_t maxE = 64000;

    // Time channel range
    Double_t binsT = 1000;
    Double_t minT = -32000;
    Double_t maxT = 10000;

    // Multiplicity range for strips
    Double_t binsM = 10;
    Double_t minM = -0.5;
    Double_t maxM = 9.5;

    // Energy range [keV]
    Double_t binsCalE = 2000;
    Double_t minCalE = 0;
    Double_t maxCalE = 10000;

    // Time range [ns]
    Double_t binsCalT = 500;
    Double_t minCalT = -3000;
    Double_t maxCalT = 3000;

    // Z of detector
    Double_t binsZ = 100;
    Double_t minZ = 0;
    Double_t maxZ = 10;

    // Pos of detector
    Double_t binsX = 32;
    Double_t minX = -50;
    Double_t maxX = 50;

    // MAIN FOLDER
    TFolder* mainfolTTTX = new TFolder("TTT10", "TTT10 info");

    // Folder for mapped data
    TFolder* mapfolTTTX = new TFolder("Map", "Map TTT10 info");

    // Folder for calibrated data
    TFolder* calfolTTTX = new TFolder("Cal", "Cal TTT10 info");

    // Folder for hit data
    TFolder* hitfolTTTX = new TFolder("Hit", "Hit TTT10 info");

    cMap = new TCanvas("Tttx_Mapped", "Tttx_Mapped", 10, 10, 800, 700);
    cMap->Divide(4, 3);
    mapfolTTTX->Add(cMap);
    cCal = new TCanvas("Tttx_Cal", "Tttx_Cal", 10, 10, 800, 700);
    cCal->Divide(4, 2);
    calfolTTTX->Add(cCal);
    cHit = new TCanvas("Tttx_Hit", "Tttx_Hit", 10, 10, 800, 700);
    cHit->Divide(4, 2);
    hitfolTTTX->Add(cHit);
    mainfolTTTX->Add(mapfolTTTX);
    mainfolTTTX->Add(calfolTTTX);
    mainfolTTTX->Add(hitfolTTTX);

    //==============================  Mapped Data  ========================================

    fh2_EnergyVsStrip_all.resize(fNbDets);
    fh2_EnergyVsStrip_recent.resize(fNbDets);
    fh2_TimeVsStrip.resize(fNbDets);
    fh1_Strip.resize(fNbDets);
    fh2_MultVsStrip.resize(fNbDets);

    for (Int_t i = 0; i < fNbDets; i++)
    {
        sprintf(Name1, "fh2_energy_vs_strip_recent_det_%d", i + 1);
        sprintf(Name2, "Mapped energy vs strip number for TTT10 Det: %d, last %d events", i + 1, fNbEntries);
        fh2_EnergyVsStrip_recent[i] =
            R3B::root_owned<TH2F>(Name1, Name2, fNbStrips, 1, fNbStrips + 1, binsE_pileup, minE_pileup, maxE_pileup);
        fh2_EnergyVsStrip_recent[i]->GetXaxis()->SetTitle("Strip number");
        fh2_EnergyVsStrip_recent[i]->GetYaxis()->SetTitle("Energy [channels]");
        fh2_EnergyVsStrip_recent[i]->GetYaxis()->SetTitleOffset(1.4);
        fh2_EnergyVsStrip_recent[i]->GetXaxis()->CenterTitle(true);
        cMap->cd(i + 1);
        fh2_EnergyVsStrip_recent[i]->Draw("colz");
    }
    for (Int_t i = 0; i < fNbDets; i++)
    {
        sprintf(Name1, "fh2_energy_vs_strip_all_det_%d", i + 1);
        sprintf(Name2, "Mapped energy vs strip number for TTT10 Det: %d, all events", i + 1);
        fh2_EnergyVsStrip_all[i] = R3B::root_owned<TH2F>(Name1, Name2, fNbStrips, 1, fNbStrips + 1, binsE, minE, maxE);
        fh2_EnergyVsStrip_all[i]->GetXaxis()->SetTitle("Strip number");
        fh2_EnergyVsStrip_all[i]->GetYaxis()->SetTitle("Energy [channels]");
        fh2_EnergyVsStrip_all[i]->GetYaxis()->SetTitleOffset(1.4);
        fh2_EnergyVsStrip_all[i]->GetXaxis()->CenterTitle(true);
        cMap->cd(i + 3);
        fh2_EnergyVsStrip_all[i]->Draw("colz");
    }
    for (Int_t i = 0; i < fNbDets; i++)
    {
        sprintf(Name1, "strip_count_det_%d", i + 1);
        sprintf(Name2, "Number of counts per strip for TTT10 Det: %d, all events", i + 1);
        fh1_Strip[i] = R3B::root_owned<TH1F>(Name1, Name2, fNbStrips, 1, fNbStrips + 1);
        fh1_Strip[i]->GetXaxis()->SetTitle("Strip number");
        fh1_Strip[i]->GetYaxis()->SetTitle("Counts");
        fh1_Strip[i]->GetYaxis()->SetTitleOffset(1.4);
        fh1_Strip[i]->GetXaxis()->CenterTitle(true);
        cMap->cd(i + 5);
        fh1_Strip[i]->Draw();
    }
    for (Int_t i = 0; i < fNbDets; i++)
    {
        sprintf(Name1, "time_vs_strip_det_%d", i + 1);
        sprintf(Name2, "Time - trigger time vs strip number for TTT10 Det: %d", i + 1);
        fh2_TimeVsStrip[i] = R3B::root_owned<TH2F>(Name1, Name2, fNbStrips, 1, fNbStrips + 1, binsT, minT, maxT);
        fh2_TimeVsStrip[i]->GetXaxis()->SetTitle("Strip number");
        fh2_TimeVsStrip[i]->GetYaxis()->SetTitle("Time - trigger time [channels]");
        fh2_TimeVsStrip[i]->GetYaxis()->SetTitleOffset(1.4);
        fh2_TimeVsStrip[i]->GetXaxis()->CenterTitle(true);
        cMap->cd(i + 7);
        fh2_TimeVsStrip[i]->Draw("colz");
    }

    sprintf(Name1, "highest_e_correlation");
    sprintf(Name2, "Highest energy of Det 2 Vs highest energy of Det 1");
    fh2_E2VsE1 = R3B::root_owned<TH2F>(Name1, Name2, binsE, minE, maxE, binsE, minE, maxE);
    fh2_E2VsE1->GetXaxis()->SetTitle("Highest energy hit of Det 1 [channels]");
    fh2_E2VsE1->GetYaxis()->SetTitle("Highest energy hit of Det 2 [channels]");
    fh2_E2VsE1->GetYaxis()->SetTitleOffset(1.4);
    fh2_E2VsE1->GetXaxis()->CenterTitle(true);
    cMap->cd(9);
    fh2_E2VsE1->Draw("colz");

    sprintf(Name1, "highest_e_strip_correlation");
    sprintf(Name2, "Strip with highest energy of Det 2 Vs strip with highest energy of Det 1");
    fh2_E2StripVsE1Strip =
        R3B::root_owned<TH2F>(Name1, Name2, fNbStrips, 1, fNbStrips + 1, fNbStrips, 1, fNbStrips + 1);
    fh2_E2StripVsE1Strip->GetXaxis()->SetTitle("strip number of Det 1");
    fh2_E2StripVsE1Strip->GetYaxis()->SetTitle("strip number of Det 2");
    fh2_E2StripVsE1Strip->GetYaxis()->SetTitleOffset(1.4);
    fh2_E2StripVsE1Strip->GetXaxis()->CenterTitle(true);
    cMap->cd(10);
    fh2_E2StripVsE1Strip->Draw("colz");

    for (Int_t i = 0; i < fNbDets; i++)
    {
        sprintf(Name1, "multiplicity_det_%d", i + 1);
        sprintf(Name2, "Multiplicity vs strip number for TTT10 Det: %d", i + 1);
        fh2_MultVsStrip[i] = R3B::root_owned<TH2F>(Name1, Name2, fNbStrips, 1, fNbStrips + 1, binsM, minM, maxM);
        fh2_MultVsStrip[i]->GetXaxis()->SetTitle("Strip number");
        fh2_MultVsStrip[i]->GetYaxis()->SetTitle("Multiplicity");
        fh2_MultVsStrip[i]->GetYaxis()->SetTitleOffset(1.4);
        fh2_MultVsStrip[i]->GetXaxis()->CenterTitle(true);
        cMap->cd(i + 11);
        fh2_MultVsStrip[i]->Draw("colz");
    }

    //============================================= CAL DATA ====================================================

    fh2_CalEnergyVsStrip.resize(fNbDets);
    fh2_CalTimeVsStrip.resize(fNbDets);

    for (Int_t i = 0; i < fNbDets; i++)
    {
        sprintf(Name1, "fh2_cal_energy_vs_strip_det_%d", i + 1);
        sprintf(Name2, "Calibrated energy vs strip number for TTT10 Det: %d", i + 1);
        fh2_CalEnergyVsStrip[i] =
            R3B::root_owned<TH2F>(Name1, Name2, fNbStrips, 1, fNbStrips + 1, binsCalE, minCalE, maxCalE);
        fh2_CalEnergyVsStrip[i]->GetXaxis()->SetTitle("Strip number");
        fh2_CalEnergyVsStrip[i]->GetYaxis()->SetTitle("Energy [keV]");
        fh2_CalEnergyVsStrip[i]->GetYaxis()->SetTitleOffset(1.4);
        fh2_CalEnergyVsStrip[i]->GetXaxis()->CenterTitle(true);
        cCal->cd(i + 1);
        fh2_CalEnergyVsStrip[i]->Draw("colz");
    }
    for (Int_t i = 0; i < fNbDets; i++)
    {
        sprintf(Name1, "fh2_cal_time_vs_strip_det_%d", i + 1);
        sprintf(Name2, "Time vs strip number for TTT10 Det: %d", i + 1);
        fh2_CalTimeVsStrip[i] =
            R3B::root_owned<TH2F>(Name1, Name2, fNbStrips, 1, fNbStrips + 1, binsCalT, minCalT, maxCalT);
        fh2_CalTimeVsStrip[i]->GetXaxis()->SetTitle("Strip number");
        fh2_CalTimeVsStrip[i]->GetYaxis()->SetTitle("Time [ns]");
        fh2_CalTimeVsStrip[i]->GetYaxis()->SetTitleOffset(1.4);
        fh2_CalTimeVsStrip[i]->GetXaxis()->CenterTitle(true);
        cCal->cd(i + 3);
        fh2_CalTimeVsStrip[i]->Draw("colz");
    }

    sprintf(Name1, "fh2_cal_mult_correlation");
    sprintf(Name2, "Multiplicity of Det 2 Vs Multiplicity of Det 1");
    fh2_CalMultCorr = R3B::root_owned<TH2F>(Name1, Name2, binsM, minM, maxM, binsM, minM, maxM);
    fh2_CalMultCorr->GetXaxis()->SetTitle("Multiplicity of Det 1");
    fh2_CalMultCorr->GetYaxis()->SetTitle("Multiplicity of Det 2");
    fh2_CalMultCorr->GetYaxis()->SetTitleOffset(1.4);
    fh2_CalMultCorr->GetXaxis()->CenterTitle(true);
    cCal->cd(5);
    fh2_CalMultCorr->Draw("colz");

    sprintf(Name1, "fh2_cal_energy_correlation");
    sprintf(Name2, "Energy of Det 2 Vs Energy of Det 1");
    fh2_CalECorr = R3B::root_owned<TH2F>(Name1, Name2, binsCalE, minCalE, maxCalE, binsCalE, minCalE, maxCalE);
    fh2_CalECorr->GetXaxis()->SetTitle("Energy of Det 1 [keV]");
    fh2_CalECorr->GetYaxis()->SetTitle("Energy of Det 2 [keV]");
    fh2_CalECorr->GetYaxis()->SetTitleOffset(1.4);
    fh2_CalECorr->GetXaxis()->CenterTitle(true);
    cCal->cd(6);
    fh2_CalECorr->Draw("colz");

    sprintf(Name1, "fh2_cal_time_correlation");
    sprintf(Name2, "Time of Det 2 Vs Time of Det 1");
    fh2_CalTCorr = R3B::root_owned<TH2F>(Name1, Name2, binsCalT, minCalT, maxCalT, binsCalT, minCalT, maxCalT);
    fh2_CalTCorr->GetXaxis()->SetTitle("Time of Det 1 [ns]");
    fh2_CalTCorr->GetYaxis()->SetTitle("Time of Det 2 [ns]");
    fh2_CalTCorr->GetYaxis()->SetTitleOffset(1.4);
    fh2_CalTCorr->GetXaxis()->CenterTitle(true);
    cCal->cd(7);
    fh2_CalTCorr->Draw("colz");

    //============================================= HIT DATA ====================================================

    fh2_HitEnergyVsPos.resize(fNbDets);
    fh2_HitTimeVsPos.resize(fNbDets);

    for (Int_t i = 0; i < fNbDets; i++)
    {
        sprintf(Name1, "fh2_hit_energy_vs_pos_det_%d", i + 1);
        sprintf(Name2, "Hit energy vs Position for TTT10 Det: %d", i + 1);
        fh2_HitEnergyVsPos[i] = R3B::root_owned<TH2F>(Name1, Name2, binsX, minX, maxX, binsCalE, minCalE, maxCalE);
        fh2_HitEnergyVsPos[i]->GetXaxis()->SetTitle("Position [mm]");
        fh2_HitEnergyVsPos[i]->GetYaxis()->SetTitle("Energy [keV]");
        fh2_HitEnergyVsPos[i]->GetYaxis()->SetTitleOffset(1.4);
        fh2_HitEnergyVsPos[i]->GetXaxis()->CenterTitle(true);
        cHit->cd(i + 1);
        fh2_HitEnergyVsPos[i]->Draw("colz");
    }
    for (Int_t i = 0; i < fNbDets; i++)
    {
        sprintf(Name1, "fh2_hit_time_vs_pos_det_%d", i + 1);
        sprintf(Name2, "Time vs Position for TTT10 Det: %d", i + 1);
        fh2_HitTimeVsPos[i] = R3B::root_owned<TH2F>(Name1, Name2, binsX, minX, maxX, binsCalT, minCalT, maxCalT);
        fh2_HitTimeVsPos[i]->GetXaxis()->SetTitle("Position [mm]");
        fh2_HitTimeVsPos[i]->GetYaxis()->SetTitle("Time [ns]");
        fh2_HitTimeVsPos[i]->GetYaxis()->SetTitleOffset(1.4);
        fh2_HitTimeVsPos[i]->GetXaxis()->CenterTitle(true);
        cHit->cd(i + 3);
        fh2_HitTimeVsPos[i]->Draw("colz");
    }

    sprintf(Name1, "fh2_hit_mult_correlation");
    sprintf(Name2, "Multiplicity of Det 2 Vs Multiplicity of Det 1");
    fh2_HitMultCorr = R3B::root_owned<TH2F>(Name1, Name2, binsM, minM, maxM, binsM, minM, maxM);
    fh2_HitMultCorr->GetXaxis()->SetTitle("Multiplicity of Det 1");
    fh2_HitMultCorr->GetYaxis()->SetTitle("Multiplicity of Det 2");
    fh2_HitMultCorr->GetYaxis()->SetTitleOffset(1.4);
    fh2_HitMultCorr->GetXaxis()->CenterTitle(true);
    cHit->cd(5);
    fh2_HitMultCorr->Draw("colz");

    sprintf(Name1, "fh2_hit_energy_correlation");
    sprintf(Name2, "Energy of Det 2 Vs Energy of Det 1");
    fh2_HitECorr = R3B::root_owned<TH2F>(Name1, Name2, binsCalE, minCalE, maxCalE, binsCalE, minCalE, maxCalE);
    fh2_HitECorr->GetXaxis()->SetTitle("Energy of Det 1 [keV]");
    fh2_HitECorr->GetYaxis()->SetTitle("Energy of Det 2 [keV]");
    fh2_HitECorr->GetYaxis()->SetTitleOffset(1.4);
    fh2_HitECorr->GetXaxis()->CenterTitle(true);
    cHit->cd(6);
    fh2_HitECorr->Draw("colz");

    sprintf(Name1, "fh2_hit_z_correlation");
    sprintf(Name2, "Z of Det 2 Vs Z of Det 1");
    fh2_HitZCorr = R3B::root_owned<TH2F>(Name1, Name2, binsZ, minZ, maxZ, binsZ, minZ, maxZ);
    fh2_HitZCorr->GetXaxis()->SetTitle("Z of Det 1");
    fh2_HitZCorr->GetYaxis()->SetTitle("Z of Det 2");
    fh2_HitZCorr->GetYaxis()->SetTitleOffset(1.4);
    fh2_HitZCorr->GetXaxis()->CenterTitle(true);
    cHit->cd(7);
    fh2_HitZCorr->Draw("colz");

    sprintf(Name1, "fh2_hit_pos_correlation");
    sprintf(Name2, "Position of Det 2 Hit Vs Position of Det 1 Hit");
    fh2_HitPosCorr = R3B::root_owned<TH2F>(Name1, Name2, binsX, minX, maxX, binsX, minX, maxX);
    fh2_HitPosCorr->GetXaxis()->SetTitle("x Pos of Det 1 [mm]");
    fh2_HitPosCorr->GetYaxis()->SetTitle("x Pos of Det 2 [mm]");
    fh2_HitPosCorr->GetYaxis()->SetTitleOffset(1.4);
    fh2_HitPosCorr->GetXaxis()->CenterTitle(true);
    cHit->cd(8);
    fh2_HitPosCorr->Draw("colz");

    auto cTttx_zcharge_tpat = new TCanvas("Tttx_Charge_vs_tpat", "Charges Z vs tpat", 10, 10, 500, 500);
    sprintf(Name1, "fh2_tttx_Charge_vs_tpat");
    sprintf(Name2, "Charges Z vs tpat");
    fh2_Zcharge_tpat = R3B::root_owned<TH2F>(Name1, Name2, 17, -0.5, 16.5, 1000, 0, 12);
    fh2_Zcharge_tpat->GetXaxis()->SetTitle("Tpat");
    fh2_Zcharge_tpat->GetYaxis()->SetTitle("Charge-Z");
    fh2_Zcharge_tpat->GetXaxis()->SetTitleOffset(1.1);
    fh2_Zcharge_tpat->GetYaxis()->SetTitleOffset(1.15);
    fh2_Zcharge_tpat->GetXaxis()->CenterTitle(true);
    fh2_Zcharge_tpat->Draw("colz");
    mainfolTTTX->Add(cTttx_zcharge_tpat);

    //===============================================================================================================
    run->AddObject(mainfolTTTX);

    // Register command to reset histograms
    run->GetHttpServer()->RegisterCommand("Reset_Tttx_HIST", Form("/Objects/%s/->Reset_Histo()", GetName()));
    R3BLOG(info, "Init Complete");

    return kSUCCESS;
}

void R3BTttxOnlineSpectra::Reset_Histo()
{
    R3BLOG(info, "Reset_Histo");
    // Map data
    for (Int_t i = 0; i < fNbDets; i++)
    {
        fh2_EnergyVsStrip_recent[i]->Reset();
        fh2_EnergyVsStrip_all[i]->Reset();
        fh1_Strip[i]->Reset();
        fh2_TimeVsStrip[i]->Reset();
        fh2_MultVsStrip[i]->Reset();
    }
    fh2_E2VsE1->Reset();
    fh2_E2StripVsE1Strip->Reset();

    // Cal data
    for (Int_t i = 0; i < fNbDets; i++)
    {
        fh2_CalEnergyVsStrip[i]->Reset();
        fh2_CalTimeVsStrip[i]->Reset();
    }
    fh2_CalMultCorr->Reset();
    fh2_CalECorr->Reset();
    fh2_CalTCorr->Reset();

    // Hit data
    for (Int_t i = 0; i < fNbDets; i++)
    {
        fh2_HitEnergyVsPos[i]->Reset();
        fh2_HitTimeVsPos[i]->Reset();
    }
    fh2_HitMultCorr->Reset();
    fh2_HitECorr->Reset();
    fh2_HitZCorr->Reset();
    fh2_HitPosCorr->Reset();
    fh2_Zcharge_tpat->Reset();
}

void R3BTttxOnlineSpectra::Exec(Option_t* option)
{

    std::vector<int> tpatindex;
    int tpatbin = 0;
    if (header)
    {
        for (Int_t i = 0; i < 16; i++)
        {
            tpatbin = (header->GetTpat() & (1 << i));
            if (tpatbin != 0)
                tpatindex.push_back(i + 1);
        }
    }

    // Fill Mapped Data
    if (!fMappedItemsTttx || fMappedItemsTttx->GetEntriesFast() == 0)
    {
        return;
    }
    UInt_t det = 0;
    UInt_t strip = 0;

    energy.clear();
    time.clear();
    mult.clear();
    highest_e.clear();
    highest_e_strip.clear();

    energy.resize(fNbDets);
    time.resize(fNbDets);
    mult.resize(fNbDets);
    for (int idet = 0; idet < fNbDets; idet++)
    {
        energy.at(idet).resize(fNbStrips);
        time.at(idet).resize(fNbStrips + fNbTref + fNbTrig);
        mult.at(idet).resize(fNbStrips);
    }
    highest_e.resize(fNbDets);
    highest_e_strip.resize(fNbDets);

    for (Int_t j = 0; j < fNbDets; j++)
    {
        highest_e[j] = 0;
        highest_e_strip[j] = 0;
        for (Int_t i = 0; i < fNbStrips; i++)
        {
            energy[j][i] = NAN; // mult=1 !!!
            mult[j][i] = 0;
        }
        for (Int_t i = 0; i < fNbStrips + fNbTref + fNbTrig; i++)
        {
            time[j][i] = NAN; // mult=1 !!!
        }
    }

    Int_t nHits = fMappedItemsTttx->GetEntriesFast();
    for (Int_t ihit = 0; ihit < nHits; ihit++)
    {
        R3BTttxMappedData* hit = dynamic_cast<R3BTttxMappedData*>(fMappedItemsTttx->At(ihit));
        if (!hit)
        {
            continue;
        }
        det = hit->GetDetID() - 1;
        strip = hit->GetStripID() - 1; // stripID - 1 --> 32 => Tref, 33 => Ttrig

        if (strip < fNbStrips)
        {
            fh2_EnergyVsStrip_recent[det]->Fill(strip + 1, hit->GetEnergy());
            fh2_EnergyVsStrip_all[det]->Fill(strip + 1, hit->GetEnergy());
            fh1_Strip[det]->Fill(strip + 1);
            energy[det][strip] = hit->GetEnergy() & 0xffff; // Take only lower 16 bits
            mult[det][strip]++;
            if (!isnan(energy[det][strip]) && energy[det][strip] > highest_e[det])
            {
                highest_e[det] = hit->GetEnergy();
                highest_e_strip[det] = hit->GetStripID();
            }
        }
        time[det][strip] = hit->GetTime();
    }
    for (Int_t i_det = 0; i_det < fNbDets; i_det++)
    {
        for (Int_t i_strip = 0; i_strip < fNbStrips; i_strip++)
        {
            if (mult[i_det][i_strip] > 0)
                fh2_MultVsStrip[i_det]->Fill(i_strip + 1, mult[i_det][i_strip]);
            if (!isnan(time[i_det][i_strip]) && !isnan(time[i_det][fNbStrips + fNbTref + fNbTrig - 1]))
                fh2_TimeVsStrip[i_det]->Fill(
                    i_strip + 1,
                    time[i_det][i_strip] - time[i_det][fNbStrips + fNbTref + fNbTrig - 1]); // mult=1 !
        }
    }
    if (!isnan(highest_e[0]) && !isnan(highest_e[1]))
    {
        fh2_E2VsE1->Fill(highest_e[0], highest_e[1]);
        if (abs(highest_e[0] - highest_e[1]) < 1000.)
        {
            auto mean_charge = TMath::Sqrt(highest_e[0] * highest_e[1]) / 1000.;
            for (const auto& itpat : tpatindex)
                fh2_Zcharge_tpat->Fill(itpat, mean_charge);
        }
    }
    fh2_E2StripVsE1Strip->Fill(highest_e_strip[0], highest_e_strip[1]);

    for (Int_t i_det = 0; i_det < fNbDets; i_det++)
    {
        if (fh2_EnergyVsStrip_recent[i_det]->GetEntries() > fNbEntries)
        {
            fh2_EnergyVsStrip_recent[i_det]->Reset();
        }
    }

    // Fill Cal Data
    if (!fCalItemsTttx || fCalItemsTttx->GetEntriesFast() == 0)
    {
        return;
    }
    cal_mult.clear();
    cal_mult.resize(fNbDets);

    for (Int_t j = 0; j < fNbDets; j++)
        cal_mult[j] = 0;

    Int_t nCalHits = fCalItemsTttx->GetEntriesFast();
    for (Int_t ihit = 0; ihit < nCalHits; ihit++)
    {
        R3BTttxCalData* hit = dynamic_cast<R3BTttxCalData*>(fCalItemsTttx->At(ihit));
        if (!hit)
        {
            continue;
        }
        det = hit->GetDetID() - 1;
        strip = hit->GetStripID() - 1; // stripID - 1 --> 32 => Tref, 33 => Ttrig
        if (strip < fNbStrips)
        {
            fh2_CalEnergyVsStrip[det]->Fill(strip + 1, hit->GetEnergy());
            fh2_CalTimeVsStrip[det]->Fill(strip + 1, hit->GetTime());
            cal_mult[det]++;
        }
        for (int jhit = ihit + 1; jhit < nCalHits; jhit++)
        {
            R3BTttxCalData* hit2 = dynamic_cast<R3BTttxCalData*>(fCalItemsTttx->At(jhit));
            if (det == hit2->GetDetID() - 1)
            {
                continue;
            }
            if (det == 0)
            {
                fh2_CalECorr->Fill(hit->GetEnergy(), hit2->GetEnergy());
                fh2_CalTCorr->Fill(hit->GetTime(), hit2->GetTime());
            }
            else
            {
                fh2_CalECorr->Fill(hit2->GetEnergy(), hit->GetEnergy());
                fh2_CalTCorr->Fill(hit2->GetTime(), hit->GetTime());
            }
        }
    }
    fh2_CalMultCorr->Fill(cal_mult[0], cal_mult[1]);

    // Fill Hit Data
    if (!fHitItemsTttx || fHitItemsTttx->GetEntriesFast() == 0)
    {
        return;
    }
    hit_mult.clear();
    hit_mult.resize(fNbDets);
    Double_t pos = 0;
    for (Int_t j = 0; j < fNbDets; j++)
        hit_mult[j] = 0;

    Int_t nHitHits = fHitItemsTttx->GetEntriesFast();
    for (Int_t ihit = 0; ihit < nHitHits; ihit++)
    {
        R3BTttxHitData* hit = dynamic_cast<R3BTttxHitData*>(fHitItemsTttx->At(ihit));
        if (!hit)
        {
            continue;
        }
        if (hit->GetDetID() == 1 || hit->GetDetID() == 2)
        { // don't consider DetID==0 (both)
            det = hit->GetDetID() - 1;
            pos = hit->GetX();

            fh2_HitEnergyVsPos[det]->Fill(pos, hit->GetEnergy());
            fh2_HitTimeVsPos[det]->Fill(pos, hit->GetTime());
            hit_mult[det]++;
            for (int jhit = ihit + 1; jhit < nHitHits; jhit++)
            {
                R3BTttxHitData* hit2 = dynamic_cast<R3BTttxHitData*>(fHitItemsTttx->At(jhit));
                if ((hit->GetDetID() == 1 && hit2->GetDetID() == 1) || (hit->GetDetID() == 2 && hit2->GetDetID() == 2))
                {
                    continue;
                }
                if (hit->GetDetID() == 1 && hit2->GetDetID() == 2)
                {
                    fh2_HitECorr->Fill(hit->GetEnergy(), hit2->GetEnergy());
                    fh2_HitZCorr->Fill(hit->GetChargeZ(), hit2->GetChargeZ());
                    fh2_HitPosCorr->Fill(hit->GetX(), hit2->GetX());
                }
                if (hit->GetDetID() == 1 && hit2->GetDetID() == 2)
                {
                    fh2_HitECorr->Fill(hit2->GetEnergy(), hit->GetEnergy());
                    fh2_HitZCorr->Fill(hit2->GetChargeZ(), hit->GetChargeZ());
                    fh2_HitPosCorr->Fill(hit2->GetX(), hit->GetX());
                }
            }
        }
    }
    fh2_HitMultCorr->Fill(hit_mult[0], hit_mult[1]);
}

void R3BTttxOnlineSpectra::FinishEvent()
{
    if (fMappedItemsTttx)
    {
        fMappedItemsTttx->Clear();
    }
    if (fCalItemsTttx)
    {
        fCalItemsTttx->Clear();
    }
    if (fHitItemsTttx)
    {
        fHitItemsTttx->Clear();
    }
}

void R3BTttxOnlineSpectra::FinishTask()
{
    if (fMappedItemsTttx)
    {
        cMap->Write();
    }
    if (fCalItemsTttx)
    {
        cCal->Write();
    }
    if (fHitItemsTttx)
    {
        cHit->Write();
        fh2_Zcharge_tpat->Write();
    }
}

ClassImp(R3BTttxOnlineSpectra);
