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
// -----              R3BGeneralOnlineSpectra             -----
// -----    Created 05/03/22 by J.L. Rodriguez-Sanchez    -----
// -----           Fill R3B online histograms             -----
// ------------------------------------------------------------

#include "R3BGeneralOnlineSpectra.h"
#include "R3BAmsOnlineSpectra.h"
#include "R3BCalifaOnlineSpectra.h"
#include "R3BEventHeader.h"
#include "R3BFootOnlineSpectra.h"
#include "R3BIncomingTrackingOnlineSpectra.h"
#include "R3BLogger.h"
#include "R3BLosOnlineSpectra.h"
#include "R3BMusicOnlineSpectra.h"
#include "R3BMwpcCorrelationOnlineSpectra.h"
#include "R3BMwpcOnlineSpectra.h"
#include "R3BTofDOnlineSpectra.h"
#include "R3BTwimOnlineSpectra.h"
#include "R3BTwimvsFootOnlineSpectra.h"
#include "R3BWRData.h"

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
#include "TLegend.h"
#include "TLegendEntry.h"
#include "TMath.h"
#include "TRandom.h"
#include "TVector3.h"

R3BGeneralOnlineSpectra::R3BGeneralOnlineSpectra()
    : R3BGeneralOnlineSpectra("GeneralOnlineSpectra", 1)
{
}

R3BGeneralOnlineSpectra::R3BGeneralOnlineSpectra(const TString& name, Int_t iVerbose)
    : FairTask(name, iVerbose)
    , fEventHeader(nullptr)
    , fMwpc0Online(NULL)
    , fMwpc01Online(NULL)
    , fMwpc02Online(NULL)
    , fMwpc12Online(NULL)
    , fMwpc1Online(NULL)
    , fMwpc2Online(NULL)
    , fTwimOnline(NULL)
    , fMusicOnline(NULL)
    , fAmsOnline(NULL)
    , fLosOnline(NULL)
    , fFootOnline(NULL)
    , fCalifaOnline(NULL)
    , fTofdOnlineSpectra(NULL)
    , fIncomingTrackingOnline(NULL)
    , fTwimvsFootOnlineSpectra(NULL)
    , fWRItemsMaster(NULL)
    , fWRItemsSofia(NULL)
    , fWRItemsCalifa(NULL)
    , fWRItemsNeuland(NULL)
    , fWRItemsS2(NULL)
    , fWRItemsS8(NULL)
    , fNEvents(0)
{
}

R3BGeneralOnlineSpectra::~R3BGeneralOnlineSpectra()
{
    R3BLOG(debug1, "");
    if (fWRItemsMaster)
        delete fWRItemsMaster;
    if (fWRItemsSofia)
        delete fWRItemsSofia;
    if (fWRItemsCalifa)
        delete fWRItemsCalifa;
    if (fWRItemsNeuland)
        delete fWRItemsNeuland;
    if (fWRItemsS2)
        delete fWRItemsS2;
    if (fWRItemsS8)
        delete fWRItemsS8;
}

InitStatus R3BGeneralOnlineSpectra::Init()
{
    R3BLOG(info, "");
    FairRootManager* mgr = FairRootManager::Instance();
    R3BLOG_IF(fatal, NULL == mgr, "FairRootManager not found");

    fEventHeader = dynamic_cast<R3BEventHeader*>(mgr->GetObject("EventHeader."));
    R3BLOG_IF(warn, NULL == fEventHeader, "EventHeader. not found");

    FairRunOnline* run = FairRunOnline::Instance();
    R3BLOG_IF(fatal, NULL == run, "FairRunOnline not found");
    run->GetHttpServer()->Register("", this);

    // get access to WR-Master data
    fWRItemsMaster = dynamic_cast<TClonesArray*>(mgr->GetObject("WRMasterData"));
    R3BLOG_IF(warn, !fWRItemsMaster, "WRMasterData not found");

    // get access to WR-Sofia data
    fWRItemsSofia = dynamic_cast<TClonesArray*>(mgr->GetObject("SofWRData"));
    R3BLOG_IF(warn, !fWRItemsSofia, "SofWRData not found");

    // get access to WR-Califa data
    fWRItemsCalifa = dynamic_cast<TClonesArray*>(mgr->GetObject("WRCalifaData"));
    R3BLOG_IF(warn, !fWRItemsCalifa, "WRCalifaData not found");

    // get access to WR-Neuland data
    fWRItemsNeuland = dynamic_cast<TClonesArray*>(mgr->GetObject("WRNeulandData"));
    R3BLOG_IF(warn, !fWRItemsNeuland, "WRNeulandData not found");

    // get access to WR-S2 data
    fWRItemsS2 = dynamic_cast<TClonesArray*>(mgr->GetObject("WRS2Data"));
    R3BLOG_IF(warn, !fWRItemsS2, "WRS2Data not found");

    // get access to WR-S8 data
    fWRItemsS8 = dynamic_cast<TClonesArray*>(mgr->GetObject("WRS8Data"));
    R3BLOG_IF(warn, !fWRItemsS8, "WRS8Data not found");

    // Looking for Mwpc0 online
    fMwpc0Online = dynamic_cast<R3BMwpcOnlineSpectra*>(FairRunOnline::Instance()->GetTask("Mwpc0OnlineSpectra"));
    R3BLOG_IF(warn, !fMwpc0Online, "Mwpc0OnlineSpectra not found");

    // Looking for Mwpc0_1 online
    fMwpc01Online = dynamic_cast<R3BMwpcCorrelationOnlineSpectra*>(
        FairRunOnline::Instance()->GetTask("Mwpc0_1CorrelationOnlineSpectra"));
    R3BLOG_IF(warn, !fMwpc01Online, "Mwpc0_1CorrelationOnlineSpectra not found");

    // Looking for Mwpc0_2 online
    fMwpc02Online = dynamic_cast<R3BMwpcCorrelationOnlineSpectra*>(
        FairRunOnline::Instance()->GetTask("Mwpc0_2CorrelationOnlineSpectra"));
    R3BLOG_IF(warn, !fMwpc02Online, "Mwpc0_2CorrelationOnlineSpectra not found");

    // Looking for Mwpc1_2 online
    fMwpc12Online = dynamic_cast<R3BMwpcCorrelationOnlineSpectra*>(
        FairRunOnline::Instance()->GetTask("Mwpc1_2CorrelationOnlineSpectra"));
    R3BLOG_IF(warn, !fMwpc12Online, "Mwpc1_2CorrelationOnlineSpectra not found");

    // Looking for Mwpc1 online
    fMwpc1Online = dynamic_cast<R3BMwpcOnlineSpectra*>(FairRunOnline::Instance()->GetTask("Mwpc1OnlineSpectra"));
    R3BLOG_IF(warn, !fMwpc1Online, "Mwpc1OnlineSpectra not found");

    // Looking for Mwpc2 online
    fMwpc2Online = dynamic_cast<R3BMwpcOnlineSpectra*>(FairRunOnline::Instance()->GetTask("Mwpc2OnlineSpectra"));
    if (!fMwpc2Online)
        LOG(warn) << "R3BGeneralOnlineSpectra::Mwpc2OnlineSpectra not found";

    // Looking for Twim online
    fTwimOnline = dynamic_cast<R3BTwimOnlineSpectra*>(FairRunOnline::Instance()->GetTask("TwimOnlineSpectra"));
    R3BLOG_IF(warn, !fTwimOnline, "TwimOnlineSpectra not found");

    // Looking for Music online
    fMusicOnline = dynamic_cast<R3BMusicOnlineSpectra*>(FairRunOnline::Instance()->GetTask("MusicOnlineSpectra"));
    R3BLOG_IF(warn, !fMusicOnline, "MusicOnlineSpectra not found");

    // Looking for AMS online
    fAmsOnline = dynamic_cast<R3BAmsOnlineSpectra*>(FairRunOnline::Instance()->GetTask("AmsOnlineSpectra"));
    R3BLOG_IF(warn, !fAmsOnline, "AmsOnlineSpectra not found");

    // Looking for LOS online
    fLosOnline = dynamic_cast<R3BLosOnlineSpectra*>(FairRunOnline::Instance()->GetTask("LosOnlineSpectra"));
    R3BLOG_IF(warn, !fLosOnline, "LosOnlineSpectra not found");

    // Looking for FOOT online
    fFootOnline = dynamic_cast<R3BFootOnlineSpectra*>(FairRunOnline::Instance()->GetTask("FootOnlineSpectra"));
    R3BLOG_IF(warn, !fFootOnline, "FootOnlineSpectra not found");

    // Looking for CALIFA online
    fCalifaOnline = dynamic_cast<R3BCalifaOnlineSpectra*>(FairRunOnline::Instance()->GetTask("CALIFAOnlineSpectra"));
    R3BLOG_IF(warn, !fCalifaOnline, "CALIFAOnlineSpectra not found");

    // Looking for TOFD online
    fTofdOnlineSpectra = dynamic_cast<R3BTofDOnlineSpectra*>(FairRunOnline::Instance()->GetTask("TofdOnlineSpectra"));
    R3BLOG_IF(warn, !fTofdOnlineSpectra, "TofdOnlineSpectra not found");

    // Looking for Incoming Tracking online
    fIncomingTrackingOnline = dynamic_cast<R3BIncomingTrackingOnlineSpectra*>(
        FairRunOnline::Instance()->GetTask("IncomingTrackingOnlineSpectra"));
    R3BLOG_IF(warn, !fIncomingTrackingOnline, "IncomingTrackingOnlineSpectra not found");

    // Looking for Incoming Tracking online
    fTwimvsFootOnlineSpectra =
        dynamic_cast<R3BTwimvsFootOnlineSpectra*>(FairRunOnline::Instance()->GetTask("TwimvsFootOnlineSpectra"));
    R3BLOG_IF(warn, !fTwimvsFootOnlineSpectra, "TwimvsFootOnlineSpectra not found");

    // Create histograms for detectors
    char Name1[255];
    char Name2[255];
    char Name3[255];

    // Tpats
    cTpat = new TCanvas("TPats", "Tpat information", 10, 10, 800, 700);
    cTpat->Divide(2, 1);
    cTpat->cd(1);
    fh1_tpat = new TH1F("fh1_tpat", "TPat information", 17, -0.5, 16.5);
    fh1_tpat->GetXaxis()->SetTitle("Trigger number (TPat)");
    fh1_tpat->GetYaxis()->SetTitle("Counts");
    fh1_tpat->GetXaxis()->CenterTitle(true);
    fh1_tpat->GetYaxis()->CenterTitle(true);
    fh1_tpat->GetXaxis()->SetLabelSize(0.04);
    fh1_tpat->GetXaxis()->SetTitleSize(0.04);
    fh1_tpat->GetYaxis()->SetTitleOffset(1.1);
    fh1_tpat->GetXaxis()->SetTitleOffset(1.1);
    fh1_tpat->GetYaxis()->SetLabelSize(0.04);
    fh1_tpat->GetYaxis()->SetTitleSize(0.04);
    fh1_tpat->SetFillColor(kBlue + 2);
    fh1_tpat->Draw("");
    cTpat->cd(2)->SetLogz();
    fh2_tpat = new TH2F(
        "fh2_tpat",
        "TPat correlation (diagonal = no correlating tpats); Lower trigger number (TPat); Higher trigger number (TPat)",
        17,
        -0.5,
        16.5,
        17,
        -0.5,
        16.5);
    fh2_tpat->GetXaxis()->CenterTitle(true);
    fh2_tpat->GetYaxis()->CenterTitle(true);
    fh2_tpat->GetXaxis()->SetLabelSize(0.04);
    fh2_tpat->GetXaxis()->SetTitleSize(0.04);
    fh2_tpat->GetYaxis()->SetTitleOffset(1.1);
    fh2_tpat->GetXaxis()->SetTitleOffset(1.1);
    fh2_tpat->GetYaxis()->SetLabelSize(0.04);
    fh2_tpat->GetYaxis()->SetTitleSize(0.04);
    fh2_tpat->Draw("colz");

    // Triggers
    cTrigger = new TCanvas("Triggers", "Trigger information", 10, 10, 800, 700);
    fh1_trigger = new TH1F("fh1_trigger", "Trigger information", 17, -0.5, 16.5);
    fh1_trigger->GetXaxis()->SetTitle("Trigger number");
    fh1_trigger->GetYaxis()->SetTitle("Counts");
    fh1_trigger->GetXaxis()->CenterTitle(true);
    fh1_trigger->GetYaxis()->CenterTitle(true);
    fh1_trigger->GetXaxis()->SetLabelSize(0.04);
    fh1_trigger->GetXaxis()->SetTitleSize(0.04);
    fh1_trigger->GetYaxis()->SetTitleOffset(1.1);
    fh1_trigger->GetXaxis()->SetTitleOffset(1.1);
    fh1_trigger->GetYaxis()->SetLabelSize(0.04);
    fh1_trigger->GetYaxis()->SetTitleSize(0.04);
    fh1_trigger->SetFillColor(kBlue + 2);
    fh1_trigger->Draw("");

    // Difference between master and sofia WRs
    cWr = new TCanvas("WR_Master_Sofia", "WR_Master_Sofia", 10, 10, 500, 500);
    fh1_wr[0] = new TH1F(
        "fh1_WRSE_Master_Sofia", "WR-Master - WRSE-Sofia(back) , WR-Master - WRME-Sofia(red) ", 1200, -4100, 4100);
    fh1_wr[0]->GetXaxis()->SetTitle("WRs difference");
    fh1_wr[0]->GetYaxis()->SetTitle("Counts");
    fh1_wr[0]->GetYaxis()->SetTitleOffset(1.3);
    fh1_wr[0]->GetXaxis()->CenterTitle(true);
    fh1_wr[0]->GetYaxis()->CenterTitle(true);
    fh1_wr[0]->SetFillColor(29);
    fh1_wr[0]->SetLineColor(1);
    fh1_wr[0]->SetLineWidth(2);
    fh1_wr[0]->Draw("");
    fh1_wr[1] = new TH1F("fh1_WRME_Master_Sofia", "WR-Master - WRME-Sofia", 1200, -4100, 4100);
    fh1_wr[1]->SetLineColor(2);
    fh1_wr[1]->SetLineWidth(2);
    fh1_wr[1]->Draw("same");

    // Difference between Califa-Sofia WRs
    sprintf(Name1, "WRs_Sofia_vs_others");
    cWrs = new TCanvas(Name1, Name1, 10, 10, 500, 500);
    sprintf(Name2, "fh1_WR_Sofia_Wixhausen");
    sprintf(Name3, "WR-Sofia - WR-Other"); // Messel (blue), Wixhausen (red)
    fh1_wrs[0] = new TH1F(Name2, Name3, 1200, -4100, 4100);
    fh1_wrs[0]->GetXaxis()->SetTitle("WRs difference");
    fh1_wrs[0]->GetYaxis()->SetTitle("Counts");
    fh1_wrs[0]->GetYaxis()->SetTitleOffset(1.3);
    fh1_wrs[0]->GetXaxis()->CenterTitle(true);
    fh1_wrs[0]->GetYaxis()->CenterTitle(true);
    fh1_wrs[0]->SetLineColor(2);
    fh1_wrs[0]->SetLineWidth(3);
    gPad->SetLogy(1);
    fh1_wrs[0]->Draw("");
    fh1_wrs[1] = new TH1F("fh1_WR_Sofia_Califa_Messel", "", 1200, -4100, 4100);
    fh1_wrs[1]->SetLineColor(4);
    fh1_wrs[1]->SetLineWidth(3);
    if (fWRItemsCalifa)
        fh1_wrs[1]->Draw("same");
    fh1_wrs[2] = new TH1F("fh1_WR_Sofia_Neuland", "", 1200, -4100, 4100);
    fh1_wrs[2]->SetLineColor(3);
    fh1_wrs[2]->SetLineWidth(3);
    if (fWRItemsNeuland)
        fh1_wrs[2]->Draw("same");
    fh1_wrs[3] = new TH1F("fh1_WR_Sofia_S2", "", 1200, -4100, 4100);
    fh1_wrs[3]->SetLineColor(1);
    fh1_wrs[3]->SetLineWidth(3);
    if (fWRItemsS2)
        fh1_wrs[3]->Draw("same");
    fh1_wrs[4] = new TH1F("fh1_WR_Sofia_S8", "", 1200, -4100, 4100);
    fh1_wrs[4]->SetLineColor(5);
    fh1_wrs[4]->SetLineWidth(3);
    if (fWRItemsS8)
        fh1_wrs[4]->Draw("same");

    TLegend* leg = new TLegend(0.05, 0.9, 0.39, 0.9999, NULL, "brNDC");
    leg->SetBorderSize(0);
    leg->SetTextFont(62);
    leg->SetTextSize(0.03);
    leg->SetLineColor(1);
    leg->SetLineStyle(1);
    leg->SetLineWidth(1);
    leg->SetFillColor(0);
    leg->SetFillStyle(0);
    TLegendEntry* entry = leg->AddEntry("null", "Califa_Wixhausen", "l");
    entry->SetLineColor(4);
    entry->SetLineStyle(1);
    entry->SetLineWidth(3);
    entry->SetTextFont(62);
    entry = leg->AddEntry("null", "Califa_Messel", "l");
    entry->SetLineColor(2);
    entry->SetLineStyle(1);
    entry->SetLineWidth(3);
    entry->SetTextFont(62);
    if (fWRItemsNeuland)
    {
        entry = leg->AddEntry("null", "Neuland", "l");
        entry->SetLineColor(3);
        entry->SetLineStyle(1);
        entry->SetLineWidth(3);
        entry->SetTextFont(62);
    }
    if (fWRItemsS2)
    {
        entry = leg->AddEntry("null", "S2", "l");
        entry->SetLineColor(1);
        entry->SetLineStyle(1);
        entry->SetLineWidth(3);
        entry->SetTextFont(62);
    }
    if (fWRItemsS8)
    {
        entry = leg->AddEntry("null", "S8", "l");
        entry->SetLineColor(5);
        entry->SetLineStyle(1);
        entry->SetLineWidth(3);
        entry->SetTextFont(62);
    }
    leg->Draw();

    // MAIN FOLDER-R3B
    TFolder* mainfol = new TFolder("R3B_General", "R3B WhiteRabbit and trigger info");
    mainfol->Add(cTrigger);
    mainfol->Add(cTpat);
    if (fWRItemsMaster && fWRItemsSofia)
        mainfol->Add(cWr);
    if (fWRItemsSofia && fWRItemsCalifa)
        mainfol->Add(cWrs);
    run->AddObject(mainfol);

    // Register command to reset histograms
    run->GetHttpServer()->RegisterCommand("Reset_GENERAL_HIST", Form("/Objects/%s/->Reset_GENERAL_Histo()", GetName()));

    return kSUCCESS;
}

void R3BGeneralOnlineSpectra::Reset_GENERAL_Histo()
{
    R3BLOG(info, "");
    fh1_tpat->Reset();
    fh2_tpat->Reset();
    fh1_trigger->Reset();
    if (fWRItemsMaster && fWRItemsSofia)
    {
        fh1_wr[0]->Reset();
        fh1_wr[1]->Reset();
    }
    if (fWRItemsCalifa && fWRItemsSofia)
    {
        fh1_wrs[0]->Reset();
        fh1_wrs[1]->Reset();
        if (fWRItemsNeuland)
            fh1_wrs[2]->Reset();
        if (fWRItemsS2)
            fh1_wrs[3]->Reset();
        if (fWRItemsS8)
            fh1_wrs[4]->Reset();
    }
    // Reset Mwpc0 histograms if they exist somewhere
    if (fMwpc0Online)
        fMwpc0Online->Reset_Histo();
    // Reset Mwpc0_1 histograms if they exist somewhere
    if (fMwpc01Online)
        fMwpc01Online->Reset_Histo();
    // Reset Mwpc0_2 histograms if they exist somewhere
    if (fMwpc02Online)
        fMwpc02Online->Reset_Histo();
    // Reset Mwpc1_2 histograms if they exist somewhere
    if (fMwpc12Online)
        fMwpc12Online->Reset_Histo();
    // Reset Mwpc1 histograms if they exist somewhere
    if (fMwpc1Online)
        fMwpc1Online->Reset_Histo();
    // Reset Mwpc2 histograms if they exist somewhere
    if (fMwpc2Online)
        fMwpc2Online->Reset_Histo();
    // Reset Twim histograms if they exist somewhere
    if (fTwimOnline)
        fTwimOnline->Reset_Histo();
    // Reset Music histograms if they exist somewhere
    if (fMusicOnline)
        fMusicOnline->Reset_Histo();
    // Reset AMS histograms if they exist somewhere
    if (fAmsOnline)
        fAmsOnline->Reset_AMS_Histo();
    // Reset Califa histograms if they exist somewhere
    if (fCalifaOnline)
        fCalifaOnline->Reset_CALIFA_Histo();
    // Reset FOOT histograms if they exist somewhere
    if (fFootOnline)
        fFootOnline->Reset_FOOT_Histo();
    // Reset LOS histograms if they exist somewhere
    if (fLosOnline)
        fLosOnline->Reset_LOS_Histo();
    // Reset Incoming Tracking histograms if they exist somewhere
    if (fIncomingTrackingOnline)
        fIncomingTrackingOnline->Reset_Histo();
    // Reset Twim vs Foot histograms if they exist somewhere
    if (fTwimvsFootOnlineSpectra)
        fTwimvsFootOnlineSpectra->Reset_Histo();
    // Reset TofD histograms if they exist somewhere
    if (fTofdOnlineSpectra)
        fTofdOnlineSpectra->Reset_Histo();
}

void R3BGeneralOnlineSpectra::Exec(Option_t* option)
{
    // Fill histogram with trigger information
    if (fEventHeader->GetTpat() > 0)
    {
        for (Int_t itpat = 0; itpat < 16; itpat++)
        {
            if ((fEventHeader->GetTpat() & (1 << itpat)) == 0)
                continue;
            fh1_tpat->Fill(itpat + 1);
            tpatbin.push_back(itpat + 1);
        }
        if (tpatbin.size() == 1)
        {
            fh2_tpat->Fill(tpatbin.back(), tpatbin.back()); // if only one tpatbit, fill diagonal
        }
        else if (tpatbin.size() > 1)
        {
            auto combination = static_cast<Double_t>(tpatbin.size() * (tpatbin.size() - 1) / 2);
            for (Int_t lower = 0; lower < tpatbin.size() - 1; lower++)
            {
                for (Int_t higher = lower + 1; higher < tpatbin.size(); higher++)
                {
                    fh2_tpat->Fill(tpatbin.at(lower), tpatbin.at(higher), 1. / combination); // Weighed by combination
                }
            }
        }
        else
        {
            R3BLOG(warn, "No tpat found below 16");
        }
    }
    else if (fEventHeader->GetTpat() == 0)
    {
        fh1_tpat->Fill(0.);
        fh2_tpat->Fill(0., 0.);
    }
    else
    {
        LOG(info) << fNEvents << " " << fEventHeader->GetTpat();
    }

    // Trigger
    if (fEventHeader)
        fh1_trigger->Fill(fEventHeader->GetTrigger());

    // WR data
    if (fWRItemsSofia && fWRItemsSofia->GetEntriesFast() > 0)
    {
        // SOFIA
        Int_t nHits = fWRItemsSofia->GetEntriesFast();
        int64_t wrs[2];
        for (Int_t ihit = 0; ihit < nHits; ihit++)
        {
            R3BWRData* hit = dynamic_cast<R3BWRData*>(fWRItemsSofia->At(ihit));
            if (!hit)
                continue;
            wrs[ihit] = hit->GetTimeStamp();
        }

        // Califa
        if (fWRItemsCalifa && fWRItemsCalifa->GetEntriesFast() > 0)
        {
            nHits = fWRItemsCalifa->GetEntriesFast();
            int64_t wr[nHits];
            for (Int_t ihit = 0; ihit < nHits; ihit++)
            {
                R3BWRData* hit = dynamic_cast<R3BWRData*>(fWRItemsCalifa->At(ihit));
                if (!hit)
                    continue;
                wr[ihit] = hit->GetTimeStamp();
            }
            fh1_wrs[0]->Fill(wrs[0] - wr[0]); // messel
            fh1_wrs[1]->Fill(wrs[0] - wr[1]); // wixhausen
        }
        // Neuland
        if (fWRItemsNeuland && fWRItemsNeuland->GetEntriesFast() > 0)
        {
            nHits = fWRItemsNeuland->GetEntriesFast();
            for (Int_t ihit = 0; ihit < nHits; ihit++)
            {
                R3BWRData* hit = dynamic_cast<R3BWRData*>(fWRItemsNeuland->At(ihit));
                if (!hit)
                    continue;
                fh1_wrs[2]->Fill(int64_t(wrs[0] - hit->GetTimeStamp()));
            }
            // fh1_wrs[4]->GetMaximum();
            fh1_wrs[0]->SetMaximum(5. * fh1_wrs[2]->GetBinContent(fh1_wrs[2]->GetMaximumBin()));
        }
        // S2
        if (fWRItemsS2 && fWRItemsS2->GetEntriesFast() > 0)
        {
            nHits = fWRItemsS2->GetEntriesFast();
            for (Int_t ihit = 0; ihit < nHits; ihit++)
            {
                R3BWRData* hit = dynamic_cast<R3BWRData*>(fWRItemsS2->At(ihit));
                if (!hit)
                    continue;
                fh1_wrs[3]->Fill(int64_t(wrs[0] - hit->GetTimeStamp()));
            }
        }
        // S8
        if (fWRItemsS8 && fWRItemsS8->GetEntriesFast() > 0)
        {
            nHits = fWRItemsS8->GetEntriesFast();
            for (Int_t ihit = 0; ihit < nHits; ihit++)
            {
                R3BWRData* hit = dynamic_cast<R3BWRData*>(fWRItemsS8->At(ihit));
                if (!hit)
                    continue;
                fh1_wrs[4]->Fill(int64_t(wrs[0] - hit->GetTimeStamp()));
            }
        }
        // Master
        if (fWRItemsMaster && fWRItemsMaster->GetEntriesFast() > 0)
        {
            nHits = fWRItemsMaster->GetEntriesFast();
            int64_t wrm = 0.;
            for (Int_t ihit = 0; ihit < nHits; ihit++)
            {
                R3BWRData* hit = dynamic_cast<R3BWRData*>(fWRItemsMaster->At(ihit));
                if (!hit)
                    continue;
                wrm = hit->GetTimeStamp();
            }
            fh1_wr[0]->Fill(wrm - wrs[0]);
            fh1_wr[1]->Fill(wrm - wrs[1]);
        }
    }
    fNEvents += 1;
    return;
}

void R3BGeneralOnlineSpectra::FinishEvent()
{
    tpatbin.clear();
    if (fWRItemsMaster)
    {
        fWRItemsMaster->Clear();
    }
    if (fWRItemsSofia)
    {
        fWRItemsSofia->Clear();
    }
    if (fWRItemsCalifa)
    {
        fWRItemsCalifa->Clear();
    }
    if (fWRItemsNeuland)
    {
        fWRItemsNeuland->Clear();
    }
    if (fWRItemsS2)
    {
        fWRItemsS2->Clear();
    }
    if (fWRItemsS8)
    {
        fWRItemsS8->Clear();
    }
}

void R3BGeneralOnlineSpectra::FinishTask()
{
    // Write canvas in the root file
    cTrigger->Write();
    cTpat->Write();
    if (fWRItemsMaster && fWRItemsSofia)
    {
        cWr->Write();
        cWrs->Write();
    }
}

ClassImp(R3BGeneralOnlineSpectra);
