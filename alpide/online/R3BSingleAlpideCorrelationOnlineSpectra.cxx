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
// -----      R3BSingleAlpideCorrelationOnlineSpectra     -----
// -----    Created 28/10/22 by J.L. Rodriguez-Sanchez    -----
// ----- Fill Single ALPIDE correlation online histograms -----
// ------------------------------------------------------------

#include "FairLogger.h"
#include "FairRootManager.h"
#include "FairRunOnline.h"
#include "FairRuntimeDb.h"

#include "R3BAlpideHitData.h"
#include "R3BEventHeader.h"
#include "R3BLogger.h"
#include "R3BSingleAlpideCorrelationOnlineSpectra.h"

#include "TCanvas.h"
#include "TClonesArray.h"
#include "TFolder.h"
#include "THttpServer.h"
#include "TMath.h"

R3BSingleAlpideCorrelationOnlineSpectra::R3BSingleAlpideCorrelationOnlineSpectra()
    : R3BSingleAlpideCorrelationOnlineSpectra("SingleAlpideCorrelationOnlineSpectra", 1, "1", 1, "2", 2)
{
}

R3BSingleAlpideCorrelationOnlineSpectra::R3BSingleAlpideCorrelationOnlineSpectra(const TString& name,
                                                                                 Int_t iVerbose,
                                                                                 const TString& namedet1,
                                                                                 const int id1,
                                                                                 const TString& namedet2,
                                                                                 const int id2)
    : FairTask(name, iVerbose)
    , fHitItems(NULL)
    , fNameDet1(namedet1)
    , fId1(id1)
    , fNameDet2(namedet2)
    , fId2(id2)
    , fNEvents(0)
{
}

R3BSingleAlpideCorrelationOnlineSpectra::~R3BSingleAlpideCorrelationOnlineSpectra() { R3BLOG(debug1, ""); }

InitStatus R3BSingleAlpideCorrelationOnlineSpectra::Init()
{
    R3BLOG(info, "");

    FairRootManager* mgr = FairRootManager::Instance();
    R3BLOG_IF(fatal, NULL == mgr, "FairRootManager not found");

    FairRunOnline* run = FairRunOnline::Instance();
    run->GetHttpServer()->Register("", this);

    fHitItems = dynamic_cast<TClonesArray*>(mgr->GetObject("AlpideHitData"));
    if (!fHitItems)
    {
        R3BLOG(fatal, "AlpideHitData not found");
        return kFATAL;
    }

    // Create histograms for detectors
    TString Name1;
    TString Name2;

    cHit2D = new TCanvas(fNameDet1 + "-" + fNameDet2 + "_hit",
                         "Alpide" + fNameDet1 + "-Alpide" + fNameDet2 + " hit info",
                         10,
                         10,
                         800,
                         700);
    cHit2D->Divide(2, 2);

    // X Hit data
    Name1 = "fh2_" + fNameDet1 + "-" + fNameDet2 + "_hitx";
    Name2 = "Position-X: Alpide" + fNameDet1 + " vs Alpide" + fNameDet2;
    fh2_hitx = new TH2F(Name1, Name2, 200, -15, 15, 200, -15, 15);
    fh2_hitx->GetXaxis()->SetTitle("X" + fNameDet1 + " [mm]");
    fh2_hitx->GetYaxis()->SetTitle("X" + fNameDet2 + " [mm]");
    fh2_hitx->GetYaxis()->SetTitleOffset(1.1);
    fh2_hitx->GetXaxis()->CenterTitle(true);
    fh2_hitx->GetYaxis()->CenterTitle(true);
    fh2_hitx->GetXaxis()->SetLabelSize(0.045);
    fh2_hitx->GetXaxis()->SetTitleSize(0.045);
    fh2_hitx->GetYaxis()->SetLabelSize(0.045);
    fh2_hitx->GetYaxis()->SetTitleSize(0.045);
    cHit2D->cd(1);
    fh2_hitx->Draw("col");

    // Y Hit data
    Name1 = "fh2_" + fNameDet1 + "-" + fNameDet2 + "_hity";
    Name2 = "Position-Y: Alpide" + fNameDet1 + " vs Alpide" + fNameDet2;
    fh2_hity = new TH2F(Name1, Name2, 100, -7.5, 7.5, 100, -7.5, 7.5);
    fh2_hity->GetXaxis()->SetTitle("Y" + fNameDet1 + " [mm]");
    fh2_hity->GetYaxis()->SetTitle("Y" + fNameDet2 + " [mm]");
    fh2_hity->GetYaxis()->SetTitleOffset(1.1);
    fh2_hity->GetXaxis()->CenterTitle(true);
    fh2_hity->GetYaxis()->CenterTitle(true);
    fh2_hity->GetXaxis()->SetLabelSize(0.045);
    fh2_hity->GetXaxis()->SetTitleSize(0.045);
    fh2_hity->GetYaxis()->SetLabelSize(0.045);
    fh2_hity->GetYaxis()->SetTitleSize(0.045);
    cHit2D->cd(2);
    fh2_hity->Draw("col");

    // X Hit diff. data
    Name1 = "fh1_" + fNameDet1 + "-" + fNameDet2 + "_hitx_dif";
    Name2 = "Dif. position-X: Alpide" + fNameDet1 + " - Alpide" + fNameDet2;
    fh1_hitxdif = new TH1F(Name1, Name2, 1000, -5, 5);
    fh1_hitxdif->GetXaxis()->SetTitle("X" + fNameDet1 + " - X" + fNameDet2 + " [mm]");
    fh1_hitxdif->GetXaxis()->CenterTitle(true);
    fh1_hitxdif->GetXaxis()->SetLabelSize(0.045);
    fh1_hitxdif->GetXaxis()->SetTitleSize(0.045);
    fh1_hitxdif->GetYaxis()->SetLabelSize(0.045);
    fh1_hitxdif->GetYaxis()->SetTitleSize(0.045);
    fh1_hitxdif->GetYaxis()->SetTitle("Counts");
    fh1_hitxdif->GetYaxis()->SetTitleOffset(1.1);
    cHit2D->cd(3);
    fh1_hitxdif->SetFillColor(31);
    fh1_hitxdif->Draw("hist");

    // Y Hit diff. data
    Name1 = "fh1_" + fNameDet1 + "-" + fNameDet2 + "_hity_dif";
    Name2 = "Dif. position-Y: Alpide" + fNameDet1 + " - Alpide" + fNameDet2;
    fh1_hitydif = new TH1F(Name1, Name2, 1000, -5, 5);
    fh1_hitydif->GetXaxis()->SetTitle("Y" + fNameDet1 + " - Y" + fNameDet2 + " [mm]");
    fh1_hitydif->GetXaxis()->CenterTitle(true);
    fh1_hitydif->GetXaxis()->SetLabelSize(0.045);
    fh1_hitydif->GetXaxis()->SetTitleSize(0.045);
    fh1_hitydif->GetYaxis()->SetLabelSize(0.045);
    fh1_hitydif->GetYaxis()->SetTitleSize(0.045);
    fh1_hitydif->GetYaxis()->SetTitle("Counts");
    fh1_hitydif->GetYaxis()->SetTitleOffset(1.1);
    cHit2D->cd(4);
    fh1_hitydif->SetFillColor(31);
    fh1_hitydif->Draw("hist");

    // MAIN FOLDER
    TFolder* mainfol = new TFolder("Alpide" + fNameDet1 + "-Alpide" + fNameDet2,
                                   "Alpide" + fNameDet1 + "-Alpide" + fNameDet2 + " info");
    if (fHitItems)
    {
        mainfol->Add(cHit2D);
    }
    run->AddObject(mainfol);

    // Register command to reset histograms
    run->GetHttpServer()->RegisterCommand("Reset_Alpide" + fNameDet1 + "_" + fNameDet2 + "_HIST",
                                          Form("/Objects/%s/->Reset_Histo()", GetName()));

    return kSUCCESS;
}

void R3BSingleAlpideCorrelationOnlineSpectra::Reset_Histo()
{
    R3BLOG(info, "");
    fh2_hitx->Reset();
    fh2_hity->Reset();
    fh1_hitxdif->Reset();
    fh1_hitydif->Reset();
}

void R3BSingleAlpideCorrelationOnlineSpectra::Exec(Option_t* option)
{
    // Fill hit data
    if (fHitItems && fHitItems->GetEntriesFast() > 0)
    {
        float x1 = NAN, x2 = NAN, y1 = NAN, y2 = NAN;
        auto nHits = fHitItems->GetEntriesFast();
        int size1 = 0, size2 = 0;
        for (Int_t ihit = 0; ihit < nHits; ihit++)
        {
            auto hit = dynamic_cast<R3BAlpideHitData*>(fHitItems->At(ihit));
            if (!hit)
                continue;
            Int_t senid = hit->GetSensorId();
            if (senid == fId1 && hit->GetClusterSize() > size1)
            {
                size1 = hit->GetClusterSize();
                x1 = hit->GetX();
                y1 = hit->GetY();
            }
            else if (senid == fId2 && hit->GetClusterSize() > size2)
            {
                size2 = hit->GetClusterSize();
                x2 = hit->GetX();
                y2 = hit->GetY();
            }
        }
        if (size1 > 0 && size2 > 0)
        {
            fh2_hitx->Fill(x1, x2);
            fh1_hitxdif->Fill(x1 - x2);
            fh2_hity->Fill(y1, y2);
            fh1_hitydif->Fill(y1 - y2);
        }
    }

    fNEvents++;
    return;
}

void R3BSingleAlpideCorrelationOnlineSpectra::FinishEvent()
{
    if (fHitItems)
    {
        fHitItems->Clear();
    }
}

void R3BSingleAlpideCorrelationOnlineSpectra::FinishTask()
{
    if (fHitItems)
    {
        cHit2D->Write();
    }
}

ClassImp(R3BSingleAlpideCorrelationOnlineSpectra);
