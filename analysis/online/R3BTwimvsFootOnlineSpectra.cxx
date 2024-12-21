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
// -----           R3BTwimvsFootOnlineSpectra             -----
// -----    Created 24/02/22 by J.L. Rodriguez-Sanchez    -----
// -----  Fill Twim-Foot or Musli-Foot online histograms  -----
// ------------------------------------------------------------

/*
 * This task should fill histograms for Twim-Foot or Musli-Foot correlations
 */

// ROOT headers
#include "TCanvas.h"
#include "TClonesArray.h"
#include "TFolder.h"
#include "TH1F.h"
#include "TH2F.h"
#include "THttpServer.h"
#include "TMath.h"
#include "TVector3.h"

// FAIR headers
#include "FairLogger.h"
#include "FairRootManager.h"
#include "FairRunOnline.h"
#include "FairRuntimeDb.h"

// R3B headers
#include "R3BEventHeader.h"
#include "R3BFootHitData.h"
#include "R3BMusliHitData.h"
#include "R3BTwimHitData.h"
#include "R3BTwimvsFootOnlineSpectra.h"

// R3BTwimvsFootOnlineSpectra::Default Constructor --------------------------
R3BTwimvsFootOnlineSpectra::R3BTwimvsFootOnlineSpectra()
    : R3BTwimvsFootOnlineSpectra("TwimvsFootOnlineSpectra", 1)
{
}

// R3BTwimvsFootOnlineSpectra::Standard Constructor --------------------------
R3BTwimvsFootOnlineSpectra::R3BTwimvsFootOnlineSpectra(const TString& name, Int_t iVerbose)
    : FairTask(name, iVerbose)
    , fHitItemsFoot(NULL)
    , fHitItemsTwim(NULL)
    , fEventHeader(NULL)
    , fNbFootDet(2)
    , fNEvents(0)
    , fMusli(kFALSE)
    , fFootDetId1(1)
    , fFootDetId2(2)
{
}

// Virtual R3BTwimvsFootOnlineSpectra::Destructor
R3BTwimvsFootOnlineSpectra::~R3BTwimvsFootOnlineSpectra()
{
    R3BLOG(debug1, "Destructor");
    if (fHitItemsFoot)
        delete fHitItemsFoot;
    if (fHitItemsTwim)
        delete fHitItemsTwim;
}

InitStatus R3BTwimvsFootOnlineSpectra::Init()
{
    R3BLOG(info, "");

    FairRootManager* mgr = FairRootManager::Instance();
    R3BLOG_IF(fatal, NULL == mgr, "FairRootManager not found");

    FairRunOnline* run = FairRunOnline::Instance();
    run->GetHttpServer()->Register("", this);

    // Look for the R3BEventHeader
    fEventHeader = dynamic_cast<R3BEventHeader*>(mgr->GetObject("EventHeader."));
    if (!fEventHeader)
    {
        R3BLOG(warn, "EventHeader. not found");
        fEventHeader = dynamic_cast<R3BEventHeader*>(mgr->GetObject("R3BEventHeader"));
    }
    else
        R3BLOG(info, "EventHeader. found");

    // get access to hit data of the FOOT detector
    fHitItemsFoot = dynamic_cast<TClonesArray*>(mgr->GetObject("FootHitData"));
    if (!fHitItemsFoot)
    {
        R3BLOG(fatal, "FootHitData not found");
        return kFATAL;
    }

    // get access to hit data of the Musli
    fHitItemsTwim = dynamic_cast<TClonesArray*>(mgr->GetObject("MusliHitData"));
    if (fHitItemsTwim)
    {
        fMusli = kTRUE;
        R3BLOG(info, "MusliHitData found");
    }
    else
    {
        // get access to hit data of the TWIM
        fHitItemsTwim = dynamic_cast<TClonesArray*>(mgr->GetObject("TwimHitData"));
        R3BLOG_IF(warn, !fHitItemsTwim, "TwimHitData not found");
    }

    // MAIN FOLDER-Twim-Foot
    TFolder* mainfolTwim;
    if (fMusli)
    {
        mainfolTwim = new TFolder("MUSLI_vs_FOOT", "Musli vs FOOT info");
    }
    else
    {
        mainfolTwim = new TFolder("TWIM_vs_FOOT", "TWIM vs FOOT info");
    }

    // Create histograms for detectors
    char Name1[255];
    char Name2[255];
    // Energy range for strips
    Double_t binsE = 1000;
    Double_t minE = 0;
    Double_t maxE = 8000;

    // Energy correlations
    fh2_hit_e.resize(fNbFootDet);
    for (Int_t i = 0; i < fNbFootDet; i++)
    { // one histo per detector
        if (fMusli)
        {
            sprintf(Name1, "fh2_energy_musli_vs_foot_%d", i + 1);
            sprintf(Name2, "Energy Musli vs energy Foot Det: %d", i + 1);
        }
        else
        {
            sprintf(Name1, "fh2_energy_twim_vs_foot_%d", i + 1);
            sprintf(Name2, "Energy Twim vs energy Foot Det: %d", i + 1);
        }
        fh2_hit_e[i] = new TH2F(Name1, Name2, binsE, 0, maxE, binsE, 0, 10);
        fh2_hit_e[i]->GetXaxis()->SetTitle("Energy Foot [channels]");
        if (fMusli)
        {
            fh2_hit_e[i]->GetYaxis()->SetTitle("Energy musli [charge-Z]");
        }
        else
        {
            fh2_hit_e[i]->GetYaxis()->SetTitle("Energy twim [charge-Z]");
        }
        fh2_hit_e[i]->GetYaxis()->SetTitleOffset(1.4);
        fh2_hit_e[i]->GetXaxis()->CenterTitle(true);
        fh2_hit_e[i]->GetYaxis()->CenterTitle(true);
        fh2_hit_e[i]->Draw("zcol");
        if (fHitItemsTwim && fHitItemsFoot)
            mainfolTwim->Add(fh2_hit_e[i]);
    }

    // Position correlations
    fh2_hit_x.resize(fNbFootDet);
    for (Int_t i = 0; i < fNbFootDet; i++)
    { // one histo per detector
        if (fMusli)
        {
            sprintf(Name1, "fh2_Xpos_musli_vs_foot_%d", i + 1);
            sprintf(Name2, "X-pos Musli vs Position Foot Det: %d", i + 1);
        }
        else
        {
            sprintf(Name1, "fh2_Xpos_twim_vs_foot_%d", i + 1);
            sprintf(Name2, "X-pos Twim vs Position Foot Det: %d", i + 1);
        }
        fh2_hit_x[i] = new TH2F(Name1, Name2, binsE, -60, 60, binsE, -60, 60);
        fh2_hit_x[i]->GetXaxis()->SetTitle("Position Foot [mm]");
        if (fMusli)
        {
            fh2_hit_x[i]->GetYaxis()->SetTitle("X-pos musli [mm]");
        }
        else
        {
            fh2_hit_x[i]->GetYaxis()->SetTitle("X-pos twim [mm]");
        }
        fh2_hit_x[i]->GetYaxis()->SetTitleOffset(1.4);
        fh2_hit_x[i]->GetXaxis()->CenterTitle(true);
        fh2_hit_x[i]->GetYaxis()->CenterTitle(true);
        fh2_hit_x[i]->Draw("zcol");
        if (fHitItemsTwim && fHitItemsFoot)
            mainfolTwim->Add(fh2_hit_x[i]);
    }

    run->AddObject(mainfolTwim);

    // Register command to reset histograms
    if (fMusli)
    {
        run->GetHttpServer()->RegisterCommand("Reset_MUSLIvsFOOT_HIST", Form("/Objects/%s/->Reset_Histo()", GetName()));
    }
    else
    {
        run->GetHttpServer()->RegisterCommand("Reset_TWIMvsFOOT_HIST", Form("/Objects/%s/->Reset_Histo()", GetName()));
    }

    return kSUCCESS;
}

void R3BTwimvsFootOnlineSpectra::Reset_Histo()
{
    R3BLOG(info, "Reset_Histo");
    if (fHitItemsTwim && fHitItemsFoot)
        for (Int_t i = 0; i < fNbFootDet; i++)
        {
            fh2_hit_e[i]->Reset();
            fh2_hit_x[i]->Reset();
        }
}

void R3BTwimvsFootOnlineSpectra::Exec(Option_t* option)
{
    // Fill hit data
    if (fHitItemsTwim && fHitItemsTwim->GetEntriesFast() > 0 && fHitItemsFoot && fHitItemsFoot->GetEntriesFast() > 0)
    {
        Float_t etwim = 0., xtwim = 0.;
        // TWIM
        Int_t nHits2 = fHitItemsTwim->GetEntriesFast();
        for (Int_t ihit = 0; ihit < nHits2; ihit++)
        {
            if (fMusli)
            {
                auto hit = dynamic_cast<R3BMusliHitData*>(fHitItemsTwim->At(ihit));
                if (!hit)
                    continue;
                etwim = hit->GetZcharge();
                xtwim = hit->GetX();
            }
            else
            {
                auto hit = dynamic_cast<R3BTwimHitData*>(fHitItemsTwim->At(ihit));
                if (!hit)
                    continue;
                etwim = hit->GetZcharge();
                xtwim = hit->GetX();
            }
        }

        Float_t efoot[fNbFootDet], xfoot[fNbFootDet];
        Int_t mulfootdet[fNbFootDet];
        for (Int_t i = 0; i < fNbFootDet; i++)
        {
            efoot[i] = 0.;
            xfoot[i] = -100.;
            mulfootdet[i] = 0;
        }
        // Foot
        Int_t nHits1 = fHitItemsFoot->GetEntriesFast();
        for (Int_t ihit = 0; ihit < nHits1; ihit++)
        {
            auto hit = dynamic_cast<R3BFootHitData*>(fHitItemsFoot->At(ihit));
            if (!hit)
                continue;

            if ((hit->GetDetId() == fFootDetId1) && mulfootdet[hit->GetDetId() - 1] == 0)
            {
                efoot[hit->GetDetId() - 1] = hit->GetEnergy();
                xfoot[hit->GetDetId() - 1] = hit->GetPos();
                mulfootdet[hit->GetDetId() - 1]++;
            }

            if ((hit->GetDetId() == fFootDetId2) && mulfootdet[hit->GetDetId() - 1] == 0)
            {
                efoot[hit->GetDetId() - 1] = hit->GetEnergy();
                xfoot[hit->GetDetId() - 1] = hit->GetPos();
                mulfootdet[hit->GetDetId() - 1]++;
            }
        }
        // Fill histograms
        for (Int_t i = 0; i < fNbFootDet; i++)
        {
            if (efoot[i] > 0.)
                fh2_hit_e[i]->Fill(efoot[i], etwim);
            if (xfoot[i] > -100.)
                fh2_hit_x[i]->Fill(xfoot[i], xtwim);
        }
    }
    fNEvents += 1;
}

void R3BTwimvsFootOnlineSpectra::FinishEvent()
{
    if (fHitItemsFoot)
    {
        fHitItemsFoot->Clear();
    }
    if (fHitItemsTwim)
    {
        fHitItemsTwim->Clear();
    }
}

void R3BTwimvsFootOnlineSpectra::FinishTask()
{
    if (fHitItemsTwim && fHitItemsFoot)
        for (Int_t i = 0; i < fNbFootDet; i++)
        {
            fh2_hit_e[i]->Write();
            fh2_hit_x[i]->Write();
        }
}

ClassImp(R3BTwimvsFootOnlineSpectra);
