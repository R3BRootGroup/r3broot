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
// -----        R3BMwpcCorrelationOnlineSpectra           -----
// -----    Created 02/12/19 by J.L. Rodriguez-Sanchez    -----
// -----           Fill MWPC online histograms            -----
// ------------------------------------------------------------

/*
 * This task should fill histograms with MWPC online data correlations
 */

#include "R3BMwpcCorrelationOnlineSpectra.h"
#include "R3BEventHeader.h"
#include "R3BLogger.h"
#include "R3BMwpcCalData.h"
#include "R3BMwpcHitData.h"
#include "R3BShared.h"

#include "FairLogger.h"
#include "FairRootManager.h"
#include "FairRunAna.h"
#include "FairRunOnline.h"
#include "FairRuntimeDb.h"

#include <TCanvas.h>
#include <TClonesArray.h>
#include <TFolder.h>
#include <TH1F.h>
#include <TH2F.h>
#include <THttpServer.h>
#include <TMath.h>
#include <TRandom.h>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdlib.h>

using namespace std;

R3BMwpcCorrelationOnlineSpectra::R3BMwpcCorrelationOnlineSpectra()
    : R3BMwpcCorrelationOnlineSpectra("MwpcCorrelationOnlineSpectra", 1, "Mwpc1", "Mwpc2")
{
}

R3BMwpcCorrelationOnlineSpectra::R3BMwpcCorrelationOnlineSpectra(const TString& name,
                                                                 Int_t iVerbose,
                                                                 const TString& namedet1,
                                                                 const TString& namedet2)
    : FairTask(name, iVerbose)
    , fNameDet1(namedet1)
    , fNameDet2(namedet2)
{
}

InitStatus R3BMwpcCorrelationOnlineSpectra::Init()
{
    LOG(info) << "R3B" + fNameDet1 + "vs" + fNameDet2 + "CorrelationOnlineSpectra::Init ";

    // try to get a handle on the EventHeader. EventHeader may not be
    // present though and hence may be null. Take care when using.

    FairRootManager* mgr = FairRootManager::Instance();
    if (NULL == mgr)
        LOG(fatal) << "R3B" + fNameDet1 + "vs" + fNameDet2 + "CorrelationOnlineSpectra::Init FairRootManager not found";

    FairRunOnline* run = FairRunOnline::Instance();
    run->GetHttpServer()->Register("", this);

    // get access to mapped data of the mwpc detectors
    fCalItemsMwpc1 = dynamic_cast<TClonesArray*>(mgr->GetObject(fNameDet1 + "CalData"));
    R3BLOG_IF(fatal, fCalItemsMwpc1 == nullptr, fNameDet1 + "CalData not found");

    fCalItemsMwpc2 = dynamic_cast<TClonesArray*>(mgr->GetObject(fNameDet2 + "CalData"));
    R3BLOG_IF(fatal, fCalItemsMwpc2 == nullptr, fNameDet2 + "CalData not found");

    // get access to hit data of mwpcs
    fHitItemsMwpc1 = dynamic_cast<TClonesArray*>(mgr->GetObject(fNameDet1 + "HitData"));
    if (!fHitItemsMwpc1)
        LOG(warn) << "R3B" + fNameDet1 + "vs" + fNameDet2 + "CorrelationOnlineSpectra: " + fNameDet1 +
                         "HitData not found";

    fHitItemsMwpc2 = dynamic_cast<TClonesArray*>(mgr->GetObject(fNameDet2 + "HitData"));
    if (!fHitItemsMwpc2)
        LOG(warn) << "R3B" + fNameDet1 + "vs" + fNameDet2 + "CorrelationOnlineSpectra: " + fNameDet1 +
                         "HitData not found";

    // Create histograms for detectors
    TString Name1;
    TString Name2;

    cMWPCCal2D =
        new TCanvas(fNameDet1 + "_" + fNameDet2 + "_cal", fNameDet1 + "-" + fNameDet2 + " cal info", 10, 10, 800, 700);
    cMWPCCal2D->Divide(2, 1);

    // MWPC: Cal data
    Name1 = "fh2_" + fNameDet1 + "_" + fNameDet2 + "_calx";
    Name2 = "Position-X: " + fNameDet1 + " vs " + fNameDet2;
    if (fNameDet2 != "Mwpc3")
        fh2_mwpc_calx = R3B::root_owned<TH2F>(Name1, Name2, 64, 0.5, 64.5, 64, 0.5, 64.5);
    else
        fh2_mwpc_calx = R3B::root_owned<TH2F>(Name1, Name2, 64, 0.5, 64.5, 288, 0.5, 288.5);
    fh2_mwpc_calx->GetXaxis()->SetTitle(fNameDet1 + "-X [pads]");
    fh2_mwpc_calx->GetYaxis()->SetTitle(fNameDet2 + "-X [pads]");
    fh2_mwpc_calx->GetYaxis()->SetTitleOffset(1.1);
    fh2_mwpc_calx->GetXaxis()->CenterTitle(true);
    fh2_mwpc_calx->GetYaxis()->CenterTitle(true);
    fh2_mwpc_calx->GetXaxis()->SetLabelSize(0.045);
    fh2_mwpc_calx->GetXaxis()->SetTitleSize(0.045);
    fh2_mwpc_calx->GetYaxis()->SetLabelSize(0.045);
    fh2_mwpc_calx->GetYaxis()->SetTitleSize(0.045);
    cMWPCCal2D->cd(1);
    fh2_mwpc_calx->Draw("col");

    Name1 = "fh2_" + fNameDet1 + "_" + fNameDet2 + "_caly";
    Name2 = "Position-Y: " + fNameDet1 + " vs " + fNameDet2;
    if (fNameDet2 != "Mwpc3")
        fh2_mwpc_caly = R3B::root_owned<TH2F>(Name1, Name2, 64, 0.5, 64.5, 64, 0.5, 64.5);
    else
        fh2_mwpc_caly = R3B::root_owned<TH2F>(Name1, Name2, 64, 0.5, 64.5, 120, 0.5, 120.5);
    fh2_mwpc_caly->GetXaxis()->SetTitle(fNameDet1 + "-Y [pads]");
    fh2_mwpc_caly->GetYaxis()->SetTitle(fNameDet2 + "-Y [pads]");
    fh2_mwpc_caly->GetYaxis()->SetTitleOffset(1.1);
    fh2_mwpc_caly->GetXaxis()->CenterTitle(true);
    fh2_mwpc_caly->GetYaxis()->CenterTitle(true);
    fh2_mwpc_caly->GetXaxis()->SetLabelSize(0.045);
    fh2_mwpc_caly->GetXaxis()->SetTitleSize(0.045);
    fh2_mwpc_caly->GetYaxis()->SetLabelSize(0.045);
    fh2_mwpc_caly->GetYaxis()->SetTitleSize(0.045);
    cMWPCCal2D->cd(2);
    fh2_mwpc_caly->Draw("col");

    cMWPCHit2D =
        new TCanvas(fNameDet1 + "_" + fNameDet2 + "_hit", fNameDet1 + "-" + fNameDet2 + " hit info", 10, 10, 800, 700);
    cMWPCHit2D->Divide(2, 1);

    // MWPC: Hit data
    Name1 = "fh2_" + fNameDet1 + "-" + fNameDet2 + "_hitx";
    Name2 = "Position-X: " + fNameDet1 + " vs " + fNameDet2;
    if (fNameDet2 != "Mwpc3")
        fh2_mwpc_hitx = R3B::root_owned<TH2F>(Name1, Name2, 500, -100., 100., 500, -100., 100.);
    else
        fh2_mwpc_hitx = R3B::root_owned<TH2F>(Name1, Name2, 600, -100., 100., 1200, -450, 450);
    fh2_mwpc_hitx->GetXaxis()->SetTitle(fNameDet1 + "-X [mm]");
    fh2_mwpc_hitx->GetYaxis()->SetTitle(fNameDet2 + "-X [mm]");
    fh2_mwpc_hitx->GetYaxis()->SetTitleOffset(1.1);
    fh2_mwpc_hitx->GetXaxis()->CenterTitle(true);
    fh2_mwpc_hitx->GetYaxis()->CenterTitle(true);
    fh2_mwpc_hitx->GetXaxis()->SetLabelSize(0.045);
    fh2_mwpc_hitx->GetXaxis()->SetTitleSize(0.045);
    fh2_mwpc_hitx->GetYaxis()->SetLabelSize(0.045);
    fh2_mwpc_hitx->GetYaxis()->SetTitleSize(0.045);
    cMWPCHit2D->cd(1);
    fh2_mwpc_hitx->Draw("col");

    Name1 = "fh2_" + fNameDet1 + "_" + fNameDet2 + "_hity";
    Name2 = "Position-Y: " + fNameDet1 + " vs " + fNameDet2;
    if (fNameDet2 != "Mwpc3")
        fh2_mwpc_hity = R3B::root_owned<TH2F>(Name1, Name2, 500, -100., 100., 500, -100., 100.);
    else
        fh2_mwpc_hity = R3B::root_owned<TH2F>(Name1, Name2, 600, -100., 100., 800, -300., 300.);
    fh2_mwpc_hity->GetXaxis()->SetTitle(fNameDet1 + "-Y [mm]");
    fh2_mwpc_hity->GetYaxis()->SetTitle(fNameDet2 + "-Y [mm]");
    fh2_mwpc_hity->GetYaxis()->SetTitleOffset(1.1);
    fh2_mwpc_hity->GetXaxis()->CenterTitle(true);
    fh2_mwpc_hity->GetYaxis()->CenterTitle(true);
    fh2_mwpc_hity->GetXaxis()->SetLabelSize(0.045);
    fh2_mwpc_hity->GetXaxis()->SetTitleSize(0.045);
    fh2_mwpc_hity->GetYaxis()->SetLabelSize(0.045);
    fh2_mwpc_hity->GetYaxis()->SetTitleSize(0.045);
    cMWPCHit2D->cd(2);
    fh2_mwpc_hity->Draw("col");

    // MAIN FOLDER
    TFolder* mainfolMW = new TFolder(fNameDet1 + "-" + fNameDet2, fNameDet1 + "-" + fNameDet2 + " info");
    mainfolMW->Add(cMWPCCal2D);
    if (fHitItemsMwpc1 && fHitItemsMwpc2)
        mainfolMW->Add(cMWPCHit2D);
    run->AddObject(mainfolMW);

    // Register command to reset histograms
    run->GetHttpServer()->RegisterCommand("Reset_" + fNameDet1 + "_" + fNameDet2 + "_HIST",
                                          Form("/Objects/%s/->Reset_Histo()", GetName()));

    return kSUCCESS;
}

void R3BMwpcCorrelationOnlineSpectra::Reset_Histo()
{
    LOG(info) << "R3B" + fNameDet1 + "vs" + fNameDet2 + "CorrelationOnlineSpectra::Reset_Histo";
    // Cal data
    fh2_mwpc_calx->Reset();
    fh2_mwpc_caly->Reset();

    // Hit data
    if (fHitItemsMwpc1 && fHitItemsMwpc2)
    {
        fh2_mwpc_hitx->Reset();
        fh2_mwpc_hity->Reset();
    }
}

void R3BMwpcCorrelationOnlineSpectra::Exec(Option_t* /*option*/)
{
    FairRootManager* mgr = FairRootManager::Instance();
    if (NULL == mgr)
        LOG(fatal) << "R3B" + fNameDet1 + "vs" + fNameDet2 + "CorrelationOnlineSpectra::Exec FairRootManager not found";

    // Fill Cal data
    if (fCalItemsMwpc1 && fCalItemsMwpc1->GetEntriesFast() > 0 && fCalItemsMwpc2 &&
        fCalItemsMwpc2->GetEntriesFast() > 0)
    {
        // First MWPC
        Int_t nHits = fCalItemsMwpc1->GetEntriesFast();
        Int_t maxpadx1 = -1, maxpady1 = -1, maxqx1 = 0, maxqy1 = 0;
        for (Int_t ihit = 0; ihit < nHits; ihit++)
        {
            R3BMwpcCalData* hit = dynamic_cast<R3BMwpcCalData*>(fCalItemsMwpc1->At(ihit));
            if (!hit)
                continue;
            if (hit->GetPlane() == 1 || hit->GetPlane() == 2)
            {
                if (hit->GetQ() > maxqx1)
                {
                    maxpadx1 = hit->GetPad();
                    maxqx1 = hit->GetQ();
                }
            }
            if (hit->GetPlane() == 3)
            {
                if (hit->GetQ() > maxqy1)
                {
                    maxpady1 = hit->GetPad();
                    maxqy1 = hit->GetQ();
                }
            }
        }
        // Second MWPC
        nHits = fCalItemsMwpc2->GetEntriesFast();
        Int_t maxpadx2 = -1, maxpady2 = -1, maxqx2 = 0, maxqy2 = 0;
        for (Int_t ihit = 0; ihit < nHits; ihit++)
        {
            R3BMwpcCalData* hit = dynamic_cast<R3BMwpcCalData*>(fCalItemsMwpc2->At(ihit));
            if (!hit)
                continue;
            if (hit->GetPlane() == 1 || hit->GetPlane() == 2)
            {
                if (hit->GetQ() > maxqx2)
                {
                    maxpadx2 = hit->GetPad();
                    maxqx2 = hit->GetQ();
                }
            }
            if (hit->GetPlane() == 3)
            {
                if (hit->GetQ() > maxqy2)
                {
                    maxpady2 = hit->GetPad();
                    maxqy2 = hit->GetQ();
                }
            }
        }
        if (maxpadx1 > -1 && maxpadx2 > -1)
            fh2_mwpc_calx->Fill(maxpadx1 + gRandom->Uniform(-0.5, 0.5), maxpadx2 + gRandom->Uniform(-0.5, 0.5));
        if (maxpady1 > -1 && maxpady2 > -1)
            fh2_mwpc_caly->Fill(maxpady1 + gRandom->Uniform(-0.5, 0.5), maxpady2 + gRandom->Uniform(-0.5, 0.5));
    }

    // Fill Hit data
    if (fHitItemsMwpc1 && fHitItemsMwpc1->GetEntriesFast() > 0 && fHitItemsMwpc2 &&
        fHitItemsMwpc2->GetEntriesFast() > 0)
    {
        Double_t mw1x = -500., mw1y = -500., mw2x = -500., mw2y = -500.;
        Int_t nHits1 = fHitItemsMwpc1->GetEntriesFast();
        for (Int_t ihit = 0; ihit < nHits1; ihit++)
        {
            R3BMwpcHitData* hit = dynamic_cast<R3BMwpcHitData*>(fHitItemsMwpc1->At(ihit));
            if (!hit)
                continue;
            mw1x = hit->GetX();
            mw1y = hit->GetY();
        }

        Int_t nHits2 = fHitItemsMwpc2->GetEntriesFast();
        for (Int_t ihit = 0; ihit < nHits2; ihit++)
        {
            R3BMwpcHitData* hit = dynamic_cast<R3BMwpcHitData*>(fHitItemsMwpc2->At(ihit));
            if (!hit)
                continue;
            mw2x = hit->GetX();
            mw2y = hit->GetY();
        }
        if (mw1x > -500. && mw1y > -500. && mw2x > -500. && mw2y > -500.)
        {
            fh2_mwpc_hitx->Fill(mw1x, mw2x);
            fh2_mwpc_hity->Fill(mw1y, mw2y);
        }
    }

    fNEvents += 1;
}

void R3BMwpcCorrelationOnlineSpectra::FinishEvent()
{
    if (fCalItemsMwpc1)
    {
        fCalItemsMwpc1->Clear();
    }
    if (fHitItemsMwpc1)
    {
        fHitItemsMwpc1->Clear();
    }
    if (fCalItemsMwpc2)
    {
        fCalItemsMwpc2->Clear();
    }
    if (fHitItemsMwpc2)
    {
        fHitItemsMwpc2->Clear();
    }
}

void R3BMwpcCorrelationOnlineSpectra::FinishTask()
{

    if (fCalItemsMwpc1 && fCalItemsMwpc2)
    {
        cMWPCCal2D->Write();
    }
    if (fHitItemsMwpc1 && fHitItemsMwpc2)
    {
        cMWPCHit2D->Write();
    }
}

ClassImp(R3BMwpcCorrelationOnlineSpectra)
