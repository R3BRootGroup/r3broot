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

// ---------------------------------------------------------------
// -----  R3BFiberMAPMTCorrelationOnlineSpectra source file  -----
// -----    Created 27/04/22 by J.L. Rodriguez-Sanchez       -----
// ---------------------------------------------------------------

#include "R3BFiberMAPMTCorrelationOnlineSpectra.h"
#include "R3BEventHeader.h"
#include "R3BFiberMAPMTHitData.h"
#include "R3BFiberMappingPar.h"
#include "R3BLogger.h"
#include "R3BShared.h"

#include "FairLogger.h"
#include "FairRootManager.h"
#include "FairRunOnline.h"
#include "FairRuntimeDb.h"

#include "TCanvas.h"
#include "TClonesArray.h"
#include "TFolder.h"
#include "TGaxis.h"
#include "TH1F.h"
#include "TH2F.h"
#include "THttpServer.h"
#include "TMath.h"

R3BFiberMAPMTCorrelationOnlineSpectra::R3BFiberMAPMTCorrelationOnlineSpectra()
    : R3BFiberMAPMTCorrelationOnlineSpectra("Fi30", "Fi31", 1)
{
}

R3BFiberMAPMTCorrelationOnlineSpectra::R3BFiberMAPMTCorrelationOnlineSpectra(const TString name1,
                                                                             const TString name2,
                                                                             Int_t iVerbose)
    : FairTask(name1 + "_" + name2 + "CorrelationOnlineSpectra", iVerbose)
    , fName1(name1)
    , fName2(name2)
    , fTrigger(-1)
    , fTpat1(-1)
    , fTpat2(-1)
    , fNEvents(0)
    , fMapPar1(NULL)
    , fMapPar2(NULL)
    , fHitItems1(NULL)
    , fHitItems2(NULL)
    , fNbfibers1(512)
    , fNbfibers2(512)
{
}

R3BFiberMAPMTCorrelationOnlineSpectra::~R3BFiberMAPMTCorrelationOnlineSpectra()
{
    if (fh_Fib_posX)
        delete fh_Fib_posX;
    if (fh_Fib_posY)
        delete fh_Fib_posY;
    if (fh_Fib_tot)
        delete fh_Fib_tot;
}

void R3BFiberMAPMTCorrelationOnlineSpectra::SetParContainers()
{
    fMapPar1 = dynamic_cast<R3BFiberMappingPar*>(FairRuntimeDb::instance()->getContainer(fName1 + "MappingPar"));
    R3BLOG_IF(error, !fMapPar1, "Couldn't get " << fName1 << "MappingPar");
    if (fMapPar1)
    {
        fNbfibers1 = fMapPar1->GetNbChannels();
        R3BLOG(info, fName1 << "MappingPar found with " << fNbfibers1 << " fibers");
    }

    fMapPar2 = dynamic_cast<R3BFiberMappingPar*>(FairRuntimeDb::instance()->getContainer(fName2 + "MappingPar"));
    R3BLOG_IF(error, !fMapPar2, "Couldn't get " << fName2 << "MappingPar");
    if (fMapPar2)
    {
        fNbfibers2 = fMapPar2->GetNbChannels();
        R3BLOG(info, fName2 << "MappingPar found with " << fNbfibers2 << " fibers");
    }
}

InitStatus R3BFiberMAPMTCorrelationOnlineSpectra::ReInit()
{
    SetParContainers();
    return kSUCCESS;
}

InitStatus R3BFiberMAPMTCorrelationOnlineSpectra::Init()
{
    R3BLOG(info, "For firbers " << fName1 << " and " << fName2);

    FairRootManager* mgr = FairRootManager::Instance();
    R3BLOG_IF(fatal, NULL == mgr, "FairRootManager not found");

    header = dynamic_cast<R3BEventHeader*>(mgr->GetObject("EventHeader."));
    R3BLOG_IF(warn, header == nullptr, "EventHeader. not found");
    R3BLOG_IF(info, header, " EventHeader. found");

    // uncomment lines below when ucesb avaliable
    FairRunOnline* run = FairRunOnline::Instance();
    run->GetHttpServer()->Register("", this);

    fHitItems1 = dynamic_cast<TClonesArray*>(mgr->GetObject(fName1 + "Hit"));
    R3BLOG_IF(warn, fHitItems1 == nullptr, fName1 + "Hit not found");

    fHitItems2 = dynamic_cast<TClonesArray*>(mgr->GetObject(fName2 + "Hit"));
    R3BLOG_IF(warn, fHitItems2 == nullptr, fName2 + "Hit not found");

    //------------------------------------------------------------------------
    // create histograms
    //------------------------------------------------------------------------

    auto mainfolder = new TFolder(fName1 + "_" + fName2 + "_Cor", fName1 + " vs " + fName2 + " info");

    auto FibCanvas =
        new TCanvas(fName1 + "_" + fName2 + "_HitCor", fName1 + "_" + fName2 + " Hit Correlations", 10, 10, 910, 910);

    if (fHitItems1 && fHitItems2)
    {
        FibCanvas->Divide(3, 1);

        // Position X
        fh_Fib_posX = R3B::root_owned<TH2F>(fName1 + "_" + fName2 + "Hit_posX",
                                            fName1 + "_" + fName2 + " Hit X position",
                                            fNbfibers1,
                                            -25.6,
                                            25.6,
                                            fNbfibers2,
                                            -25.6,
                                            25.6);
        fh_Fib_posX->GetXaxis()->SetTitle(fName1 + " X position [cm]");
        fh_Fib_posX->GetYaxis()->SetTitle(fName2 + " X position [cm]");
        fh_Fib_posX->GetYaxis()->SetTitleOffset(1.02);
        fh_Fib_posX->GetXaxis()->CenterTitle(true);
        fh_Fib_posX->GetYaxis()->CenterTitle(true);

        // Position Y
        fh_Fib_posY = R3B::root_owned<TH2F>(fName1 + "_" + fName2 + "Hit_posY",
                                            fName1 + "_" + fName2 + " Hit Y position",
                                            fNbfibers1,
                                            -25.6,
                                            25.6,
                                            fNbfibers2,
                                            -25.6,
                                            25.6);
        fh_Fib_posY->GetXaxis()->SetTitle(fName1 + " Y position [cm]");
        fh_Fib_posY->GetYaxis()->SetTitle(fName2 + " Y position [cm]");
        fh_Fib_posY->GetYaxis()->SetTitleOffset(1.02);
        fh_Fib_posY->GetXaxis()->CenterTitle(true);
        fh_Fib_posY->GetYaxis()->CenterTitle(true);

        // ToT
        fh_Fib_tot = R3B::root_owned<TH2F>(
            fName1 + "_" + fName2 + "Hit_ToT", fName1 + "_" + fName2 + " Hit ToT", 200, 0., 100., 200, 0., 100.);
        fh_Fib_tot->GetXaxis()->SetTitle(fName1 + " ToT");
        fh_Fib_tot->GetYaxis()->SetTitle(fName2 + " ToT");
        fh_Fib_tot->GetYaxis()->SetTitleOffset(1.02);
        fh_Fib_tot->GetXaxis()->CenterTitle(true);
        fh_Fib_tot->GetYaxis()->CenterTitle(true);

        FibCanvas->cd(1);
        fh_Fib_posX->Draw("colz");

        FibCanvas->cd(2);
        fh_Fib_posY->Draw("colz");

        FibCanvas->cd(3);
        fh_Fib_tot->Draw("colz");

        mainfolder->Add(FibCanvas);
    }

    run->AddObject(mainfolder);
    run->GetHttpServer()->RegisterCommand("Reset_" + fName1 + "_Cor_" + fName2,
                                          Form("/Objects/%s/->Reset_Histo()", GetName()));

    return kSUCCESS;
}
void R3BFiberMAPMTCorrelationOnlineSpectra::Reset_Histo()
{
    R3BLOG(info, "For firbers " << fName1 << " and " << fName2);
    if (fHitItems1 && fHitItems2)
    {
        fh_Fib_posX->Reset();
        fh_Fib_posY->Reset();
        fh_Fib_tot->Reset();
    }
}
void R3BFiberMAPMTCorrelationOnlineSpectra::Exec(Option_t* option)
{
    fNEvents += 1;

    if ((fTrigger >= 0) && header && (header->GetTrigger() != fTrigger))
        return;

    // fTpat = 1-16; fTpat_bit = 0-15
    if (fTpat1 > -1 && fTpat2 > -1)
    {
        Int_t fTpat_bit1 = fTpat1 - 1;
        Int_t fTpat_bit2 = fTpat2 - 1;
        Int_t tpatbin;
        for (int i = 0; i < 16; i++)
        {
            tpatbin = (header->GetTpat() & (1 << i));
            if (tpatbin != 0 && (i < fTpat_bit1 || i > fTpat_bit2))
            {
                return;
            }
        }
    }

    if (fHitItems1 && fHitItems1->GetEntriesFast() > 0 && fHitItems2 && fHitItems2->GetEntriesFast() > 0)
    {
        Double_t xpos1 = 0. / 0., xpos2 = 0. / 0.;
        Double_t ypos1 = 0. / 0., ypos2 = 0. / 0.;
        Double_t tot1 = 0., tot2 = 0.;
        Int_t nHits1 = fHitItems1->GetEntriesFast();

        for (Int_t ihit = 0; ihit < nHits1; ihit++)
        {
            auto hit = dynamic_cast<R3BFiberMAPMTHitData*>(fHitItems1->At(ihit));
            if (!hit)
                continue;
            // Looking for the maximum
            if (hit->GetEloss() > tot1)
            {
                tot1 = hit->GetEloss();
                xpos1 = hit->GetX();
                ypos1 = hit->GetY();
            }
        }

        Int_t nHits2 = fHitItems2->GetEntriesFast();
        for (Int_t ihit = 0; ihit < nHits2; ihit++)
        {
            auto hit = dynamic_cast<R3BFiberMAPMTHitData*>(fHitItems2->At(ihit));
            if (!hit)
                continue;
            // Looking for the maximum
            if (hit->GetEloss() > tot2)
            {
                tot2 = hit->GetEloss();
                xpos2 = hit->GetX();
                ypos2 = hit->GetY();
            }
        }

        if (tot1 > 0 && tot2 > 0)
        {
            fh_Fib_posX->Fill(xpos1, xpos2);
            fh_Fib_posY->Fill(ypos1, ypos2);
            fh_Fib_tot->Fill(tot1, tot2);
        }
    }
    return;
}

void R3BFiberMAPMTCorrelationOnlineSpectra::FinishEvent()
{
    if (fHitItems1)
    {
        fHitItems1->Clear();
    }
    if (fHitItems2)
    {
        fHitItems2->Clear();
    }
}

void R3BFiberMAPMTCorrelationOnlineSpectra::FinishTask()
{
    if (fHitItems1 && fHitItems2)
    {
        fh_Fib_posX->Write();
        fh_Fib_posY->Write();
        fh_Fib_tot->Write();
    }
}

ClassImp(R3BFiberMAPMTCorrelationOnlineSpectra)
