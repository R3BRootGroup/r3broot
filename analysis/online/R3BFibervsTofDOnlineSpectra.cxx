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
// -----      R3BFibervsTofDOnlineSpectra source file        -----
// -----    Created 29/04/22 by J.L. Rodriguez-Sanchez       -----
// ---------------------------------------------------------------

#include "R3BFibervsTofDOnlineSpectra.h"
#include "R3BEventHeader.h"
#include "R3BFiberMAPMTHitData.h"
#include "R3BFiberMappingPar.h"
#include "R3BLogger.h"
#include "R3BTofdHitData.h"

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

R3BFibervsTofDOnlineSpectra::R3BFibervsTofDOnlineSpectra()
    : R3BFibervsTofDOnlineSpectra("Fi30", 1)
{
}

R3BFibervsTofDOnlineSpectra::R3BFibervsTofDOnlineSpectra(const TString name, Int_t iVerbose)
    : FairTask(name + "_TofDCorrelationOnlineSpectra", iVerbose)
    , fName(name)
    , fTrigger(-1)
    , fTpat1(-1)
    , fTpat2(-1)
    , fNEvents(0)
    , fMapPar(NULL)
    , fHitItems(NULL)
    , fHitTofdItems(NULL)
    , fNbfibers(512)
    , fNbTofdPlanes(4)
    , fNbTofdPaddlesPerPlane(44)
{
}

R3BFibervsTofDOnlineSpectra::~R3BFibervsTofDOnlineSpectra()
{
    for (Int_t i = 0; i < fNbTofdPlanes; i++)
    {
        if (fh2_cor_posX[i])
            delete fh2_cor_posX[i];
        if (fh2_cor_posY[i])
            delete fh2_cor_posY[i];
    }
}

void R3BFibervsTofDOnlineSpectra::SetParContainers()
{
    fMapPar = dynamic_cast<R3BFiberMappingPar*>(FairRuntimeDb::instance()->getContainer(fName + "MappingPar"));
    R3BLOG_IF(error, !fMapPar, "Couldn't get " << fName << "MappingPar");
    if (fMapPar)
    {
        fNbfibers = fMapPar->GetNbChannels();
        R3BLOG(info, fName << "MappingPar found with " << fNbfibers << " fibers");
    }
}

InitStatus R3BFibervsTofDOnlineSpectra::ReInit()
{
    SetParContainers();
    return kSUCCESS;
}

InitStatus R3BFibervsTofDOnlineSpectra::Init()
{
    R3BLOG(info, "For firber " << fName);

    FairRootManager* mgr = FairRootManager::Instance();
    R3BLOG_IF(fatal, NULL == mgr, "FairRootManager not found");

    header = dynamic_cast<R3BEventHeader*>(mgr->GetObject("EventHeader."));
    if (!header)
    {
        R3BLOG(warn, "EventHeader. not found");
        header = dynamic_cast<R3BEventHeader*>(mgr->GetObject("R3BEventHeader"));
    }
    else
        R3BLOG(info, " EventHeader. found");

    // uncomment lines below when ucesb avaliable
    FairRunOnline* run = FairRunOnline::Instance();
    run->GetHttpServer()->Register("", this);

    fHitItems = dynamic_cast<TClonesArray*>(mgr->GetObject(fName + "Hit"));
    R3BLOG_IF(warn, NULL == fHitItems, fName + "Hit not found");

    fHitTofdItems = dynamic_cast<TClonesArray*>(mgr->GetObject("TofdHit"));
    R3BLOG_IF(warn, NULL == fHitTofdItems, "TofdHit not found");

    //------------------------------------------------------------------------
    // create histograms
    //------------------------------------------------------------------------

    auto mainfolder = new TFolder(fName + "_TofD_Cor", fName + " vs TofD info");

    auto canv = new TCanvas(fName + "_TofD_HitCor", fName + "_TofD Hit Correlations", 10, 10, 910, 910);

    if (fHitItems && fHitTofdItems)
    {
        fh2_cor_posX.resize(fNbTofdPlanes);
        fh2_cor_posY.resize(fNbTofdPlanes);
        canv->Divide(2, fNbTofdPlanes);
        char str1[256];
        char str2[256];
        for (Int_t i = 0; i < fNbTofdPlanes; i++)
        {
            // Position X
            snprintf(str1, sizeof str1, "_TofDplane%d_Hit_posX", i + 1);
            snprintf(str2, sizeof str2, "_TofDplane%d Hit X position", i + 1);
            fh2_cor_posX[i] =
                new TH2F(fName + str1, fName + str2, fNbfibers, -25.6, 25.6, fNbTofdPaddlesPerPlane, -60., 60.);
            fh2_cor_posX[i]->GetXaxis()->SetTitle(fName + " X position [cm]");
            fh2_cor_posX[i]->GetYaxis()->SetTitle("TofD X position [cm]");
            fh2_cor_posX[i]->GetYaxis()->SetTitleOffset(1.02);
            fh2_cor_posX[i]->GetXaxis()->CenterTitle(true);
            fh2_cor_posX[i]->GetYaxis()->CenterTitle(true);

            // Position Y
            snprintf(str1, sizeof str1, "_TofDplane%d_Hit_posY", i + 1);
            snprintf(str2, sizeof str2, "_TofDplane%d Hit Y position", i + 1);
            fh2_cor_posY[i] = new TH2F(fName + str1, fName + str2, fNbfibers, -25.6, 25.6, 800, -60., 60.);
            fh2_cor_posY[i]->GetXaxis()->SetTitle(fName + " Y position [cm]");
            fh2_cor_posY[i]->GetYaxis()->SetTitle("TofD Y position [cm]");
            fh2_cor_posY[i]->GetYaxis()->SetTitleOffset(1.02);
            fh2_cor_posY[i]->GetXaxis()->CenterTitle(true);
            fh2_cor_posY[i]->GetYaxis()->CenterTitle(true);

            canv->cd(2 * i + 1);
            fh2_cor_posX[i]->Draw("colz");

            canv->cd(2 * (i + 1));
            fh2_cor_posY[i]->Draw("colz");
        }
        mainfolder->Add(canv);
    }

    run->AddObject(mainfolder);
    run->GetHttpServer()->RegisterCommand("Reset_" + fName + "_Cor_TofD",
                                          Form("/Objects/%s/->Reset_Histo()", GetName()));

    return kSUCCESS;
}
void R3BFibervsTofDOnlineSpectra::Reset_Histo()
{
    R3BLOG(info, "For firber " << fName);
    if (fHitItems && fHitTofdItems)
    {
        for (Int_t i = 0; i < fNbTofdPlanes; i++)
        {
            fh2_cor_posX[i]->Reset();
            fh2_cor_posY[i]->Reset();
        }
    }
}
void R3BFibervsTofDOnlineSpectra::Exec(Option_t* option)
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

    if (fHitItems && fHitItems->GetEntriesFast() > 0 && fHitTofdItems && fHitTofdItems->GetEntriesFast() > 0)
    {
        Double_t xpos1 = 0. / 0., xpos2[fNbTofdPlanes];
        Double_t ypos1 = 0. / 0., ypos2[fNbTofdPlanes];
        Double_t tot = 0., tofdq[fNbTofdPlanes];

        for (Int_t i = 0; i < fNbTofdPlanes; i++)
        {
            xpos2[i] = 0. / 0.;
            ypos2[i] = 0. / 0.;
            tofdq[i] = 0.;
        }

        Int_t nHits1 = fHitItems->GetEntriesFast();
        for (Int_t ihit = 0; ihit < nHits1; ihit++)
        {
            auto hit = dynamic_cast<R3BFiberMAPMTHitData*>(fHitItems->At(ihit));
            if (!hit)
                continue;
            // Looking for the maximum
            if (hit->GetEloss() > tot)
            {
                tot = hit->GetEloss();
                xpos1 = hit->GetX();
                ypos1 = hit->GetY();
            }
        }

        Int_t nHits2 = fHitTofdItems->GetEntriesFast();
        for (Int_t ihit = 0; ihit < nHits2; ihit++)
        {
            auto hit = dynamic_cast<R3BTofdHitData*>(fHitTofdItems->At(ihit));
            if (!hit)
                continue;
            // Looking for the maximum
            Int_t iPlane = hit->GetDetId() - 1;
            if (hit->GetEloss() > tofdq[iPlane])
            {
                tofdq[iPlane] = hit->GetEloss();
                xpos2[iPlane] = hit->GetX();
                ypos2[iPlane] = hit->GetY();
            }
        }

        for (Int_t i = 0; i < fNbTofdPlanes; i++)
            if (tot > 0 && tofdq[i] > 0)
            {
                fh2_cor_posX[i]->Fill(xpos1, xpos2[i]);
                fh2_cor_posY[i]->Fill(ypos1, ypos2[i]);
            }
    }
    return;
}

void R3BFibervsTofDOnlineSpectra::FinishEvent()
{
    if (fHitItems)
    {
        fHitItems->Clear();
    }
    if (fHitTofdItems)
    {
        fHitTofdItems->Clear();
    }
}

void R3BFibervsTofDOnlineSpectra::FinishTask()
{
    if (fHitItems && fHitTofdItems)
    {
        for (Int_t i = 0; i < fNbTofdPlanes; i++)
        {
            fh2_cor_posX[i]->Write();
            fh2_cor_posY[i]->Write();
        }
    }
}

ClassImp(R3BFibervsTofDOnlineSpectra);
