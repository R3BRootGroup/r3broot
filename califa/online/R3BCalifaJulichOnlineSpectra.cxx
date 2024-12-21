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
// -----             R3BCalifaJulichOnlineSpectra                 -----
// -----    Created 16/07/21  by J.L. Rodriguez-Sanchez   -----
// -----          Fill CalifaJulich online histograms             -----
// ------------------------------------------------------------

/*
 * This task should fill histograms with CalifaJulich online data
 */

#include "R3BCalifaJulichOnlineSpectra.h"
#include "R3BAmsMappedData.h"

#include "R3BEventHeader.h"
#include "THttpServer.h"

#include "FairLogger.h"
#include "FairRootManager.h"
#include "FairRunAna.h"
#include "FairRunOnline.h"
#include "FairRuntimeDb.h"
#include "TCanvas.h"
#include "TFolder.h"
#include "TH1F.h"
#include "TH2F.h"

#include "TClonesArray.h"
#include <cstdlib>
#include <ctime>
#include <fstream>
#include <iostream>
#include <sstream>

R3BCalifaJulichOnlineSpectra::R3BCalifaJulichOnlineSpectra()
    : FairTask("CalifaJulichOnlineSpectra", 1)
    , fMappedItemsCalifa(NULL)
    , fMappedItemsSi(NULL)
    // , fCalItems(NULL)
    // , fHitItems(NULL)
    , fTrigger(-1)
    , fNEvents(0)
    , fNbDet(1)
{
}

R3BCalifaJulichOnlineSpectra::R3BCalifaJulichOnlineSpectra(const TString& name, Int_t iVerbose)
    : FairTask(name, iVerbose)
    , fMappedItemsSi(NULL)
    // , fCalItems(NULL)
    // , fHitItems(NULL)
    , fTrigger(-1)
    , fNEvents(0)
    , fNbDet(1)
{
}

R3BCalifaJulichOnlineSpectra::~R3BCalifaJulichOnlineSpectra()
{
    LOG(debug) << "R3BCalifaJulichOnlineSpectra::Delete instance";
    if (fMappedItemsCalifa)
        delete fMappedItemsCalifa;
    if (fMappedItemsSi)
        delete fMappedItemsSi;
    // if (fCalItems)
    //     delete fCalItems;
    // if (fHitItems)
    //     delete fHitItems;
}

InitStatus R3BCalifaJulichOnlineSpectra::Init()
{
    LOG(info) << "R3BCalifaJulichOnlineSpectra::Init()";

    // Looking for FairRootManager
    FairRootManager* mgr = FairRootManager::Instance();
    if (NULL == mgr)
        LOG(fatal) << "R3BCalifaJulichOnlineSpectra::FairRootManager not found";

    // Create histograms for all the detectors

    // Get access to Mapped data
    fMappedItemsSi = dynamic_cast<TClonesArray*>(mgr->GetObject("AmsMappedData"));
    if (!fMappedItemsSi)
    {
        LOG(fatal) << "R3BCalifaJulichOnlineSpectra::AmsMappedData not found";
        return kFATAL;
    }

    // Get access to Mapped data
    fMappedItemsCalifa = dynamic_cast<TClonesArray*>(mgr->GetObject("CalifaMappedData"));
    if (!fMappedItemsCalifa)
    {
        LOG(fatal) << "R3BCalifaJulichOnlineSpectra::CalifaMappedData not found";
        return kFATAL;
    }
    // // Get access to Cal data
    // fCalItems = (TClonesArray*)mgr->GetObject("CalifaJulichSiCalData");
    // if (!fCalItems)
    //     LOG(warn) << "R3BCalifaJulichOnlineSpectra::CalifaJulichSiCalData not found";
    //
    // // Get access to Hit data
    // fHitItems = (TClonesArray*)mgr->GetObject("CalifaJulichSiHitData");
    // if (!fHitItems)
    //     LOG(warn) << "R3BCalifaJulichOnlineSpectra::CalifaJulichSiHitData not found";

    // Energy range for strips
    Double_t binsE = 200;
    Double_t minE = 0;
    Double_t maxE = 3500;

    char Name1[255];
    char Name2[255];

    // MAIN FOLDER-CalifaJulich
    TFolder* mainfolCalifa = new TFolder("Califa", "Califa info");
    TFolder* mainfolSi = new TFolder("Si", "Si info");
    // Folder for mapped data
    TFolder* mapfolCalifa = new TFolder("Mapped", "Califa Mapped");
    TFolder* mapfolSi = new TFolder("Mapped", "Si Mapped");

    // Mapped data
    fh2_EnergyVsStrip.resize(fNbDet);
    for (Int_t i = 0; i < fNbDet; i++)
    { // one histo per detector
        sprintf(Name1, "fh2_energy_vs_strip_det_%d", i + 1);
        sprintf(Name2, "Energy vs strip number for Si Det: %d", i + 1);
        fh2_EnergyVsStrip[i] = new TH2F(Name1, Name2, 64, 1, 65, binsE, minE, maxE);
        fh2_EnergyVsStrip[i]->GetXaxis()->SetTitle("Strip number");
        fh2_EnergyVsStrip[i]->GetYaxis()->SetTitle("Energy [channels]");
        fh2_EnergyVsStrip[i]->GetYaxis()->SetTitleOffset(1.4);
        fh2_EnergyVsStrip[i]->GetXaxis()->CenterTitle(true);
        fh2_EnergyVsStrip[i]->GetYaxis()->CenterTitle(true);
        fh2_EnergyVsStrip[i]->Draw("col");
        mapfolSi->Add(fh2_EnergyVsStrip[i]);
    }
    mainfolSi->Add(mapfolSi);

    // // Folder for cal data
    // if (fCalItems)
    // {
    //     TFolder* calfol = new TFolder("Cal", "Cal CalifaJulich info");
    //     mainfol->Add(calfol);
    // }
    //
    // // Folder for hit data
    // if (fHitItems)
    // {
    //     TFolder* hitfol = new TFolder("Hit", "Hit CalifaJulich info");
    //     mainfol->Add(hitfol);
    // }

    // Looking for FairRunOnline
    FairRunOnline* run = FairRunOnline::Instance();
    run->GetHttpServer()->Register("", this);
    run->AddObject(mainfolSi);
    run->AddObject(mainfolCalifa);

    // Register command to reset histograms
    run->GetHttpServer()->RegisterCommand("Reset_CalifaJulich",
                                          Form("/Objects/%s/->Reset_CalifaJulich_Histo()", GetName()));

    return kSUCCESS;
}

void R3BCalifaJulichOnlineSpectra::Reset_CalifaJulich_Histo()
{
    LOG(info) << "R3BCalifaJulichOnlineSpectra::Reset_CalifaJulich_Histo";

    // Mapped data
    for (Int_t i = 0; i < fNbDet; i++)
    {
        fh2_EnergyVsStrip[i]->Reset();
    }
}

void R3BCalifaJulichOnlineSpectra::Exec(Option_t* option)
{
    // Fill mapped data
    if (fMappedItemsSi && fMappedItemsSi->GetEntriesFast() > 0)
    {
        auto nHits = fMappedItemsSi->GetEntriesFast();
        for (Int_t ihit = 0; ihit < nHits; ihit++)
        {
            R3BAmsMappedData* hit = dynamic_cast<R3BAmsMappedData*>(fMappedItemsSi->At(ihit));
            if (!hit)
                continue;
            fh2_EnergyVsStrip[hit->GetDetectorId() - 1]->Fill(hit->GetStripId(), hit->GetEnergy());
        }
    }

    fNEvents += 1;
}

void R3BCalifaJulichOnlineSpectra::FinishEvent()
{
    if (fMappedItemsSi)
    {
        fMappedItemsSi->Clear();
    }
    if (fMappedItemsCalifa)
    {
        fMappedItemsCalifa->Clear();
    }
    // if (fCalItems)
    // {
    //     fCalItems->Clear();
    // }
    // if (fHitItems)
    // {
    //     fHitItems->Clear();
    // }
}

void R3BCalifaJulichOnlineSpectra::FinishTask()
{
    if (fMappedItemsSi)
    {
        for (Int_t i = 0; i < fNbDet; i++)
            fh2_EnergyVsStrip[i]->Write();
    }
}

ClassImp(R3BCalifaJulichOnlineSpectra);
