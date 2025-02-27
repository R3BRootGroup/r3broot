/******************************************************************************
 *   Copyright (C) 2025 GSI Helmholtzzentrum für Schwerionenforschung GmbH    *
 *   Copyright (C) 2025 Members of R3B Collaboration                          *
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
// -----             R3BActafOnlineSpectra                -----
// -----    Created 27/02/25 by J.L. Rodriguez-Sanchez    -----
// -----          Fill Actaf online histograms            -----
// ------------------------------------------------------------

// ROOT headers
#include <TCanvas.h>
#include <TClonesArray.h>
#include <TFolder.h>
#include <TGaxis.h>
#include <TH1F.h>
#include <TH2F.h>
#include <THttpServer.h>
#include <TMath.h>

// FAIR headers
#include <FairLogger.h>
#include <FairRootManager.h>
#include <FairRunOnline.h>
#include <FairRuntimeDb.h>

// R3B headers
// #include "R3BActafCalData.h"
// #include "R3BActafHitData.h"
#include "R3BActafMappedData.h"
// #include "R3BActafMappingPar.h"
#include "R3BActafOnlineSpectra.h"
#include "R3BEventHeader.h"
#include "R3BLogger.h"
#include "R3BShared.h"

// R3BActafOnlineSpectra::Default Constructor --------------------------
R3BActafOnlineSpectra::R3BActafOnlineSpectra()
    : R3BActafOnlineSpectra("ActafOnlineSpectra", 1)
{
}

// R3BActafOnlineSpectra::Standard Constructor --------------------------
R3BActafOnlineSpectra::R3BActafOnlineSpectra(const TString& name, Int_t iVerbose)
    : FairTask(name, iVerbose)
{
}

void R3BActafOnlineSpectra::SetParContainers()
{
    // Parameter Container
    FairRuntimeDb* rtdb = FairRuntimeDb::instance();
    R3BLOG_IF(fatal, rtdb == nullptr, "FairRuntimeDb not found");

    // fMap_Par = dynamic_cast<R3BActafMappingPar*>(rtdb->getContainer("ActafMappingPar"));
    // R3BLOG_IF(fatal, fMap_Par == nullptr, "Container ActafMappingPar not found");
}

void R3BActafOnlineSpectra::SetParameter()
{
    //--- Parameter Container ---
}

InitStatus R3BActafOnlineSpectra::Init()
{
    R3BLOG(info, "");
    auto* mgr = FairRootManager::Instance();
    R3BLOG_IF(fatal, nullptr == mgr, "FairRootManager not found");

    header = dynamic_cast<R3BEventHeader*>(mgr->GetObject("EventHeader."));
    R3BLOG_IF(error, header == nullptr, "EventHeader. not found");
    R3BLOG_IF(info, header, "EventHeader. found");

    auto* run = FairRunOnline::Instance();
    run->GetHttpServer()->Register("", this);

    fMappedItems = dynamic_cast<TClonesArray*>(mgr->GetObject("ActafMappedData"));
    R3BLOG_IF(fatal, fMappedItems == nullptr, "ActafMappedData not found");

    // fCalItems = dynamic_cast<TClonesArray*>(mgr->GetObject("ActafCalData"));
    // R3BLOG_IF(warn, fCalItems == nullptr, "ActafCalData not found");

    // fHitItems = dynamic_cast<TClonesArray*>(mgr->GetObject("ActafHitData"));
    // R3BLOG_IF(warn, fHitItems == nullptr, "ActafHitData not found");

    // MAIN FOLDER-Actaf
    auto* mainfol = new TFolder("Actaf", "Actaf info");
    // Folder for mapped data
    auto* mapfol = new TFolder("Map", "Map Actaf info");

    //
    // Create histograms
    //
    SetParameter();

    fh2_RawTraces.resize(fPads);

    for (int adc = 0; adc < fFadcs; adc++)
    {
        std::string nameCanvas = "Pad_" + std::to_string(adc + 1) + "_traces_map";
        auto* cMap = new TCanvas(nameCanvas.c_str(), "mapped info", 10, 10, 500, 500);

        for (int chn = 0; chn < fChn; ++chn)
        {
            auto index = adc * fChn + chn + 1;
            std::string nameHist = "fh2_Pad_" + std::to_string(index) + "_trace";
            std::string titleHist = "Raw trace: Pad " + std::to_string(index);

            fh2_RawTraces[index] =
                R3B::root_owned<TH2F>(nameHist.c_str(), titleHist.c_str(), 1346, 1, 2692, 500, 6000, 13000);

            fh2_RawTraces[index]->GetXaxis()->SetTitle("Time [Channels]");
            fh2_RawTraces[index]->GetYaxis()->SetTitle("A");
            fh2_RawTraces[index]->GetYaxis()->SetTitleOffset(1.1);
            fh2_RawTraces[index]->GetXaxis()->CenterTitle(true);
            fh2_RawTraces[index]->GetYaxis()->CenterTitle(true);

            cMap->cd();
            fh2_RawTraces[index]->Draw("colz");

            mapfol->Add(cMap);
        }
    }
    mainfol->Add(mapfol);

    run->AddObject(mainfol);

    // Register command to reset histograms
    run->GetHttpServer()->RegisterCommand("Reset_Actaf_HIST", Form("/Objects/%s/->Reset_Histo()", GetName()));

    return kSUCCESS;
}

// -----   Public method ReInit   ----------------------------------------------
InitStatus R3BActafOnlineSpectra::ReInit()
{
    SetParContainers();
    SetParameter();
    return kSUCCESS;
}

void R3BActafOnlineSpectra::Reset_Histo()
{
    R3BLOG(info, "");

    if (fMappedItems)
    {
        for (const auto& hist : fh2_RawTraces)
        {
            hist->Reset();
        }
    }

    return;
}

void R3BActafOnlineSpectra::Exec(Option_t* /*option*/)
{
    // Check for requested trigger (Todo: should be done globablly / somewhere else)
    if ((fTrigger >= 0) && (header) && (header->GetTrigger() != fTrigger))
        return;

    if (fTpat1 > 0 && fTpat2 > 0 && (header))
    {
        // fTpat = 1-16; fTpat_bit = 0-15
        Int_t fTpat_bit1 = fTpat1 - 1;
        Int_t fTpat_bit2 = fTpat2 - 1;
        Int_t tpatbin = 0;
        for (int i = 0; i < 16; i++)
        {
            tpatbin = (header->GetTpat() & (1 << i));
            if (tpatbin != 0 && (i < fTpat_bit1 || i > fTpat_bit2))
            {
                return;
            }
        }
    }

    // Fill mapped data
    if (fMappedItems && fMappedItems->GetEntriesFast() > 0)
    {
        auto nHits = fMappedItems->GetEntriesFast();
        for (int ihit = 0; ihit < nHits; ihit++)
        {
            auto* hit = dynamic_cast<R3BActafMappedData*>(fMappedItems->At(ihit));
            if (!hit)
                continue;
            auto chn = hit->GetChannel() - 1;

            auto vec = hit->GetTrace(); // std::vector
            std::size_t index = 0;
            for (const auto& value : vec)
                fh2_RawTraces[chn]->Fill(index++, value);
        }
    }

    fNEvents++;
    return;
}

void R3BActafOnlineSpectra::FinishEvent()
{
    R3BLOG(debug, "Cleaning data structures");
    if (fMappedItems)
    {
        fMappedItems->Clear();
    }
}

void R3BActafOnlineSpectra::FinishTask()
{
    for (const auto& hist : fh2_RawTraces)
    {
        hist->Write();
    }
}
ClassImp(R3BActafOnlineSpectra)
