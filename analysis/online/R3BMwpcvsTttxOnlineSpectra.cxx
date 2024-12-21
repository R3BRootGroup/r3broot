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

#include "R3BMwpcvsTttxOnlineSpectra.h"
#include "R3BEventHeader.h"
#include "R3BLogger.h"
#include "R3BMwpcHitData.h"
#include "R3BShared.h"
#include "R3BTttxHitData.h"

#include <FairLogger.h>
#include <FairRootManager.h>
#include <FairRunOnline.h>
#include <FairRuntimeDb.h>

#include <TCanvas.h>
#include <TFolder.h>
#include <TH2F.h>
#include <THttpServer.h>
#include <TVector3.h>

#include <TClonesArray.h>
#include <sstream>

R3BMwpcvsTttxOnlineSpectra::R3BMwpcvsTttxOnlineSpectra()
    : R3BMwpcvsTttxOnlineSpectra("R3BMwpcvsTttxOnlineSpectra", 1)
{
}

R3BMwpcvsTttxOnlineSpectra::R3BMwpcvsTttxOnlineSpectra(const TString& name, int iVerbose)
    : FairTask(name, iVerbose)
{
}

void R3BMwpcvsTttxOnlineSpectra::SetParContainers()
{
    // Parameter Container
    // Reading amsStripCalPar from FairRuntimeDb
    FairRuntimeDb* rtdb = FairRuntimeDb::instance();
    R3BLOG_IF(fatal, !rtdb, "FairRuntimeDb not found");
}

InitStatus R3BMwpcvsTttxOnlineSpectra::Init()
{
    R3BLOG(info, "");
    FairRootManager* mgr = FairRootManager::Instance();
    R3BLOG_IF(fatal, mgr == nullptr, "FairRootManager not found");

    fHeader = dynamic_cast<R3BEventHeader*>(mgr->GetObject("EventHeader."));
    R3BLOG_IF(fatal, fHeader == nullptr, "EventHeader. not found");

    // get access to Hit data
    fHitItemsMw = dynamic_cast<TClonesArray*>(mgr->GetObject("Mwpc1HitData"));
    R3BLOG_IF(fatal, !fHitItemsMw, "Mwpc1HitData not found");

    fHitItemsTttx = dynamic_cast<TClonesArray*>(mgr->GetObject("tttxHitData"));
    R3BLOG_IF(fatal, !fHitItemsTttx, "TttxHitData not found");

    // Create histograms for detectors

    cCMwTx = std::make_unique<TCanvas>("Mwpc_Tttx_Xcor", "Mwpc-X vs Tttx-X", 10, 10, 500, 500);
    cCMwTx->Divide(fNbHist, 1);
    fh2_xcor.resize(fNbHist);

    cCMwTx_mul = std::make_unique<TCanvas>("Mwpc_Tttx_mulcor", "Mwpc vs Tttx multiplicities", 10, 10, 500, 500);
    cCMwTx_mul->Divide(fNbHist, 1);
    fh2_multcor.resize(fNbHist);

    for (int i = 0; i < fNbHist; i++)
    {
        cCMwTx->cd(i + 1);
        gPad->SetLogz();
        std::stringstream ss1;
        ss1 << "fh2_Mwpc_Tttx_Xcor_" << i;
        fh2_xcor[i] = R3B::root_owned<TH2F>(ss1.str().c_str(), ss1.str().c_str(), 400, -100, 100, 67, -50, 50);
        fh2_xcor[i]->GetXaxis()->SetTitle("Mwpc-X [mm]");
        fh2_xcor[i]->GetYaxis()->SetTitle("Tttx-X [mm]");
        fh2_xcor[i]->GetYaxis()->SetTitleOffset(1.2);
        fh2_xcor[i]->GetXaxis()->CenterTitle(true);
        fh2_xcor[i]->GetYaxis()->CenterTitle(true);
        fh2_xcor[i]->Draw("COLZ");
    }

    for (int i = 0; i < fNbHist; i++)
    {
        cCMwTx_mul->cd(i + 1);
        std::stringstream ss1;
        ss1 << "fh2_Mwpc_Tttx_mul_" << i;
        fh2_multcor[i] = R3B::root_owned<TH2F>(ss1.str().c_str(), ss1.str().c_str(), 10, 0, 10, 15, 0, 15);
        fh2_multcor[i]->GetXaxis()->SetTitle("Mwpc-mul");
        fh2_multcor[i]->GetYaxis()->SetTitle("Tttx-mul");
        fh2_multcor[i]->GetYaxis()->SetTitleOffset(1.2);
        fh2_multcor[i]->GetXaxis()->CenterTitle(true);
        fh2_multcor[i]->GetYaxis()->CenterTitle(true);
        fh2_multcor[i]->Draw("COLZ");
    }

    // MAIN FOLDER-Mwpc-Tttx
    auto mainfol = new TFolder("Mwpc_vs_Tttx", "MWPC vs Tttx info");

    if (fHitItemsMw && fHitItemsTttx)
    {
        mainfol->Add(cCMwTx_mul.get());
        mainfol->Add(cCMwTx.get());
    }

    auto run = FairRunOnline::Instance();
    run->GetHttpServer()->Register("", this);

    run->AddObject(mainfol);

    // Register command to reset histograms
    run->GetHttpServer()->RegisterCommand("Reset_Mwpc_vs_Tttx", Form("/Objects/%s/->Reset_Histo()", GetName()));

    return kSUCCESS;
}

// -----   Public method ReInit   ----------------------------------------------
InitStatus R3BMwpcvsTttxOnlineSpectra::ReInit()
{
    SetParContainers();
    // SetParameter();
    return kSUCCESS;
}

void R3BMwpcvsTttxOnlineSpectra::Reset_Histo()
{
    R3BLOG(info, "");
    for (const auto& hist : fh2_xcor)
    {
        hist->Reset();
    }
    for (const auto& hist : fh2_multcor)
    {
        hist->Reset();
    }
}

void R3BMwpcvsTttxOnlineSpectra::Exec(Option_t* /*option*/)
{
    if ((fTpat > 0) && (fHeader) && ((fHeader->GetTpat() & fTpat) != fTpat))
        return;

    int nHitsMw = fHitItemsMw->GetEntriesFast();
    int nHitsTx = fHitItemsTttx->GetEntriesFast();

    if (nHitsMw == 0 || nHitsTx == 0)
        return;

    std::vector<TVector3> VMw;
    for (int ihit = 0; ihit < nHitsMw; ihit++)
    {
        auto hit = dynamic_cast<R3BMwpcHitData*>(fHitItemsMw->At(ihit));
        VMw.push_back(TVector3(hit->GetX(), hit->GetY(), 0));
    }

    std::vector<int> mul_tttx(fNbHist, 0);
    for (int ihit = 0; ihit < nHitsTx; ihit++)
    {
        auto hit = dynamic_cast<R3BTttxHitData*>(fHitItemsTttx->At(ihit));
        auto idet = hit->GetDetID();
        auto posX = hit->GetX();
        if (idet < 0)
            continue;
        // Multiplicity per tttx detector: 0-all
        mul_tttx[idet]++;

        for (const auto& vtemp : VMw)
        {
            if (abs(vtemp.X() - posX) < fMinDist)
                fh2_xcor[idet]->Fill(vtemp.X(), posX);
        }
    }

    for (int ihist = 0; ihist < fNbHist; ihist++)
    {
        if (mul_tttx[ihist] > 0)
            fh2_multcor[ihist]->Fill(nHitsMw, mul_tttx[ihist]);
    }

    fNEvents += 1;
}

void R3BMwpcvsTttxOnlineSpectra::FinishEvent()
{
    if (fHitItemsMw)
    {
        fHitItemsMw->Clear();
    }
    if (fHitItemsTttx)
    {
        fHitItemsTttx->Clear();
    }
}

void R3BMwpcvsTttxOnlineSpectra::FinishTask()
{
    // Write canvas for Hit data
    if (fHitItemsMw && fHitItemsTttx)
    {
        for (const auto& hist : fh2_xcor)
        {
            hist->Write();
        }
        for (const auto& hist : fh2_multcor)
        {
            hist->Write();
        }
    }
}

ClassImp(R3BMwpcvsTttxOnlineSpectra)
