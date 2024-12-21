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

#include "R3BTofDvsTttxOnlineSpectra.h"
#include "R3BEventHeader.h"
#include "R3BLogger.h"
#include "R3BShared.h"
#include "R3BTofdHitData.h"
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

R3BTofDvsTttxOnlineSpectra::R3BTofDvsTttxOnlineSpectra()
    : R3BTofDvsTttxOnlineSpectra("R3BTofDvsTttxOnlineSpectra", 1)
{
}

R3BTofDvsTttxOnlineSpectra::R3BTofDvsTttxOnlineSpectra(const TString& name, int iVerbose)
    : FairTask(name, iVerbose)
{
}

void R3BTofDvsTttxOnlineSpectra::SetParContainers()
{
    // Parameter Container
    // Reading amsStripCalPar from FairRuntimeDb
    FairRuntimeDb* rtdb = FairRuntimeDb::instance();
    R3BLOG_IF(fatal, !rtdb, "FairRuntimeDb not found");
}

InitStatus R3BTofDvsTttxOnlineSpectra::Init()
{
    R3BLOG(info, "");
    FairRootManager* mgr = FairRootManager::Instance();
    R3BLOG_IF(fatal, mgr == nullptr, "FairRootManager not found");

    fHeader = dynamic_cast<R3BEventHeader*>(mgr->GetObject("EventHeader."));
    R3BLOG_IF(fatal, fHeader == nullptr, "EventHeader. not found");

    // get access to Hit data
    fHitItemsTofd = dynamic_cast<TClonesArray*>(mgr->GetObject("TofdHit"));
    R3BLOG_IF(fatal, !fHitItemsTofd, "TofdHit not found");

    fHitItemsTttx = dynamic_cast<TClonesArray*>(mgr->GetObject("tttxHitData"));
    R3BLOG_IF(fatal, !fHitItemsTttx, "TttxHitData not found");

    // Create histograms for detectors
    cCTofDTx = std::make_unique<TCanvas>("TofD_Tttx_charges", "TofD vs Tttx-Charges", 10, 10, 500, 500);
    cCTofDTx->Divide(fNbHist_tofd_tttx, 2);
    fh2_charges.resize(fNbHist_tofd_tttx * 2);

    for (int i = 0; i < 2 * fNbHist_tofd_tttx; i++)
    {
        cCTofDTx->cd(i + 1);
        gPad->SetLogz();
        std::stringstream ss1;
        if (i < fNbHist_tofd_tttx)
            ss1 << "fh2_TofD_Tttx_charge_" << i;
        else
            ss1 << "fh2_TofD_Tttx_charge_highest_" << i;
        fh2_charges[i] = R3B::root_owned<TH2F>(ss1.str().c_str(), ss1.str().c_str(), 400, 0, 10, 400, 0, 10);
        fh2_charges[i]->GetXaxis()->SetTitle("TofD-Charge");
        fh2_charges[i]->GetYaxis()->SetTitle("Tttx-Charge");
        fh2_charges[i]->GetYaxis()->SetTitleOffset(1.2);
        fh2_charges[i]->GetXaxis()->CenterTitle(true);
        fh2_charges[i]->GetYaxis()->CenterTitle(true);
        fh2_charges[i]->Draw("COLZ");
    }

    // MAIN FOLDER-TofD-Tttx
    auto mainfol = new TFolder("TofD_vs_Tttx", "TofD vs Tttx info");

    if (fHitItemsTofd && fHitItemsTttx)
    {
        mainfol->Add(cCTofDTx.get());
    }

    auto run = FairRunOnline::Instance();
    run->GetHttpServer()->Register("", this);

    run->AddObject(mainfol);

    // Register command to reset histograms
    run->GetHttpServer()->RegisterCommand("Reset_TofD_vs_Tttx", Form("/Objects/%s/->Reset_Histo()", GetName()));

    return kSUCCESS;
}

// -----   Public method ReInit   ----------------------------------------------
InitStatus R3BTofDvsTttxOnlineSpectra::ReInit()
{
    SetParContainers();
    // SetParameter();
    return kSUCCESS;
}

void R3BTofDvsTttxOnlineSpectra::Reset_Histo()
{
    R3BLOG(info, "");
    for (const auto& hist : fh2_charges)
    {
        hist->Reset();
    }
}

void R3BTofDvsTttxOnlineSpectra::Exec(Option_t* /*option*/)
{
    if ((fTpat > 0) && (fHeader) && ((fHeader->GetTpat() & fTpat) != fTpat))
        return;

    int nHitsTofd = fHitItemsTofd->GetEntriesFast();
    int nHitsTx = fHitItemsTttx->GetEntriesFast();

    if (nHitsTofd == 0 || nHitsTx == 0)
        return;

    double chargetofd = 0.;
    for (int ihit = 0; ihit < nHitsTofd; ihit++)
    {
        auto hit = dynamic_cast<R3BTofdHitData*>(fHitItemsTofd->At(ihit));
        auto plane = hit->GetDetId();
        if (plane != fPlaneTofd)
            continue;
        // Looking for the max. charge
        if (hit->GetEloss() > chargetofd)
            chargetofd = hit->GetEloss();
    }

    double chargetttx[fNbHist_tofd_tttx] = { 0. };
    for (int ihit = 0; ihit < nHitsTx; ihit++)
    {
        auto hit = dynamic_cast<R3BTttxHitData*>(fHitItemsTttx->At(ihit));
        auto idet = hit->GetDetID();
        // if (idet > 0)
        //   continue;
        fh2_charges[idet]->Fill(chargetofd, hit->GetChargeZ());
        // Looking for the max. charge
        if (hit->GetChargeZ() > chargetttx[idet])
            chargetttx[idet] = hit->GetChargeZ();
    }

    for (int i = 0; i < fNbHist_tofd_tttx; i++)
    {
        if (chargetofd > fMinZTofd && chargetofd < fMaxZTofd && chargetttx[i] > 0)
            fh2_charges[i + fNbHist_tofd_tttx]->Fill(chargetofd, chargetttx[i]);
    }

    fNEvents += 1;
}

void R3BTofDvsTttxOnlineSpectra::FinishEvent()
{
    if (fHitItemsTofd)
    {
        fHitItemsTofd->Clear();
    }
    if (fHitItemsTttx)
    {
        fHitItemsTttx->Clear();
    }
}

void R3BTofDvsTttxOnlineSpectra::FinishTask()
{
    // Write canvas for Hit data
    if (fHitItemsTofd && fHitItemsTttx)
    {
        for (const auto& hist : fh2_charges)
        {
            hist->Write();
        }
    }
}

ClassImp(R3BTofDvsTttxOnlineSpectra)
