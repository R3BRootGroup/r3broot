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

#include "R3BCalifavsTofDOnlineSpectra.h"
#include "R3BCalifaClusterData.h"
#include "R3BEventHeader.h"
#include "R3BLogger.h"
#include "R3BShared.h"
#include "R3BTofdHitData.h"
#include "R3BWRData.h"

#include <FairLogger.h>
#include <FairRootManager.h>
#include <FairRunOnline.h>
#include <FairRuntimeDb.h>

#include <TFolder.h>
#include <TH1F.h>
#include <TH2F.h>
#include <THttpServer.h>
#include <TVector3.h>

#include <TClonesArray.h>
#include <TMath.h>
#include <TRandom.h>

R3BCalifavsTofDOnlineSpectra::R3BCalifavsTofDOnlineSpectra()
    : R3BCalifavsTofDOnlineSpectra("CALIFAvsTofDOnlineSpectra", 1)
{
}

R3BCalifavsTofDOnlineSpectra::R3BCalifavsTofDOnlineSpectra(const TString& name, Int_t iVerbose)
    : FairTask(name, iVerbose)
{
}

void R3BCalifavsTofDOnlineSpectra::SetParContainers()
{
    // Parameter Container
    // Reading amsStripCalPar from FairRuntimeDb
    FairRuntimeDb* rtdb = FairRuntimeDb::instance();
    R3BLOG_IF(fatal, !rtdb, "FairRuntimeDb not found");
}

InitStatus R3BCalifavsTofDOnlineSpectra::Init()
{
    R3BLOG(info, "");
    FairRootManager* mgr = FairRootManager::Instance();

    R3BLOG_IF(fatal, nullptr == mgr, "FairRootManager not found");

    header = dynamic_cast<R3BEventHeader*>(mgr->GetObject("EventHeader."));

    FairRunOnline* run = FairRunOnline::Instance();
    run->GetHttpServer()->Register("", this);

    // get access to Hit data
    fHitItemsCalifa = dynamic_cast<TClonesArray*>(mgr->GetObject("CalifaClusterData"));
    R3BLOG_IF(fatal, !fHitItemsCalifa, "CalifaClusterData not found");

    fHitItemsTofd = dynamic_cast<TClonesArray*>(mgr->GetObject("TofdHit"));
    R3BLOG_IF(warn, !fHitItemsTofd, "TofdHit not found");

    // Create histograms for detectors

    // CANVAS Theta vs Phi
    cCalifa_angles = new TCanvas("Califa_Theta_vs_Phi_andTofD", "Theta vs Phi", 10, 10, 500, 500);
    cCalifa_angles->Divide(2, 2);
    fh2_Califa_theta_phi.resize(2);

    for (int i = 0; i < 2; i++)
    {
        cCalifa_angles->cd(i + 1);
        char buf[512];
        snprintf(buf, 512, "%s%s", "fh2_CalifaTofd_theta_vs_phi", (i == 0) ? "" : "_withTofD");
        fh2_Califa_theta_phi[i] = R3B::root_owned<TH2F>(buf, buf, 50, 0, 90, 180, -180, 180);
        fh2_Califa_theta_phi[i]->GetXaxis()->SetTitle("Theta [degrees]");
        fh2_Califa_theta_phi[i]->GetYaxis()->SetTitle("Phi [degrees]");
        fh2_Califa_theta_phi[i]->GetYaxis()->SetTitleOffset(1.2);
        fh2_Califa_theta_phi[i]->GetXaxis()->CenterTitle(true);
        fh2_Califa_theta_phi[i]->GetYaxis()->CenterTitle(true);
        fh2_Califa_theta_phi[i]->Draw("COLZ");
    }

    cCalifa_angles->cd(3);
    fh2_Califa_coinTheta = R3B::root_owned<TH2F>(
        "fh2_CalifaTofd_theta_correlations", "Califa theta correlations", 500, 0, 100, 500, 0, 100);
    fh2_Califa_coinTheta->GetXaxis()->SetTitle("Theta [degrees]");
    fh2_Califa_coinTheta->GetYaxis()->SetTitle("Theta [degrees]");
    fh2_Califa_coinTheta->GetYaxis()->SetTitleOffset(1.2);
    fh2_Califa_coinTheta->GetXaxis()->CenterTitle(true);
    fh2_Califa_coinTheta->GetYaxis()->CenterTitle(true);
    fh2_Califa_coinTheta->Draw("COLZ");

    cCalifa_angles->cd(4);
    fh2_Califa_coinPhi = R3B::root_owned<TH2F>(
        "fh2_CalifaTofd_phi_correlations", "Califa phi correlations", 600, -190, 190, 600, -190, 190);
    fh2_Califa_coinPhi->GetXaxis()->SetTitle("Phi [degrees]");
    fh2_Califa_coinPhi->GetYaxis()->SetTitle("Phi [degrees]");
    fh2_Califa_coinPhi->GetYaxis()->SetTitleOffset(1.2);
    fh2_Califa_coinPhi->GetXaxis()->CenterTitle(true);
    fh2_Califa_coinPhi->GetYaxis()->CenterTitle(true);
    fh2_Califa_coinPhi->Draw("COLZ");

    // MAIN FOLDER-Califa
    auto* mainfolCalifa = new TFolder("CALIFAvsTofD", "CALIFA vs TofD info");

    if (fHitItemsCalifa && fHitItemsTofd)
    {
        mainfolCalifa->Add(cCalifa_angles);
    }
    run->AddObject(mainfolCalifa);

    // Register command to reset histograms
    run->GetHttpServer()->RegisterCommand("Reset_CalifavsTofD", Form("/Objects/%s/->Reset_Histo()", GetName()));

    return kSUCCESS;
}

// -----   Public method ReInit   ----------------------------------------------
InitStatus R3BCalifavsTofDOnlineSpectra::ReInit()
{
    SetParContainers();
    // SetParameter();
    return kSUCCESS;
}

void R3BCalifavsTofDOnlineSpectra::Reset_Histo()
{
    R3BLOG(info, "");
    for (const auto& hist : fh2_Califa_theta_phi)
    {
        hist->Reset();
    }
    fh2_Califa_coinPhi->Reset();
    fh2_Califa_coinTheta->Reset();
}

void R3BCalifavsTofDOnlineSpectra::Exec(Option_t* /*option*/)
{
    if ((fTpat >= 0) && (header) && ((header->GetTpat() & fTpat) != fTpat))
        return;

    bool fZminus1 = false;
    for (Int_t ihit = 0; ihit < fHitItemsTofd->GetEntriesFast(); ihit++)
    {
        auto hit = dynamic_cast<R3BTofdHitData*>(fHitItemsTofd->At(ihit));
        if (!hit)
            continue;
        if (hit->GetDetId() == 1 && hit->GetEloss() > (fZselection - 0.5) && hit->GetEloss() < (fZselection + 0.5))
            fZminus1 = true;
    }

    Int_t nHits = fHitItemsCalifa->GetEntriesFast();
    for (Int_t ihit = 0; ihit < nHits; ihit++)
    {
        auto hit = dynamic_cast<R3BCalifaClusterData*>(fHitItemsCalifa->At(ihit));
        if (hit->GetEnergy() < fMinProtonE)
            continue;

        double theta = hit->GetTheta() * TMath::RadToDeg();
        double phi = hit->GetPhi() * TMath::RadToDeg();
        fh2_Califa_theta_phi[0]->Fill(theta, phi); // always
        if (fZminus1)
            fh2_Califa_theta_phi[1]->Fill(theta, phi); // only with TofD
    }

    // Hit data
    if (fHitItemsCalifa && fZminus1 && fHitItemsCalifa->GetEntriesFast() > 0)
    {
        Double_t theta = 0., phi = 0.;
        Double_t califa_theta[nHits];
        Double_t califa_phi[nHits];
        Double_t califa_e[nHits];
        for (Int_t ihit = 0; ihit < nHits; ihit++)
        {
            auto hit = dynamic_cast<R3BCalifaClusterData*>(fHitItemsCalifa->At(ihit));
            if (!hit)
                continue;
            theta = hit->GetTheta() * TMath::RadToDeg();
            phi = hit->GetPhi() * TMath::RadToDeg();
            califa_theta[ihit] = theta;
            califa_phi[ihit] = phi;
            califa_e[ihit] = hit->GetEnergy();
        }

        TVector3 master[2];
        Double_t maxEL = 0., maxER = 0.;
        for (Int_t i1 = 0; i1 < nHits; i1++)
        {

            if (califa_e[i1] > maxER && TMath::Abs(califa_phi[i1]) > 90.) // wixhausen
            {
                master[0].SetMagThetaPhi(1., califa_theta[i1] * TMath::DegToRad(), califa_phi[i1] * TMath::DegToRad());
                maxER = califa_e[i1];
            }
            if (califa_e[i1] > maxEL && TMath::Abs(califa_phi[i1]) < 90.)
            { // messel
                master[1].SetMagThetaPhi(1., califa_theta[i1] * TMath::DegToRad(), califa_phi[i1] * TMath::DegToRad());
                maxEL = califa_e[i1];
            }
        }
        if (maxEL > fMinProtonE && maxER > fMinProtonE &&
            TMath::Abs(master[0].Phi() * TMath::RadToDeg() - master[1].Phi() * TMath::RadToDeg()) > 150.0)
        {
            if (gRandom->Uniform(0., 1.0) < 0.5)
            {
                fh2_Califa_coinTheta->Fill(master[0].Theta() * TMath::RadToDeg(),
                                           master[1].Theta() * TMath::RadToDeg());
                fh2_Califa_coinPhi->Fill(master[0].Phi() * TMath::RadToDeg(), master[1].Phi() * TMath::RadToDeg());
            }
            else
            {
                fh2_Califa_coinTheta->Fill(master[1].Theta() * TMath::RadToDeg(),
                                           master[0].Theta() * TMath::RadToDeg());
                fh2_Califa_coinPhi->Fill(master[1].Phi() * TMath::RadToDeg(), master[0].Phi() * TMath::RadToDeg());
            }
        }
    }

    fNEvents += 1;
}

void R3BCalifavsTofDOnlineSpectra::FinishEvent()
{
    if (fHitItemsCalifa)
    {
        fHitItemsCalifa->Clear();
    }
    if (fHitItemsTofd)
    {
        fHitItemsTofd->Clear();
    }
}

void R3BCalifavsTofDOnlineSpectra::FinishTask()
{
    // Write canvas for Hit data
    if (fHitItemsCalifa)
    {
        for (const auto& hist : fh2_Califa_theta_phi)
        {
            hist->Write();
        }
        fh2_Califa_coinPhi->Write();
        fh2_Califa_coinTheta->Write();
    }
}

ClassImp(R3BCalifavsTofDOnlineSpectra)
