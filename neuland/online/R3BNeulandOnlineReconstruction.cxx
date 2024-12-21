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

#include "R3BNeulandOnlineReconstruction.h"
#include "FairRunOnline.h"
#include "R3BEventHeader.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TH3D.h"
#include "THttpServer.h"
#include "TStyle.h"
#include <FairRootManager.h>
#include <numeric>

/* This function is required to suppress boxes for empty bins - make them transparent.*/
static Double_t gEmptyBinSupressor(const Double_t* x, const Double_t*)
{
    if (x)
    {
        return *x > 0 ? 1. : 0.;
    }
    return 0.;
}

R3BNeulandOnlineReconstruction::R3BNeulandOnlineReconstruction()
    : FairTask("R3BNeulandOnlineReconstruction", 0)
    , fEventHeader(nullptr)
    , fNeulandHits("NeulandHits")
    , fNeulandClusters("NeulandClusters")
{
}

InitStatus R3BNeulandOnlineReconstruction::Init()
{
    auto run = FairRunOnline::Instance();
    if (run)
    {
        run->GetHttpServer()->Register("/Tasks", this);
        run->GetHttpServer()->RegisterCommand("Reset_Neuland_Reco", Form("/Tasks/%s/->ResetHistos()", GetName()));
    }

    auto ioman = FairRootManager::Instance();
    if (ioman == nullptr)
    {
        throw std::runtime_error("R3BNeulandOnlineReconstruction: No FairRootManager");
    }

    fEventHeader = dynamic_cast<R3BEventHeader*>(ioman->GetObject("EventHeader."));
    if (fEventHeader == nullptr)
    {
        throw std::runtime_error("R3BNeulandOnlineReconstruction: No R3BEventHeader");
    }

    fNeulandHits.Init();
    fNeulandClusters.Init();

    gStyle->SetCanvasPreferGL(kTRUE);
    gStyle->SetPalette(kViridis);
    gStyle->SetOptStat(0);

    auto canvasHits = new TCanvas("canvasHits", "Neuland Hits", 10, 10, 850, 850);
    canvasHits->Divide(3, 2);
    canvasHits->cd(1);
    hHitX = new TH1D("hHitX", "Hit X", 300, -150, 150);
    hHitX->Draw();
    canvasHits->cd(2);
    hHitY = new TH1D("hHitY", "Hit Y", 300, -150, 150);
    hHitY->Draw();
    canvasHits->cd(3);
    hHitZ = new TH1D("hHitZ", "Hit Z", 16, fDistanceToTarget, fDistanceToTarget + 16 * 5);
    hHitZ->Draw();
    canvasHits->cd(4);
    hHitT = new TH1D("hHitT", "Hit T", 100000, -10000, 10000);
    hHitTadj = new TH1D("hHitTadj", "Hit Tadj", 100000, -10000, 10000);
    hHitT->Draw();
    hHitTadj->Draw("same");
    canvasHits->cd(5);
    gPad->SetLogy();
    hHitE = new TH1D("hHitE", "Hit E", 1000, 0, 100);
    hHitE->Draw();
    canvasHits->cd(6);
    gPad->SetLogy();
    hHitMult = new TH1D("hHitMult", "Hit Multiplicity", 100, 0, 100);
    hHitMult->Draw();
    canvasHits->cd(0);
    if (run)
    {
        run->AddObject(canvasHits);
    }

    auto canvas3D = new TCanvas("canvas3D", "Neuland 3D Hits & Clusters", 10, 10, 850, 850);
    canvas3D->Divide(2, 2);
    canvas3D->cd(1);
    hHits3D = new TH3D("hHits3D", "Hits3D, nHits > 10", 16, 0, 16 * 5, 50, -125, 125, 50, -125, 125);
    hHits3D->GetXaxis()->SetTitle("Z");
    hHits3D->GetYaxis()->SetTitle("X");
    hHits3D->GetZaxis()->SetTitle("Y");
    auto lofH = hHits3D->GetListOfFunctions();
    lofH->Add(new TF1("TransferFunction", gEmptyBinSupressor, 0., 1000., 0));
    hHits3D->Draw("glcolz");
    canvas3D->cd(2);
    gPad->SetLogy();
    hClusterSize = new TH1D("hClusterSize", "Cluster Size", 100, 0, 100);
    hClusterSize->Draw();
    canvas3D->cd(3);
    hClusters3D = new TH3D("hClusters3D", "Hits 3D, nHits > 10", 16, 0, 16 * 5, 50, -125, 125, 50, -125, 125);
    hClusters3D->GetXaxis()->SetTitle("Z");
    hClusters3D->GetYaxis()->SetTitle("X");
    hClusters3D->GetZaxis()->SetTitle("Y");
    auto lofC = hClusters3D->GetListOfFunctions();
    lofC->Add(new TF1("TransferFunction", gEmptyBinSupressor, 0., 1000., 0));
    hClusters3D->Draw("glcolz");
    canvas3D->cd(4);
    gPad->SetLogy();
    hClusterMult = new TH1D("hClusterMult", "Cluster Multiplicity", 50, 0, 50);
    hClusterMult->Draw();
    canvas3D->cd(0);
    if (run)
    {
        run->AddObject(canvas3D);
    }

    auto canvasCalibr = new TCanvas("canvasCalibr", "canvasCalibr", 10, 10, 850, 850);
    canvasCalibr->Divide(2, 2);
    canvasCalibr->cd(1);
    hEtotVSNClusters = new TH2D("hEtotVSNClusters", "Calibr", 300, 0, 3000, 50, 0, 50);
    hEtotVSNClusters->GetXaxis()->SetTitle("Total Energy [MeV]");
    hEtotVSNClusters->GetYaxis()->SetTitle("Number of Clusters");
    hEtotVSNClusters->Draw("colz");
    canvasCalibr->cd(2);
    hEtot = new TH1D("hEtot", "Etot", 300, 0, 3000);
    hEtot->Draw();
    canvasCalibr->cd(0);
    if (run)
    {
        run->AddObject(canvasCalibr);
    }

    return kSUCCESS;
}

void R3BNeulandOnlineReconstruction::Exec(Option_t*)
{
    if (std::isnan(fEventHeader->GetTStart()))
    {
        return;
    }

    const auto hits = fNeulandHits.Retrieve();
    const auto nHits = hits.size();
    if (hits.empty())
    {
        return;
    }
    hHitMult->Fill(nHits);

    for (const auto& hit : hits)
    {
        hHitX->Fill(hit->GetPosition().X());
        hHitY->Fill(hit->GetPosition().Y());
        hHitZ->Fill(hit->GetPosition().Z());
        hHitT->Fill(hit->GetT());
        hHitTadj->Fill(fDistanceToTarget / hit->GetPosition().Mag() * hit->GetT());
        hHitE->Fill(hit->GetE());
    }

    const auto clusters = fNeulandClusters.Retrieve();
    const auto nClusters = clusters.size();
    hClusterMult->Fill(nClusters);
    for (const auto& cluster : clusters)
    {
        hClusterSize->Fill(cluster->GetSize());
    }

    if (hits.size() > 10)
    {
        hHits3D->Reset("ICES");
        for (const auto& hit : hits)
        {
            hHits3D->Fill(hit->GetPosition().Z() - fDistanceToTarget,
                          hit->GetPosition().X(),
                          hit->GetPosition().Y(),
                          hit->GetE());
        }
        hClusters3D->Reset("ICES");
        for (const auto& cluster : clusters)
        {
            hClusters3D->Fill(cluster->GetPosition().Z() - fDistanceToTarget,
                              cluster->GetPosition().X(),
                              cluster->GetPosition().Y(),
                              cluster->GetE());
        }
    }

    const double etot = std::accumulate(clusters.begin(),
                                        clusters.end(),
                                        0.,
                                        [](const double sum, const R3BNeulandCluster* c) { return sum + c->GetE(); });
    if (etot > 10)
    {
        hEtotVSNClusters->Fill(etot, nClusters);
    }
}

void R3BNeulandOnlineReconstruction::FinishTask()
{
    TDirectory* tmp = gDirectory;
    FairRootManager::Instance()->GetOutFile()->cd();

    gDirectory->mkdir("R3BNeulandOnlineReconstruction");
    gDirectory->cd("R3BNeulandOnlineReconstruction");

    hHitX->Write();
    hHitY->Write();
    hHitZ->Write();
    hHitT->Write();
    hHitTadj->Write();
    hHitE->Write();

    hHitMult->Write();
    hClusterMult->Write();
    hClusterSize->Write();
    hEtot->Write();
    hEtotVSNClusters->Write();

    gDirectory = tmp;
}

void R3BNeulandOnlineReconstruction::ResetHistos()
{
    hHitX->Reset();
    hHitY->Reset();
    hHitZ->Reset();
    hHitT->Reset();
    hHitTadj->Reset();
    hHitE->Reset();

    hHitMult->Reset();
    hClusterMult->Reset();
    hClusterSize->Reset();
    hEtot->Reset();
    hEtotVSNClusters->Reset();
}

ClassImp(R3BNeulandOnlineReconstruction)
