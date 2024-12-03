/******************************************************************************
 *   Copyright (C) 2019 GSI Helmholtzzentrum für Schwerionenforschung GmbH    *
 *   Copyright (C) 2019-2024 Members of R3B Collaboration                     *
 *                                                                            *
 *             This software is distributed under the terms of the            *
 *                 GNU General Public Licence (GPL) version 3,                *
 *                    copied verbatim in the file "LICENSE".                  *
 *                                                                            *
 * In applying this license GSI does not waive the privileges and immunities  *
 * granted to it by virtue of its status as an Intergovernmental Organization *
 * or submit itself to any jurisdiction.                                      *
 ******************************************************************************/

#include "R3BNeulandHitMon.h"
#include "FairLogger.h"
#include "FairRootManager.h"
#include "R3BNeulandHit.h"
#include "R3BShared.h"
#include "TDirectory.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TH3D.h"
#include <TFile.h>
#include <algorithm>
#include <iostream>
#include <numeric>
#include <utility>

R3BNeulandHitMon::R3BNeulandHitMon(TString input, TString output, const Option_t* option)
    : FairTask("R3B NeuLAND NeulandHit Monitor")
    , fOutput(std::move(output))
    , fHits(std::move(input))
{
    LOG(info) << "Using R3B NeuLAND NeulandHit Monitor";

    TString opt = option;
    opt.ToUpper();
    if (opt.Contains("3DTRACK"))
    {
        fIs3DTrackEnabled = true;
        LOG(info) << "... with 3D track visualization";
    }
}

InitStatus R3BNeulandHitMon::Init()
{
    fHits.Init();

    if (fIs3DTrackEnabled)
    {
        // XYZ -> ZXY (side view)
        const auto xbinN = 60;
        const auto ybinN = 50;
        const auto zbinN = 50;
        fh3 = R3B::root_owned<TH3D>("hHits", "hHits", xbinN, 1400., 1700., ybinN, -125., 125., zbinN, -125., 125.);
        fh3->SetTitle("NeuLAND Hits");
        fh3->GetXaxis()->SetTitle("Z");
        fh3->GetYaxis()->SetTitle("X");
        fh3->GetZaxis()->SetTitle("Y");
        FairRootManager::Instance()->Register("NeulandHitMon", "Hits in NeuLAND", fh3, kTRUE);
    }

    // define number of bins for histograms
    const auto maxHitNum = 200;
    const auto timeBinN = 30000;
    const auto zDepBinN = 60;
    const auto energyBinN = 100;
    const auto totenergyBinN = 1000;
    const auto posXYBinN = 300;
    const auto velocityBinN = 200;

    hTime = R3B::root_owned<TH1D>("hTime", "Hit time", timeBinN, -1000., 1000.);
    hTimeAdj = R3B::root_owned<TH1D>("hTimeAdj", "Hit Time adjusted for flight path", timeBinN, -1000., 1000.);

    hMult = R3B::root_owned<TH1I>("hMult", "Hit Multiplicity", maxHitNum, 0, maxHitNum);

    hDepth = R3B::root_owned<TH1D>("hDepth", "Maxial penetration depth", zDepBinN, 1400., 1700.);
    hDepthVSForemostEnergy = R3B::root_owned<TH2D>(
        "hDepthVSFrontEnergy", "Depth vs Foremost Energy", zDepBinN, 1400., 1700., energyBinN, 0., 100.);
    hDepthVSSternmostEnergy = R3B::root_owned<TH2D>(
        "hDepthVSSternmostEnergy", "Depth vs Sternmost Energy", zDepBinN, 1400., 1700., energyBinN, 0, 100.);
    hDepthVSEtot =
        R3B::root_owned<TH2D>("hDepthVSEtot", "Depth vs Total Energy", zDepBinN, 1400., 1700., totenergyBinN, 0, 1000.);
    hForemostEnergy = R3B::root_owned<TH1D>("hForemostEnergy", "Foremost energy deposition", energyBinN, 0, 100.);
    hSternmostEnergy = R3B::root_owned<TH1D>("hSternmostEnergy", "Sternmost energy deposition", energyBinN, 0, 100.);
    hEtot = R3B::root_owned<TH1D>("hEtot", "Total Energy", totenergyBinN, 0, 10000.);
    hdeltaEE = R3B::root_owned<TH2D>(
        "hdeltaEE", "Energy in Foremost Plane vs Etot", energyBinN, 0, 2000., energyBinN, 0, 250.);
    hPosVSEnergy = R3B::root_owned<TH2D>(
        "hPosVSEnergy", "Position vs Energy deposition", zDepBinN, 1400., 1700., totenergyBinN, 0, 1000.);
    hBeta = R3B::root_owned<TH1D>("hBeta", "Velocity", velocityBinN, 0., 1.);
    hE = R3B::root_owned<TH1D>("hE", "Hit Energy", energyBinN, 0., 100.);
    hX = R3B::root_owned<TH1D>("hX", "Hit X", posXYBinN, -150., 150.);
    hY = R3B::root_owned<TH1D>("hY", "Hit Y", posXYBinN, -150., 150.);
    hT = R3B::root_owned<TH1D>("hT", "Hit Delta T", timeBinN, -15., -15.);
    hTNeigh = R3B::root_owned<TH1D>("hTNeigh", "Hit Neigh Delta T", timeBinN, -15., -15.);

    return kSUCCESS;
}

void R3BNeulandHitMon::Exec(Option_t* /*option*/)
{
    const auto hits = fHits.Retrieve();

    // checking paddle multihits
    std::map<Int_t, Int_t> paddlenum;
    for (const auto& hit : hits)
    {
        auto result = paddlenum.insert(std::pair<Int_t, Int_t>(hit->GetPaddle(), 1));
        if (!result.second)
        {
            result.first->second++;
        }
    }
    auto max = std::max_element(paddlenum.begin(),
                                paddlenum.end(),
                                [](std::pair<Int_t, Int_t> lhs, std::pair<Int_t, Int_t> rhs)
                                { return (lhs.second < rhs.second); });
    LOG(debug) << "max dupli: " << max->second;

    if (fIs3DTrackEnabled)
    {
        fh3->Reset("ICES");
        for (const auto& hit : hits)
        {
            fh3->Fill(hit->GetPosition().Z(), hit->GetPosition().X(), hit->GetPosition().Y(), hit->GetE());
        }
    }

    hMult->Fill(static_cast<double>(hits.size()));
    for (const auto& hit : hits)
    {
        hPosVSEnergy->Fill(hit->GetPosition().Z(), hit->GetE());
        hTime->Fill(hit->GetT());
        hTimeAdj->Fill(fDistanceToTarget / hit->GetPosition().Mag() * hit->GetT());
        hBeta->Fill(hit->GetBeta());
        hE->Fill(hit->GetE());
        hX->Fill(hit->GetPosition().X());
        hY->Fill(hit->GetPosition().Y());
    }

    for (auto it1 = hits.begin(); it1 != hits.end(); it1++)
    {
        for (auto it2 = it1 + 1; it2 != hits.end(); it2++)
        {
            if (std::abs((*it1)->GetPosition().X() - (*it2)->GetPosition().X()) < 7.5 &&
                std::abs((*it1)->GetPosition().Y() - (*it2)->GetPosition().Y()) < 7.5 &&
                std::abs((*it1)->GetPosition().Z() - (*it2)->GetPosition().Z()) < 7.5)
            {
                hTNeigh->Fill((*it1)->GetT() - (*it2)->GetT());
            }
            hT->Fill((*it1)->GetT() - (*it2)->GetT());
        }
    }

    auto maxDepthHit = std::max_element(hits.begin(),
                                        hits.end(),
                                        [](R3BNeulandHit* one, R3BNeulandHit* another)
                                        { return one->GetPosition().Z() < another->GetPosition().Z(); });
    if (maxDepthHit != hits.end())
    {
        hDepth->Fill((*maxDepthHit)->GetPosition().Z());
        hSternmostEnergy->Fill((*maxDepthHit)->GetE());
        hDepthVSSternmostEnergy->Fill((*maxDepthHit)->GetPosition().Z(), (*maxDepthHit)->GetE());
    }

    auto minDepthHit = std::min_element(hits.begin(),
                                        hits.end(),
                                        [](R3BNeulandHit* one, R3BNeulandHit* another)
                                        { return one->GetPosition().Z() < another->GetPosition().Z(); });
    auto Etot =
        std::accumulate(hits.begin(), hits.end(), 0., [](double init, const auto* hit) { return init + hit->GetE(); });

    if (minDepthHit != hits.end())
    {
        hForemostEnergy->Fill((*minDepthHit)->GetE());
        hDepthVSForemostEnergy->Fill((*maxDepthHit)->GetPosition().Z(), (*minDepthHit)->GetE());
        hdeltaEE->Fill(Etot, (*minDepthHit)->GetE());
    }

    hEtot->Fill(Etot);
    if (maxDepthHit != hits.end())
    {
        hDepthVSEtot->Fill((*maxDepthHit)->GetPosition().Z(), Etot);
    }
}

void R3BNeulandHitMon::Finish()
{
    TDirectory* tmp = gDirectory;
    FairRootManager::Instance()->GetOutFile()->cd();

    gDirectory->mkdir(fOutput);
    gDirectory->cd(fOutput);

    hDepth->Write();
    hMult->Write();
    hTime->Write();
    hTimeAdj->Write();
    hForemostEnergy->Write();
    hSternmostEnergy->Write();
    hDepthVSForemostEnergy->Write();
    hDepthVSSternmostEnergy->Write();
    hEtot->Write();
    hDepthVSEtot->Write();
    hPosVSEnergy->Write();
    hdeltaEE->Write();
    hBeta->Write();
    hE->Write();
    hX->Write();
    hY->Write();
    hT->Write();
    hTNeigh->Write();

    gDirectory = tmp;
}

ClassImp(R3BNeulandHitMon) // NOLINT
