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

R3BNeulandHitMon::R3BNeulandHitMon(const Option_t* option)
    : FairTask("R3B NeuLAND NeulandHit Monitor")
{
    LOG(info) << "Using R3B NeuLAND NeulandHit Monitor";

    TString opt = option;
    opt.ToUpper();
    if (opt.Contains("3DTRACK"))
    {
        is_3d_track_enabled_ = true;
        LOG(info) << "... with 3D track visualization";
    }
}

auto R3BNeulandHitMon::Init() -> InitStatus
{
    neuland_hits_.init();

    if (is_3d_track_enabled_)
    {
        // XYZ -> ZXY (side view)
        const auto xbinN = 60;
        const auto ybinN = 50;
        const auto zbinN = 50;

        // define additional histogram parameters
        const auto x_bounds_lower = 1400.;
        const auto x_bounds_upper = 1700.;
        const auto yz_bounds_lower = -125.;
        const auto yz_bounds_upper = 125.;

        hist_3_ = data_monitor_.add_hist<TH3D>("hHits",
                                               "hHits",
                                               xbinN,
                                               x_bounds_lower,
                                               x_bounds_upper,
                                               ybinN,
                                               yz_bounds_lower,
                                               yz_bounds_upper,
                                               zbinN,
                                               yz_bounds_lower,
                                               yz_bounds_upper);
        hist_3_->SetTitle("NeuLAND Hits");
        hist_3_->GetXaxis()->SetTitle("Z");
        hist_3_->GetYaxis()->SetTitle("X");
        hist_3_->GetZaxis()->SetTitle("Y");
        FairRootManager::Instance()->Register("NeulandHitMon", "Hits in NeuLAND", hist_3_, kTRUE);
    }

    // define number of bins for histograms
    const auto maxHitNum = 200;
    const auto timeBinN = 30000;
    const auto zDepBinN = 60;
    const auto energyBinN = 500;
    const auto totenergyBinN = 1000;
    const auto posXYBinN = 300;
    const auto velocityBinN = 200;

    // define additional histogram parameters
    const auto z_dep_lower = 1400.;
    const auto z_dep_upper = 1700.;
    const auto total_energy_upper = 10000.;
    const auto pos_xy_lower = -150.;
    const auto pos_xy_upper = 150.;
    const auto time_para = -15.;
    const auto energy_upper_foremost_vs_e_tot_a = 2000.;
    const auto energy_upper_foremost_vs_e_tot_b = 250.;

    hist_time_ = data_monitor_.add_hist<TH1D>("hTime", "Hit time", timeBinN, -1000., 1000.);
    hist_time_adj_ =
        data_monitor_.add_hist<TH1D>("hTimeAdj", "Hit Time adjusted for flight path", timeBinN, -1000., 1000.);

    hist_mult_ = data_monitor_.add_hist<TH1I>("hMult", "Hit Multiplicity", maxHitNum, 0, maxHitNum);

    hist_depth_ =
        data_monitor_.add_hist<TH1D>("hDepth", "Maxial penetration depth", zDepBinN, z_dep_lower, z_dep_upper);
    hist_depth_vs_foremost_energy_ = data_monitor_.add_hist<TH2D>(
        "hDepthVSFrontEnergy", "Depth vs Foremost Energy", zDepBinN, z_dep_lower, z_dep_upper, energyBinN, 0., 100.);
    hist_depth_vs_sternmost_energy_ = data_monitor_.add_hist<TH2D>("hDepthVSSternmostEnergy",
                                                                   "Depth vs Sternmost Energy",
                                                                   zDepBinN,
                                                                   z_dep_lower,
                                                                   z_dep_upper,
                                                                   energyBinN,
                                                                   0,
                                                                   100.);
    hist_depth_vs_energy_tot_ = data_monitor_.add_hist<TH2D>(
        "hDepthVSEtot", "Depth vs Total Energy", zDepBinN, z_dep_lower, z_dep_upper, totenergyBinN, 0, 1000.);
    hist_foremost_energy_ =
        data_monitor_.add_hist<TH1D>("hForemostEnergy", "Foremost energy deposition", energyBinN, 0, 100.);
    hist_sternmost_energy_ =
        data_monitor_.add_hist<TH1D>("hSternmostEnergy", "Sternmost energy deposition", energyBinN, 0, 100.);
    hist_energy_tot_ = data_monitor_.add_hist<TH1D>("hEtot", "Total Energy", totenergyBinN, 0, total_energy_upper);
    hdeltaEE = data_monitor_.add_hist<TH2D>("hdeltaEE",
                                            "Energy in Foremost Plane vs Etot",
                                            energyBinN,
                                            0,
                                            energy_upper_foremost_vs_e_tot_a,
                                            energyBinN,
                                            0,
                                            energy_upper_foremost_vs_e_tot_b);
    hist_pos_vs_energy_ = data_monitor_.add_hist<TH2D>(
        "hPosVSEnergy", "Position vs Energy deposition", zDepBinN, z_dep_lower, z_dep_upper, totenergyBinN, 0, 1000.);
    hist_beta_ = data_monitor_.add_hist<TH1D>("hBeta", "Velocity", velocityBinN, 0., 1.);
    hist_energy_ = data_monitor_.add_hist<TH1D>("hE", "Hit Energy", energyBinN, 0., 100.);
    hist_x_ = data_monitor_.add_hist<TH1D>("hX", "Hit X", posXYBinN, pos_xy_lower, pos_xy_upper);
    hist_y_ = data_monitor_.add_hist<TH1D>("hY", "Hit Y", posXYBinN, pos_xy_lower, pos_xy_upper);
    hT = data_monitor_.add_hist<TH1D>("hT", "Hit Delta T", timeBinN, time_para, time_para);
    hTNeigh = data_monitor_.add_hist<TH1D>("hTNeigh", "Hit Neigh Delta T", timeBinN, time_para, time_para);

    return kSUCCESS;
}

void R3BNeulandHitMon::Exec(Option_t* /*option*/)
{
    const auto hits = neuland_hits_.get();

    // checking paddle multihits
    std::map<Int_t, Int_t> paddlenum;
    for (const auto& hit : hits)
    {
        auto result = paddlenum.insert(std::pair<Int_t, Int_t>(hit.GetPaddle(), 1));
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

    if (is_3d_track_enabled_)
    {
        hist_3_->Reset("ICES");
        for (const auto& hit : hits)
        {
            hist_3_->Fill(hit.GetPosition().Z(), hit.GetPosition().X(), hit.GetPosition().Y(), hit.GetE());
        }
    }

    hist_mult_->Fill(static_cast<double>(hits.size()));
    for (const auto& hit : hits)
    {
        hist_pos_vs_energy_->Fill(hit.GetPosition().Z(), hit.GetE());
        hist_time_->Fill(hit.GetT());
        hist_time_adj_->Fill(distance_to_target_ / hit.GetPosition().Mag() * hit.GetT());
        hist_beta_->Fill(hit.GetBeta());
        hist_energy_->Fill(hit.GetE());
        hist_x_->Fill(hit.GetPosition().X());
        hist_y_->Fill(hit.GetPosition().Y());
    }

    for (auto it1 = hits.begin(); it1 != hits.end(); it1++)
    {
        for (auto it2 = it1 + 1; it2 != hits.end(); it2++)
        {
            const auto pos_delta_cutoff = 7.5;
            if (std::abs((*it1).GetPosition().X() - (*it2).GetPosition().X()) < pos_delta_cutoff &&
                std::abs((*it1).GetPosition().Y() - (*it2).GetPosition().Y()) < pos_delta_cutoff &&
                std::abs((*it1).GetPosition().Z() - (*it2).GetPosition().Z()) < pos_delta_cutoff)
            {
                hTNeigh->Fill((*it1).GetT() - (*it2).GetT());
            }
            hT->Fill((*it1).GetT() - (*it2).GetT());
        }
    }

    auto maxDepthHit = std::max_element(hits.begin(),
                                        hits.end(),
                                        [](const R3BNeulandHit& one, const R3BNeulandHit& another)
                                        { return one.GetPosition().Z() < another.GetPosition().Z(); });
    if (maxDepthHit != hits.end())
    {
        hist_depth_->Fill((*maxDepthHit).GetPosition().Z());
        hist_sternmost_energy_->Fill((*maxDepthHit).GetE());
        hist_depth_vs_sternmost_energy_->Fill((*maxDepthHit).GetPosition().Z(), (*maxDepthHit).GetE());
    }

    auto minDepthHit = std::min_element(hits.begin(),
                                        hits.end(),
                                        [](const R3BNeulandHit& one, const R3BNeulandHit& another)
                                        { return one.GetPosition().Z() < another.GetPosition().Z(); });
    auto Etot =
        std::accumulate(hits.begin(), hits.end(), 0., [](double init, const auto& hit) { return init + hit.GetE(); });

    if (minDepthHit != hits.end())
    {
        hist_foremost_energy_->Fill((*minDepthHit).GetE());
        hist_depth_vs_foremost_energy_->Fill((*maxDepthHit).GetPosition().Z(), (*minDepthHit).GetE());
        hdeltaEE->Fill(Etot, (*minDepthHit).GetE());
    }

    hist_energy_tot_->Fill(Etot);
    if (maxDepthHit != hits.end())
    {
        hist_depth_vs_energy_tot_->Fill((*maxDepthHit).GetPosition().Z(), Etot);
    }
}

ClassImp(R3BNeulandHitMon) // NOLINT
