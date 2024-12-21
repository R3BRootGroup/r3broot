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

#include "R3BNeulandClusterMon.h"
#include "ElasticScattering.h"
#include "FairLogger.h"
#include "FairRootManager.h"
#include "R3BNeulandCluster.h"
#include "TClonesArray.h"
#include "TDirectory.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TH3D.h"
#include <TFile.h>
#include <algorithm>
#include <iostream>
#include <numeric>
#include <utility>

constexpr double rad2deg = 180. / 3.141592653589793238463;

inline Double_t GetTheta(const R3BNeulandCluster* cluster)
{
    const auto direction = cluster->GetLastHit().GetPosition() - cluster->GetFirstHit().GetPosition();
    const auto x = std::acos(direction.Y() / direction.Mag()) * rad2deg;
    // Not sure, but Kondos Theta is -90:90
    return x - 90.;
}

R3BNeulandClusterMon::R3BNeulandClusterMon(TString input, TString output, const Option_t* option)
    : FairTask("R3B NeuLAND NeulandCluster Monitor")
    , fNeulandClusters(input)
    , fOutput(std::move(output))
    , fBeta(0.793) // 600 Mev?
{
    LOG(info) << "Using R3B NeuLAND NeulandCluster Monitor";

    TString opt = option;
    opt.ToUpper();

    if (opt.Contains("3DTRACK"))
    {
        fIs3DTrackEnabled = true;
        LOG(info) << "... with 3D track visualization";
    }
    else
    {
        fIs3DTrackEnabled = false;
    }
}

InitStatus R3BNeulandClusterMon::Init()
{
    fNeulandClusters.Init();

    FairRootManager* ioman = FairRootManager::Instance();
    if (!ioman)
    {
        LOG(fatal) << "R3BNeulandClusterMon::Init: No FairRootManager";
        return kFATAL;
    }

    if (fIs3DTrackEnabled)
    {
        // XYZ -> ZXY (side view)
        fh3 = new TH3D("hClusters", "hClusters", 60, 1400, 1700, 50, -125, 125, 50, -125, 125);
        fh3->SetTitle("NeuLAND Clusters");
        fh3->GetXaxis()->SetTitle("Z");
        fh3->GetYaxis()->SetTitle("X");
        fh3->GetZaxis()->SetTitle("Y");

        ioman->Register(fOutput, "Cluster TH3Ds in NeuLAND", fh3, kTRUE);
    }

    TH1::AddDirectory(kFALSE);

    fhClusters = new TH1D("nClusters", "Number of clusters in one event", 100, 0, 100);
    fhClusterSize = new TH1D("ClusterSize", "Number of Digis for each Cluster", 100, 0, 100);
    fhClusterEnergy = new TH1D("ClusterEnergy", "Energy for each Cluster", 2000, 0., 2000.);
    fhClusterEnergyVSSize =
        new TH2D("ClusterEnergyVSSize", "Energy for each Cluster vs Cluster Size", 2000, 0., 2000., 100, 0., 100.);
    fhClusterTime = new TH1D("ClusterTime", "Time for each Cluster", 5000, 0., 500.);
    fhClusterEToF = new TH1D("ClusterEToF", "Cluster EToF", 2000, 0, 2000);
    fhClusterRValue = new TH1D("ClusterRValue", "Log Cluster R Value", 20000, -10, 1);

    fhClusterSizeVSEToF =
        new TH2D("ClusterSizeVSEToF", "Number of Digis for each Cluster vs EToF", 100, 0, 100, 2000, 0, 2000);
    fhClusterEnergyVSEToF =
        new TH2D("ClusterEnergyVSEToF", "Energy for each Cluster vs EToF", 2000, 0., 2000., 2000, 0, 2000);
    fhClusterEnergyVSSizeVSEToF = new TH3D("ClusterEnergyVSSizeVSEToF",
                                           "Energy for each Cluster vs Cluster Size vs EToF",
                                           1000,
                                           0.,
                                           1000.,
                                           100,
                                           0.,
                                           100.,
                                           1200,
                                           0,
                                           1200.);

    fhENFromScatterVSEToF = new TH2D("ENFromScatterVSEToF", "ENFromScatterVSEToF", 1000, 0, 1000, 1000, 0, 1000);

    fhClusterNumberVSEnergy =
        new TH2D("ClusterNumberEtot", "Number of Clusters vs. Total Energy", 200, 0, 2000, 50, 0, 50);
    fhClusterNumberVSEnergy->GetXaxis()->SetTitle("Total energy [MeV]");
    fhClusterNumberVSEnergy->GetYaxis()->SetTitle("Number of Clusters");

    fhClusterEToFVSEnergy =
        new TH2D("ClusterEToFVSEnergy", "Cluster E_{ToF} vs. Cluster Energy", 100, 0, 1000, 100, 0, 1000);
    fhClusterEToFVSEnergy->GetXaxis()->SetTitle("Cluster E_{ToF} [MeV]");
    fhClusterEToFVSEnergy->GetYaxis()->SetTitle("Cluster E [MeV]");

    fhClusterEToFVSTime = new TH2D("ClusterEToFVSTime", "Cluster E_{ToF} vs. Cluster Time", 100, 0, 1000, 500, 0, 500);
    fhClusterEToFVSTime->GetXaxis()->SetTitle("Cluster E_{ToF} [MeV]");
    fhClusterEToFVSTime->GetYaxis()->SetTitle("Cluster t [ns]");

    fhClusterEVSTime = new TH2D("ClusterEVSTime", "Cluster E vs. Cluster Time", 250, 0, 250, 500, 0, 250);
    fhClusterEVSTime->GetXaxis()->SetTitle("Cluster E [MeV]");
    fhClusterEVSTime->GetYaxis()->SetTitle("Cluster t [ns]");

    fhClusterLastMinusFirstDigiMagVSEnergy = new TH2D(
        "ClusterLastMinusFirstDigiMagVSEnergy", "ClusterLastMinusFirstDigiMagVSEnergy", 101, 0, 100, 60, 0, 600);

    fhClusterForemostMinusCentroidVSEnergy =
        new TH2D("ClusterForemostMinusCentroidVSEnergy",
                 "Distance between Foremost Digi and Energy Centroid vs Cluster Energy",
                 101,
                 0,
                 100,
                 60,
                 0,
                 600);
    fhClusterForemostMinusCentroidVSEnergy->GetXaxis()->SetTitle("|#vec{r_{z_{min}}}-#vec{r_{EC}}| [cm]");
    fhClusterForemostMinusCentroidVSEnergy->GetYaxis()->SetTitle("Cluster Energy [MeV]");

    fhClusterForemostMinusMaxEnergyDigiPosVSEnergy =
        new TH2D("ClusterForemostMinusMaxEnergyDigiPosVSEnergy",
                 "Distance between Foremost Digi and Max Energy Digi vs Cluster Energy",
                 101,
                 0,
                 100,
                 60,
                 0,
                 600);
    fhClusterForemostMinusMaxEnergyDigiPosVSEnergy->GetXaxis()->SetTitle("|#vec{r_{z_{min}}}-#vec{r_{E_{max}}}| [cm]");
    fhClusterForemostMinusMaxEnergyDigiPosVSEnergy->GetYaxis()->SetTitle("Cluster Energy [MeV]");

    fhClusterCentroidMinusFirstDigiPosVSEnergy =
        new TH2D("ClusterCentroidMinusFirstDigiPosVSEnergy",
                 "Distance between Cluster Energy Centroid and First Hit vs Cluster Energy",
                 101,
                 0,
                 100,
                 60,
                 0,
                 600);
    fhClusterCentroidMinusFirstDigiPosVSEnergy->GetXaxis()->SetTitle("|#vec{r_{EC}}-#vec{r_{t_{min}}}| [cm]");
    fhClusterCentroidMinusFirstDigiPosVSEnergy->GetYaxis()->SetTitle("Cluster Energy [MeV]");

    fhClusterMaxEnergyDigiMinusFirstDigiPosVSEnergy =
        new TH2D("ClusterMaxEnergyDigiMinusFirstDigiPosVSEnergy",
                 "Distance between Max Energy Digi and First Digi vs Cluster Energy",
                 101,
                 0,
                 100,
                 60,
                 0,
                 600);
    fhClusterMaxEnergyDigiMinusFirstDigiPosVSEnergy->GetXaxis()->SetTitle("|#vec{r_{E_{max}}}-#vec{r_{t_{min}}}| [cm]");
    fhClusterMaxEnergyDigiMinusFirstDigiPosVSEnergy->GetYaxis()->SetTitle("Cluster Energy [MeV]");

    fhClusterMaxEnergyDigiMinusCentroidVSEnergy =
        new TH2D("ClusterMaxEnergyDigiMinusCentroidVSEnergy",
                 "Distance between Max Energy Digi and Energy Centroid vs Cluster Energy",
                 101,
                 0,
                 100,
                 60,
                 0,
                 600);
    fhClusterMaxEnergyDigiMinusCentroidVSEnergy->GetXaxis()->SetTitle("|#vec{r_{E_{max}}}-#vec{r_{EC}}| [cm]");
    fhClusterMaxEnergyDigiMinusCentroidVSEnergy->GetYaxis()->SetTitle("Cluster Energy [MeV]");

    fhClusterEnergyMomentVSEnergy =
        new TH2D("ClusterEnergyMomentVSEnergy", "Cluster Energy Moment VS Energy", 100, 0, 100, 60, 0, 600);
    fhClusterEnergyMomentVSEnergy->GetXaxis()->SetTitle("|#vec{EM}| [cm]");
    fhClusterEnergyMomentVSEnergy->GetYaxis()->SetTitle("Cluster Energy [MeV]");

    fhClusterEnergyMoment = new TH1D("ClusterEnergyMoment", "Cluster Energy Moment", 100, 0, 100);
    fhClusterMaxEnergyDigiMinusFirstDigiMag =
        new TH1D("ClusterMaxEnergyDigiMinusFirstDigiMag", "ClusterMaxEnergyDigiMinusFirstDigiMag", 100, 0, 1000);

    fhClusterEnergyMomentVSClusterSize =
        new TH2D("ClusterEnergyMomentVSClusterSize", "Cluster Energy Moment VS Cluster Size", 100, 0, 100, 30, 0, 30);

    fhEToFVSEelastic = new TH2D("EToFVSEelastic", "EToFVSEelastic", 50, 0, 1000, 50, 10, 1000);
    fhScatteredNEnergyVSAngle = new TH2D("ScatteredNEnergyVSAngle", "ScatteredNEnergyVSAngle", 50, 0, 1000, 50, 0, 3.2);
    fhScatteredNEnergyVSEdep = new TH2D("ScatteredNEnergyVSEdep", "ScatteredNEnergyVSEdep", 100, 0, 1000, 100, 0, 1000);

    fhScatterAngleVSRecoilAngle =
        new TH2D("ScatterAngleVSRecoilAngle", "ScatterAngleVSRecoilAngle", 50, 0, 3.2, 50, 0, 3.2);

    fhSumAngleVSRatioErecoEtof = new TH2D("SumAngleVSRatioErecoEtof", "SumAngleVSRatioErecoEtof", 50, 0, 3.2, 50, 0, 2);
    fhClusterEnergyVSScatteredRecoilAngle =
        new TH2D("ClusterEnergyVSScatteredRecoilAngle", "ClusterEnergyVSScatteredRecoilAngle", 50, 0, 1000, 50, 0, 3.2);
    fhClusterEnergyVSScatteredNeutronAngle = new TH2D(
        "ClusterEnergyVSScatteredNeutronAngle", "ClusterEnergyVSScatteredNeutronAngle", 50, 0, 1000, 50, 0, 3.2);

    fhElasticTargetMass = new TH2D("ElasticTargetMass", "ElasticTargetMass", 200, -50000, 50000, 100, 0, 1000);

    /*fhtheh = new TH2D("theh", "", 100, 0, 1000, 100, 0, 2000);
    ediff = new TH2D("ediff", "Reconstructed Proton Energy VS Measured Cluster Energy", 100, 0, 300, 100, 0, 300);
    ediff->GetXaxis()->SetTitle("E_{p'}");
    ediff->GetYaxis()->SetTitle("E_{cluster}");

    etheta = new TH2D("etheta", "", 100, 0, 1000, 100, -1, 1);

    fhProtonTrack*/

    fhZ = new TH1D("Z", "Z", 5000, 0., 5000.);
    fhZVSEToF = new TH2D("ZVSEToF", "Z vs EToF", 4000, 1000, 5000, 2000, 0, 2000);
    fhDistFromCenter = new TH1D("DistFromCenter", "DistFromCenter", 500, 0., 500.);
    fhDistFromCenterVSEToF = new TH2D("DistFromCenterVSEToF", "DistFromCenter vs EToF", 500, 0, 500, 2000, 0, 2000);

    fhDeltaT = new TH1D("fhDeltaT", "T_{Last Digi} - T_{First Digi}", 1000, 0, 100);
    fhForemostMinusFirstDigiTime =
        new TH1D("fhForemostMinusFirstDigiTime", "T_{Foremost} - T_{First Digi}", 1000, -10, 10);

    fhThetaEDigi = new TH2D("fhThetaEDigi", "fhThetaEDigi", 1000, -500, 500, 400, 0, 400);
    fhThetaEDigiCosTheta = new TH2D("fhThetaEDigiCosTheta", "fhThetaEDigiCosTheta", 1000, -500, 500, 400, 0, 400);

    hT = new TH1D("hT", "Cluster Digi Delta T", 3000, 0., 15.);
    hTNeigh = new TH1D("hTNeigh", "Cluster Digi Neigh Delta T", 3000, 0., 15.);

    return kSUCCESS;
}

void R3BNeulandClusterMon::Exec(Option_t*)
{
    auto clusters = fNeulandClusters.Retrieve();
    clusters.erase(
        std::remove_if(
            clusters.begin(), clusters.end(), [&](R3BNeulandCluster* c) { return !(fClusterFilters.IsValid(c)); }),
        clusters.end());

    const auto nClusters = clusters.size();

    if (fIs3DTrackEnabled)
    {
        fh3->Reset("ICES");
        for (const auto& cluster : clusters)
        {
            const auto start = cluster->GetFirstHit().GetPosition();
            // XYZ -> ZXY (side view)
            fh3->Fill(start.Z(), start.X(), start.Y(), cluster->GetE());
        }
    }

    const Double_t etot = std::accumulate(
        clusters.begin(), clusters.end(), 0., [](Double_t sum, const R3BNeulandCluster* c) { return sum + c->GetE(); });

    fhClusterNumberVSEnergy->Fill(etot, nClusters);

    fhClusters->Fill(nClusters);
    for (const auto cluster : clusters)
    {
        fhClusterTime->Fill(cluster->GetT());
        fhClusterSize->Fill(cluster->GetSize());
        fhClusterEnergy->Fill(cluster->GetE());
        fhClusterRValue->Fill(std::log10(cluster->GetRCluster(fBeta)));
        fhClusterEnergyVSSize->Fill(cluster->GetE(), cluster->GetSize());
        fhClusterEToF->Fill(cluster->GetFirstHit().GetEToF());
        fhClusterEToFVSEnergy->Fill(cluster->GetFirstHit().GetEToF(), cluster->GetE());
        fhClusterEToFVSTime->Fill(cluster->GetFirstHit().GetEToF(), cluster->GetT());
        fhClusterEVSTime->Fill(cluster->GetFirstHit().GetE(), cluster->GetT());

        fhClusterEnergyVSEToF->Fill(cluster->GetE(), cluster->GetFirstHit().GetEToF());
        fhClusterSizeVSEToF->Fill(cluster->GetSize(), cluster->GetFirstHit().GetEToF());
        fhClusterEnergyVSSizeVSEToF->Fill(cluster->GetE(), cluster->GetSize(), cluster->GetFirstHit().GetEToF());
        if (cluster->GetSize() > 2)
        {
            fhClusterForemostMinusCentroidVSEnergy->Fill(
                (cluster->GetForemostHit().GetPosition() - cluster->GetEnergyCentroid()).Mag(), cluster->GetE());

            fhClusterForemostMinusMaxEnergyDigiPosVSEnergy->Fill(
                (cluster->GetForemostHit().GetPosition() - cluster->GetMaxEnergyHit().GetPosition()).Mag(),
                cluster->GetE());

            fhClusterCentroidMinusFirstDigiPosVSEnergy->Fill(
                (cluster->GetEnergyCentroid() - cluster->GetFirstHit().GetPosition()).Mag(), cluster->GetE());

            fhClusterMaxEnergyDigiMinusFirstDigiPosVSEnergy->Fill(
                (cluster->GetMaxEnergyHit().GetPosition() - cluster->GetFirstHit().GetPosition()).Mag(),
                cluster->GetE());
            fhClusterMaxEnergyDigiMinusCentroidVSEnergy->Fill(
                (cluster->GetMaxEnergyHit().GetPosition() - cluster->GetEnergyCentroid()).Mag(), cluster->GetE());
            fhClusterEnergyMomentVSEnergy->Fill(cluster->GetEnergyMoment(), cluster->GetE());
            fhClusterEnergyMomentVSClusterSize->Fill(cluster->GetEnergyMoment(), cluster->GetSize());

            fhClusterLastMinusFirstDigiMagVSEnergy->Fill(
                (cluster->GetLastHit().GetPosition() - cluster->GetFirstHit().GetPosition()).Mag(), cluster->GetE());

            fhClusterEnergyMoment->Fill(cluster->GetEnergyMoment());
            fhClusterMaxEnergyDigiMinusFirstDigiMag->Fill(
                (cluster->GetMaxEnergyHit().GetPosition() - cluster->GetFirstHit().GetPosition()).Mag());
        }

        fhZ->Fill(cluster->GetFirstHit().GetPosition().Z());
        fhZVSEToF->Fill(cluster->GetFirstHit().GetPosition().Z(), cluster->GetEToF());
        fhDistFromCenter->Fill(std::sqrt(std::pow(cluster->GetFirstHit().GetPosition().X(), 2) +
                                         std::pow(cluster->GetFirstHit().GetPosition().Y(), 2)));
        fhDistFromCenterVSEToF->Fill(std::sqrt(std::pow(cluster->GetFirstHit().GetPosition().X(), 2) +
                                               std::pow(cluster->GetFirstHit().GetPosition().Y(), 2)),
                                     cluster->GetEToF());
        fhDeltaT->Fill(cluster->GetLastHit().GetT() - cluster->GetFirstHit().GetT());

        fhForemostMinusFirstDigiTime->Fill(cluster->GetForemostHit().GetT() - cluster->GetFirstHit().GetT());

        if (cluster->GetSize() > 4)
        {
            const auto theta = GetTheta(cluster);
            for (const auto& digi : cluster->GetHits())
            {
                fhThetaEDigi->Fill(theta, digi.GetE());
                fhThetaEDigiCosTheta->Fill(theta, digi.GetE() * std::cos(theta / rad2deg));
            }
        }
    }

    std::sort(clusters.begin(),
              clusters.end(),
              [](const R3BNeulandCluster* a, const R3BNeulandCluster* b) { return a->GetT() < b->GetT(); });

    for (auto cluster : clusters)
    {
        if (cluster->GetSize() >= 3)
        {
            fhENFromScatterVSEToF->Fill(Neuland::NeutronEnergyFromElasticProtonScattering(cluster), cluster->GetEToF());

            fhClusterEnergyVSScatteredRecoilAngle->Fill(cluster->GetE(),
                                                        std::acos(Neuland::RecoilScatteringAngle(cluster)));
        }
    }

    for (auto ita = clusters.cbegin(); ita != clusters.cend(); ita++)
    {
        for (auto itb = ita + 1; itb != clusters.cend(); itb++)
        {

            fhScatteredNEnergyVSAngle->Fill(Neuland::ScatteredNeutronEnergy(*ita, *itb),
                                            std::acos(Neuland::ScatteredNeutronAngle(*ita, *itb)));

            fhScatteredNEnergyVSEdep->Fill(Neuland::ScatteredNeutronEnergy(*ita, *itb), (*ita)->GetE());

            const Double_t EelasticHeavy = Neuland::NeutronEnergyFromElasticScattering(*ita, *itb, 11000);
            fhEToFVSEelastic->Fill((*ita)->GetFirstHit().GetEToF(), EelasticHeavy);
            // const Double_t EelasticProton = Neuland::NeutronEnergyFromElasticScattering(*ita, *itb, 1000);
            // fhEToFVSEelastic->Fill((*ita)->GetFirstHit().GetEToF(), EelasticProton);

            if (Neuland::ScatteredNeutronEnergy(*ita, *itb) > 10.)
            {
                fhElasticTargetMass->Fill(Neuland::ElasticScatteringTargetMass(*ita, *itb), (*ita)->GetE());
                fhClusterEnergyVSScatteredNeutronAngle->Fill((*ita)->GetE(),
                                                             std::acos(Neuland::ScatteredNeutronAngle(*ita, *itb)));

                if ((*ita)->GetSize() >= 3)
                {
                    fhScatterAngleVSRecoilAngle->Fill(std::acos(Neuland::ScatteredNeutronAngle(*ita, *itb)),
                                                      std::acos(Neuland::RecoilScatteringAngle(*ita)));
                    fhSumAngleVSRatioErecoEtof->Fill(
                        std::acos(Neuland::ScatteredNeutronAngle(*ita, *itb)) +
                            std::acos(Neuland::RecoilScatteringAngle(*ita)),
                        Neuland::NeutronEnergyFromElasticProtonScattering(*ita) / (*ita)->GetFirstHit().GetEToF());
                }
            }
        }
    }

    for (const auto& cluster : clusters)
    {
        const auto& digis = cluster->GetHits();
        for (auto it1 = digis.begin(); it1 != digis.end(); it1++)
        {
            for (auto it2 = it1 + 1; it2 != digis.end(); it2++)
            {
                if (std::abs(it1->GetPosition().X() - it2->GetPosition().X()) < 7.5 &&
                    std::abs(it1->GetPosition().Y() - it2->GetPosition().Y()) < 7.5 &&
                    std::abs(it1->GetPosition().Z() - it2->GetPosition().Z()) < 7.5)
                {
                    hTNeigh->Fill(std::abs(it1->GetT() - it2->GetT()));
                }
                hT->Fill(std::abs(it1->GetT() - it2->GetT()));
            }
        }
    }
}

void R3BNeulandClusterMon::Finish()
{
    TDirectory* tmp = gDirectory;
    FairRootManager::Instance()->GetOutFile()->cd();

    gDirectory->mkdir(fOutput);
    gDirectory->cd(fOutput);

    fhClusters->Write();
    fhClusterSize->Write();
    fhClusterEnergyVSSize->Write();
    fhClusterEnergy->Write();
    fhClusterRValue->Write();
    fhClusterTime->Write();
    fhClusterNumberVSEnergy->Write();
    fhClusterEToF->Write();
    fhClusterEToFVSEnergy->Write();
    fhClusterEToFVSTime->Write();
    fhElasticTargetMass->Write();

    fhClusterSizeVSEToF->Write();
    fhClusterEnergyVSEToF->Write();
    // fhClusterEnergyVSSizeVSEToF->Write();

    fhClusterForemostMinusCentroidVSEnergy->Write();
    fhClusterForemostMinusMaxEnergyDigiPosVSEnergy->Write();
    fhClusterCentroidMinusFirstDigiPosVSEnergy->Write();
    fhClusterMaxEnergyDigiMinusFirstDigiPosVSEnergy->Write();
    fhClusterMaxEnergyDigiMinusCentroidVSEnergy->Write();
    fhClusterEnergyMomentVSEnergy->Write();
    fhClusterEnergyMomentVSClusterSize->Write();
    fhClusterEVSTime->Write();
    fhClusterEnergyMoment->Write();
    fhClusterMaxEnergyDigiMinusFirstDigiMag->Write();
    fhClusterLastMinusFirstDigiMagVSEnergy->Write();

    fhENFromScatterVSEToF->Write();

    fhEToFVSEelastic->Write();
    fhScatteredNEnergyVSAngle->Write();
    fhScatteredNEnergyVSEdep->Write();
    fhScatterAngleVSRecoilAngle->Write();
    fhSumAngleVSRatioErecoEtof->Write();
    fhClusterEnergyVSScatteredRecoilAngle->Write();
    fhClusterEnergyVSScatteredNeutronAngle->Write();

    fhZ->Write();
    fhZVSEToF->Write();
    fhDistFromCenter->Write();
    fhDistFromCenterVSEToF->Write();
    fhDeltaT->Write();
    fhForemostMinusFirstDigiTime->Write();

    fhThetaEDigi->Write();
    fhThetaEDigiCosTheta->Write();

    hT->Write();
    hTNeigh->Write();

    gDirectory = tmp;
}

ClassImp(R3BNeulandClusterMon)
