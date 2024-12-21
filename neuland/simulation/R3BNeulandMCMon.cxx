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

#include "R3BNeulandMCMon.h"
#include "FairLogger.h"
#include "FairRootManager.h"
#include "FairRun.h"
#include "Math/Vector4D.h"
#include "TDirectory.h"
#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TH3D.h"
#include <algorithm>
#include <iostream>
#include <numeric>
#include <string>

inline Bool_t IsPrimaryNeutron(const R3BMCTrack* mcTrack)
{
    return (mcTrack->GetPdgCode() == 2112 && mcTrack->GetMotherId() == -1);
}

inline Bool_t IsMotherPrimaryNeutron(const R3BMCTrack* mcTrack, const std::vector<R3BMCTrack*>& tracks)
{
    if (mcTrack->GetMotherId() < 0 || size_t(mcTrack->GetMotherId()) >= tracks.size())
    {
        return false;
    }
    return IsPrimaryNeutron(tracks.at(mcTrack->GetMotherId()));
}

inline Double_t GetKineticEnergy(const R3BMCTrack* mcTrack)
{
    return (mcTrack->GetEnergy() - mcTrack->GetMass()) * 1000.;
}

inline Double_t GetTheta(const R3BMCTrack* mcTrack) { return std::acos(mcTrack->GetPz() / mcTrack->GetP()); }

static std::map<int, std::string> lookup_table = {
    { -211, "pion" }, { 22, "gamma" },        { 111, "pion" },        { 211, "pion" },         { 2112, "n" },
    { 2212, "p" },    { 1000010020, "d, t" }, { 1000010030, "d, t" }, { 1000020030, "alpha" }, { 1000020040, "alpha" },
};

std::string lookup(int pid)
{
    const auto l = lookup_table.find(pid);
    if (l == lookup_table.end())
    {
        return std::string("heavy");
    }
    else
    {
        return l->second;
    }
}

template <typename T>
void writeout(std::map<T, TH1D*>& map, const std::string& what, const int nEvents = 0)
{
    if (!map.empty())
    {
        std::ostream* out = &std::cout;
        std::ofstream of;
        if (FairRun::Instance() != nullptr)
        {
            const std::string f =
                std::string{ FairRun::Instance()->GetUserOutputFileName().Data() } + std::string(".") + what + ".dat";
            std::cout << "Writing to file " << f << std::endl;
            of.open(f);
            out = &of;
        }
        gDirectory->mkdir(what.c_str());
        gDirectory->cd(what.c_str());
        for (const auto& kv : map)
        {
            auto h = kv.second;
            const auto c = h->GetEntries();
            if (nEvents > 0 && nEvents > c)
            {
                h->Fill(0., nEvents - c);
            }
            *out << kv.first << "\t" << c << "\t" << h->GetMean() << "\t" << h->GetMeanError() << std::endl;
            kv.second->Write();
        }
        gDirectory->cd("..");
    }
}

R3BNeulandMCMon::R3BNeulandMCMon(const Option_t* option)
    : FairTask("R3B NeuLAND Neuland Monte Carlo Monitor")
    , fIs3DTrackEnabled(false)
    , fIsFullSimAnaEnabled(false)
    , fPrimaryNeutronInteractionPoints("NeulandPrimaryPoints")
    , fMCTracks("MCTrack")
    , fNeulandPoints("NeulandPoints")
    , nEvents(0)
{
    LOG(info) << "Using R3B NeuLAND Neuland Monte Carlo Monitor";

    TString opt = option;
    opt.ToUpper();
    if (opt.Contains("3DTRACK"))
    {
        fIs3DTrackEnabled = true;
        LOG(info) << "... with 3D track visualization";
    }
    if (opt.Contains("FULLSIMANA"))
    {
        fIsFullSimAnaEnabled = true;
        LOG(info) << "... with full simulation neutron reaction product analysis";
    }
}

InitStatus R3BNeulandMCMon::Init()
{
    fPrimaryNeutronInteractionPoints.Init();
    fMCTracks.Init();
    fNeulandPoints.Init();

    fhNPNIPsEToFVSTime = new TH2D("NPNIPEToFVSTime", "NPNIP E_{ToF} vs. NPNIP Time", 100, 0, 1000, 500, 0, 500);
    fhNPNIPsEToFVSTime->GetXaxis()->SetTitle("NPNIP E_{ToF} [MeV]");
    fhNPNIPsEToFVSTime->GetYaxis()->SetTitle("NPNIP t [ns]");

    fhPDG = new TH1D("hPDG",
                     "Number of particles by PDG code created by primary neutron interaction resulting in a point",
                     8000,
                     -4000,
                     4000);
    fhEPrimarys = new TH1D("hE_primarys", "Energy of primary particles", 10000, 0, 10000);
    fhEPrimarys2 = new TH1D("hE_primarys2", "Tot. Energy of primary particles", 100000, 0, 100000);
    fhEPrimaryNeutrons = new TH1D("hE_primary_neutrons", "Energy of primary Neutrons", 10000, 0, 10000);
    fhErelMC = new TH1D("hErelMC", "Erel", 5000, 0, 5000);
    fhEtotPrim = new TH1D("hE_tot_prim",
                          "Total Light Yield of non-neutron LandPoints created by primary neutron interaction(s)",
                          6000,
                          0,
                          6000);
    fhEtot = new TH1D("hE_tot", "Total Light Yield of non-neutron LandPoints", 6000, 0, 6000);
    fhESecondaryNeutrons = new TH1D(
        "hE_secondary_neutrons", "Energy of neutron tracks created by primary neutron interaction", 6000, 0, 6000);
    fhMotherIDs = new TH1D("hmotherIDs", "MotherIDs", 6010, -10, 6000);
    fhPrimaryDaughterIDs = new TH1D("hprimary_daughter_IDs", "IDs of tracks with a primary mother", 6001, -1, 6000);
    fhMCToF = new TH1D("fhMCToF", "Energy of primary Neutron - ToF Energy from PNIPS", 2001, -1000, 1000);

    fhNPNIPSrvsz = new TH2D("fhNPNIPSrvsz", "NPNIPS R = #sqrt(x**2+y**2) vs z", 600, 1000, 4000, 40, 0, 200);
    fhNPNIPSrvsz->GetXaxis()->SetTitle("NPNIP Z [cm]");
    fhNPNIPSrvsz->GetYaxis()->SetTitle("NPNIP R [cm]");

    fhNPNIPSxy = new TH2D("fhNPNIPSxy", "NPNIPS xy", 101, -200, 200, 101, -200, 200);
    fhNPNIPSxy->GetXaxis()->SetTitle("NPNIP X [cm]");
    fhNPNIPSxy->GetYaxis()->SetTitle("NPNIP Y [cm]");

    fhnNPNIPs = new TH1D("fhnNPNIPs", "Number of NPNIPs per event", 10, -0.5, 9.5);

    fhnProducts = new TH1D("nProducts", "Number of products created by primary neutron interaction", 20, -0.5, 19.5);
    fhnProductsCharged =
        new TH1D("nProductsCharged",
                 "Number of products created by primary neutron interaction w/o gammas, neutrons, pion0",
                 20,
                 -0.5,
                 19.5);
    fhSumProductEnergy = new TH1D("SumProductEnergy", "Sum product energy", 1000, 0, 1000);
    fhnSecondaryNeutrons = new TH1D("NumberSecondaryNeutrons", "Number of Secondary Neutrons", 10, -0.5, 9.5);
    fhnSecondaryProtons = new TH1D("NumberSecondaryProtons", "Number of Secondary Protons", 10, -0.5, 9.5);

    // bins = 10**numpy.arange(-2,2.42,0.03)
    std::array<double, 148> bins = {
        1.00000000e-02, 1.07151931e-02, 1.14815362e-02, 1.23026877e-02, 1.31825674e-02, 1.41253754e-02, 1.51356125e-02,
        1.62181010e-02, 1.73780083e-02, 1.86208714e-02, 1.99526231e-02, 2.13796209e-02, 2.29086765e-02, 2.45470892e-02,
        2.63026799e-02, 2.81838293e-02, 3.01995172e-02, 3.23593657e-02, 3.46736850e-02, 3.71535229e-02, 3.98107171e-02,
        4.26579519e-02, 4.57088190e-02, 4.89778819e-02, 5.24807460e-02, 5.62341325e-02, 6.02559586e-02, 6.45654229e-02,
        6.91830971e-02, 7.41310241e-02, 7.94328235e-02, 8.51138038e-02, 9.12010839e-02, 9.77237221e-02, 1.04712855e-01,
        1.12201845e-01, 1.20226443e-01, 1.28824955e-01, 1.38038426e-01, 1.47910839e-01, 1.58489319e-01, 1.69824365e-01,
        1.81970086e-01, 1.94984460e-01, 2.08929613e-01, 2.23872114e-01, 2.39883292e-01, 2.57039578e-01, 2.75422870e-01,
        2.95120923e-01, 3.16227766e-01, 3.38844156e-01, 3.63078055e-01, 3.89045145e-01, 4.16869383e-01, 4.46683592e-01,
        4.78630092e-01, 5.12861384e-01, 5.49540874e-01, 5.88843655e-01, 6.30957344e-01, 6.76082975e-01, 7.24435960e-01,
        7.76247117e-01, 8.31763771e-01, 8.91250938e-01, 9.54992586e-01, 1.02329299e+00, 1.09647820e+00, 1.17489755e+00,
        1.25892541e+00, 1.34896288e+00, 1.44543977e+00, 1.54881662e+00, 1.65958691e+00, 1.77827941e+00, 1.90546072e+00,
        2.04173794e+00, 2.18776162e+00, 2.34422882e+00, 2.51188643e+00, 2.69153480e+00, 2.88403150e+00, 3.09029543e+00,
        3.31131121e+00, 3.54813389e+00, 3.80189396e+00, 4.07380278e+00, 4.36515832e+00, 4.67735141e+00, 5.01187234e+00,
        5.37031796e+00, 5.75439937e+00, 6.16595002e+00, 6.60693448e+00, 7.07945784e+00, 7.58577575e+00, 8.12830516e+00,
        8.70963590e+00, 9.33254301e+00, 1.00000000e+01, 1.07151931e+01, 1.14815362e+01, 1.23026877e+01, 1.31825674e+01,
        1.41253754e+01, 1.51356125e+01, 1.62181010e+01, 1.73780083e+01, 1.86208714e+01, 1.99526231e+01, 2.13796209e+01,
        2.29086765e+01, 2.45470892e+01, 2.63026799e+01, 2.81838293e+01, 3.01995172e+01, 3.23593657e+01, 3.46736850e+01,
        3.71535229e+01, 3.98107171e+01, 4.26579519e+01, 4.57088190e+01, 4.89778819e+01, 5.24807460e+01, 5.62341325e+01,
        6.02559586e+01, 6.45654229e+01, 6.91830971e+01, 7.41310241e+01, 7.94328235e+01, 8.51138038e+01, 9.12010839e+01,
        9.77237221e+01, 1.04712855e+02, 1.12201845e+02, 1.20226443e+02, 1.28824955e+02, 1.38038426e+02, 1.47910839e+02,
        1.58489319e+02, 1.69824365e+02, 1.81970086e+02, 1.94984460e+02, 2.08929613e+02, 2.23872114e+02, 2.39883292e+02,
        2.57039578e+02,
    };

    fhElossVSLight = new TH2D("fhElossVSLight",
                              "Deposited Energy vs Generated Light",
                              bins.size() - 1,
                              bins.data(),
                              bins.size() - 1,
                              bins.data());
    fhElossVSLightLog = new TH2D("fhElossVSLightLog", "Deposited Energy vs Generated Light", 600, -3, 3, 600, -3, 3);

    fhThetaLight = new TH2D("fhThetaLight", "fhThetaLight", 200, -100, 100, 400, 0, 400);

    if (fIs3DTrackEnabled)
    {
        // XYZ -> ZXY (side view)

        fh3 = new TH3D("hMCTracks", "hMCTracks", 60, 1400, 1700, 50, -125, 125, 50, -125, 125);
        fh3->SetTitle("NeuLAND MCTracks");
        fh3->GetXaxis()->SetTitle("Z");
        fh3->GetYaxis()->SetTitle("X");
        fh3->GetZaxis()->SetTitle("Y");
        FairRootManager::Instance()->Register("NeulandMCMon", "MC Tracks in NeuLAND", fh3, kTRUE);

        fh3PNIP = new TH3D("h3DNPNIPS", "h3DNPNIPS", 60, 1400, 1700, 50, -125, 125, 50, -125, 125);
        fh3PNIP->SetTitle("NeuLAND 3D NPNIPS");
        fh3PNIP->GetXaxis()->SetTitle("Z");
        fh3PNIP->GetYaxis()->SetTitle("X");
        fh3PNIP->GetZaxis()->SetTitle("Y");
        FairRootManager::Instance()->Register("Neuland3DPNIP", "First interactions in NeuLAND", fh3PNIP, kTRUE);
    }

    return kSUCCESS;
}

void R3BNeulandMCMon::Exec(Option_t*)
{
    nEvents++;

    const auto npnips = fPrimaryNeutronInteractionPoints.Retrieve();
    const auto mcTracks = fMCTracks.Retrieve();
    const auto points = fNeulandPoints.Retrieve();

    for (const auto& mcTrack : mcTracks)
    {
        // Distribution of MC Track mother id's
        fhMotherIDs->Fill(mcTrack->GetMotherId());

        // Energy of primary Particles
        if (mcTrack->GetMotherId() == -1)
        {
            fhEPrimarys->Fill(GetKineticEnergy(mcTrack));
            fhEPrimarys2->Fill(mcTrack->GetEnergy() * 1000.);
        }

        // Energy of primary Neutrons
        if (IsPrimaryNeutron(mcTrack))
        {
            fhEPrimaryNeutrons->Fill(GetKineticEnergy(mcTrack));
        }

        // Remaining energy of neutrons after first interaction
        if (IsMotherPrimaryNeutron(mcTrack, mcTracks) && mcTrack->GetPdgCode() == 2112)
        {
            fhESecondaryNeutrons->Fill(GetKineticEnergy(mcTrack));
        }
    }

    fhnNPNIPs->Fill(npnips.size());
    for (const auto& npnip : npnips)
    {
        // WIP: ToF Calculation -> Should respect other origin than 0,0,0,0.
        const Double_t s2 = std::pow(npnip->GetX(), 2) + std::pow(npnip->GetY(), 2) + std::pow(npnip->GetZ(), 2); // cm²
        const Double_t v2 = s2 / std::pow(npnip->GetTime(), 2); // ns²

        const Double_t c2 = 898.75517873681758374898; // cm²/ns²
        const Double_t massNeutron = 939.565379;      // MeV/c²
        const Double_t ETimeOfFlight = massNeutron * ((1. / std::sqrt(1 - (v2 / c2))) - 1);

        auto mcTrack = mcTracks.at(npnip->GetTrackID());
        fhNPNIPsEToFVSTime->Fill(ETimeOfFlight, npnip->GetTime());
        fhMCToF->Fill(GetKineticEnergy(mcTrack) - ETimeOfFlight);
        fhNPNIPSrvsz->Fill(npnip->GetZ(), std::sqrt(std::pow(npnip->GetX(), 2) + std::pow(npnip->GetY(), 2)));
        fhNPNIPSxy->Fill(npnip->GetX(), npnip->GetY());
    }

    {
        std::map<Int_t, Double_t> EtotPDG;
        Double_t Etot = 0.;
        Double_t EtotPrim = 0.;

        for (const auto& point : points)
        {
            const R3BMCTrack* mcTrack = mcTracks.at(point->GetTrackID());

            Etot += point->GetLightYield() * 1000.;

            // Select tracks with a primary neutron mother
            if (IsMotherPrimaryNeutron(mcTrack, mcTracks))
            {

                // Total energy of non-neutron secondary particles where mother is a primary neutron
                if (mcTrack->GetPdgCode() != 2112)
                {
                    EtotPrim += point->GetLightYield() * 1000.;
                }

                // Distribution of secondary particles
                fhPDG->Fill(mcTrack->GetPdgCode());
                fhPrimaryDaughterIDs->Fill(point->GetTrackID());

                // Buld Histograms for each particle PDG if it doesn't exist
                if (!fhmEPdg[mcTrack->GetPdgCode()])
                {
                    TString name = TString("Light Yield of PID ") + TString::Itoa(mcTrack->GetPdgCode(), 10) +
                                   TString(" with a primary neutron mother");
                    fhmEPdg[mcTrack->GetPdgCode()] =
                        new TH1D("hE_PDG_" + TString::Itoa(mcTrack->GetPdgCode(), 10), name, 3000, 0, 3000);
                }
                // Get Energy py particle where the mother is a primary neutron
                fhmEPdg[mcTrack->GetPdgCode()]->Fill(point->GetLightYield() * 1000.); // point->GetEnergyLoss()*1000.);
            }                                                                         // end primary neutron mother

            // Sum energy per particle type per event
            if (!EtotPDG[mcTrack->GetPdgCode()])
            {
                EtotPDG[mcTrack->GetPdgCode()] = 0.;
            }
            EtotPDG[mcTrack->GetPdgCode()] += point->GetLightYield() * 1000.; // point->GetEnergyLoss()*1000.;

            fhElossVSLight->Fill(point->GetEnergyLoss() * 1000., point->GetLightYield() * 1000.);
            fhElossVSLightLog->Fill(std::log10(point->GetEnergyLoss() * 1000.),
                                    std::log10(point->GetLightYield() * 1000.));

            if (!fhmElossVSLightLogPdg[mcTrack->GetPdgCode()])
            {
                fhmElossVSLightLogPdg[mcTrack->GetPdgCode()] =
                    new TH2D("fhElossVSLightLog_" + TString::Itoa(mcTrack->GetPdgCode(), 10),
                             "Deposited Energy vs Generated Light",
                             600,
                             -3,
                             3,
                             600,
                             -3,
                             3);
            }
            fhmElossVSLightLogPdg[mcTrack->GetPdgCode()]->Fill(std::log10(point->GetEnergyLoss() * 1000.),
                                                               std::log10(point->GetLightYield() * 1000.));

            fhThetaLight->Fill(GetTheta(mcTrack), point->GetLightYield() * 1000.);
        }

        fhEtot->Fill(Etot);
        fhEtotPrim->Fill(EtotPrim);

        for (const auto& kv : EtotPDG)
        {
            if (!fhmEtotPdg[kv.first])
            {
                TString name = TString("Sum Light Yield of PID ") + TString::Itoa(kv.first, 10);
                fhmEtotPdg[kv.first] = new TH1D("fhEtotPDG_" + TString::Itoa(kv.first, 10), name, 3000, 0, 3000);
            }
            fhmEtotPdg[kv.first]->Fill(kv.second);

            if (!fhmEtotPdgRel[kv.first])
            {
                TString name = TString("Percent Light Yield of PID ") + TString::Itoa(kv.first, 10) +
                               TString(" to total Light Yield.");
                fhmEtotPdgRel[kv.first] = new TH1D("fhEtotPDGRel_" + TString::Itoa(kv.first, 10), name, 110, 0, 110);
            }
            if (Etot == 0)
            {
                fhmEtotPdgRel[kv.first]->Fill(0);
            }
            else
            {
                fhmEtotPdgRel[kv.first]->Fill(kv.second / Etot * 100.);
            }
        }
    }

    {
        std::vector<int> primaryNeutronIDs;
        std::map<int, std::vector<R3BMCTrack*>> tracksByPrimaryNeutronID;

        const Int_t nTracks = mcTracks.size();
        for (Int_t j = 0; j < nTracks; j++)
        {
            R3BMCTrack* track = mcTracks.at(j);
            if (track->GetMotherId() == -1 && track->GetPdgCode() == 2112)
            {
                primaryNeutronIDs.push_back(j);
            }
        }

        for (Int_t j = 0; j < nTracks; j++)
        {
            R3BMCTrack* track = mcTracks.at(j);
            for (const Int_t primaryNeutronID : primaryNeutronIDs)
            {
                if (track->GetMotherId() == primaryNeutronID)
                {
                    tracksByPrimaryNeutronID[primaryNeutronID].push_back(track);
                }
            }
        }

        for (const Int_t primaryNeutronID : primaryNeutronIDs)
        {
            auto& tracks = tracksByPrimaryNeutronID[primaryNeutronID];
            fhnProducts->Fill(tracks.size());

            fhnProductsCharged->Fill(std::accumulate(tracks.begin(),
                                                     tracks.end(),
                                                     0,
                                                     [](const int a, const R3BMCTrack* b)
                                                     {
                                                         if (b->GetPdgCode() == 2112 || b->GetPdgCode() == 22 ||
                                                             b->GetPdgCode() == 111)
                                                         {
                                                             return a;
                                                         }
                                                         else
                                                         {
                                                             return a + 1;
                                                         }
                                                     }));

            fhnSecondaryNeutrons->Fill(std::accumulate(tracks.begin(),
                                                       tracks.end(),
                                                       0,
                                                       [](const int c, const R3BMCTrack* t)
                                                       {
                                                           if (t->GetPdgCode() == 2112)
                                                           {
                                                               return c + 1;
                                                           }
                                                           return c;
                                                       }));

            fhnSecondaryProtons->Fill(std::accumulate(tracks.begin(),
                                                      tracks.end(),
                                                      0,
                                                      [](const int c, const R3BMCTrack* t)
                                                      {
                                                          if (t->GetPdgCode() == 2212)
                                                          {
                                                              return c + 1;
                                                          }
                                                          return c;
                                                      }));

            const Double_t sumProductEnergy = std::accumulate(tracks.begin(),
                                                              tracks.end(),
                                                              0.,
                                                              [](const Double_t a, const R3BMCTrack* b)
                                                              { return a + 1000. * (b->GetEnergy() - b->GetMass()); });

            fhSumProductEnergy->Fill(sumProductEnergy);

            if (fIsFullSimAnaEnabled)
            {
                std::map<int, int> mCountByProductPdg;

                for (const auto& track : tracks)
                {
                    if (mCountByProductPdg.find(track->GetPdgCode()) == mCountByProductPdg.end())
                    {
                        mCountByProductPdg[track->GetPdgCode()] = 1;
                    }
                    else
                    {
                        mCountByProductPdg[track->GetPdgCode()]++;
                    }

                    auto& hEnergy = fhmEnergyByProductPdg[track->GetPdgCode()];
                    if (hEnergy == nullptr)
                    {
                        hEnergy = new TH1D(TString::Itoa(track->GetPdgCode(), 10),
                                           TString::Itoa(track->GetPdgCode(), 10),
                                           1000,
                                           0,
                                           1000);
                    }
                    hEnergy->Fill(1000. * (track->GetEnergy() - track->GetMass()));

                    const auto l = lookup(track->GetPdgCode());
                    auto& hEnergyRed = fhmEnergyByProductPdgReduced[l];
                    if (hEnergyRed == nullptr)
                    {
                        hEnergyRed = new TH1D(l.c_str(), l.c_str(), 1000, 0, 1000);
                    }
                    hEnergyRed->Fill(1000. * (track->GetEnergy() - track->GetMass()));
                }

                for (const auto& kv : mCountByProductPdg)
                {
                    const int pdg = kv.first;
                    auto& hCount = fhmCountByProductPdg[pdg];
                    if (hCount == nullptr)
                    {
                        hCount = new TH1D(TString::Itoa(pdg, 10), TString::Itoa(pdg, 10), 10, 0, 10);
                    }
                    hCount->Fill(kv.second);

                    const auto l = lookup(kv.first);
                    auto& hCountRed = fhmCountByProductPdgReduced[l];
                    if (hCountRed == nullptr)
                    {
                        hCountRed = new TH1D(l.c_str(), l.c_str(), 10, 0, 10);
                    }
                    hCountRed->Fill(kv.second);
                }

                std::sort(tracks.begin(),
                          tracks.end(),
                          [](const R3BMCTrack* a, const R3BMCTrack* b) { return a->GetPdgCode() < b->GetPdgCode(); });

                // std::vector<int> filteredPdgCodes{ 22, -211, 211, 111 };
                std::vector<int> filteredPdgCodes{ 22, 111 };
                const TString reaction = std::accumulate(
                    tracks.begin(),
                    tracks.end(),
                    TString(),
                    [&](TString s, const R3BMCTrack* b) -> TString
                    {
                        if (std::find(filteredPdgCodes.begin(), filteredPdgCodes.end(), b->GetPdgCode()) !=
                            filteredPdgCodes.end())
                        {
                            return s;
                        }
                        return s += "_" + TString::Itoa(b->GetPdgCode(), 10);
                    });

                if (fhmEnergyByReaction[reaction] == nullptr)
                {
                    fhmEnergyByReaction[reaction] = new TH1D(reaction, reaction, 1000, 0, 1000);
                }
                fhmEnergyByReaction[reaction]->Fill(sumProductEnergy);

                for (const auto& track : tracks)
                {
                    auto& hEnergy = fhmmEnergyByReactionByProductPdg[reaction][track->GetPdgCode()];
                    if (hEnergy == nullptr)
                    {
                        hEnergy = new TH1D(TString::Itoa(track->GetPdgCode(), 10) + reaction,
                                           TString::Itoa(track->GetPdgCode(), 10) + reaction,
                                           1000,
                                           0,
                                           1000);
                    }
                    hEnergy->Fill(1000. * (track->GetEnergy() - track->GetMass()));
                }
            }
        }
    }

    {
        ROOT::Math::PxPyPzEVector p4;
        double m0 = 0;
        for (const auto track : mcTracks)
        {
            if (track->GetMotherId() != -1)
            {
                break;
            }
            p4 += track->GetFourMomentum();
            m0 += track->GetMass();
        }
        fhErelMC->Fill((p4.M() - m0) * 1e6);
    }

    if (fIs3DTrackEnabled)
    {
        fh3->Reset("ICES");
        for (const auto& point : points)
        {
            if (point->GetLightYield() > 0)
            {
                const auto pos = point->GetPosition();
                fh3->Fill(pos.Z(), pos.X(), pos.Y(), point->GetEnergyLoss() * 1000.);
            }
        }

        fh3PNIP->Reset("ICES");
        for (const auto& npnip : npnips)
        {
            fh3PNIP->Fill(npnip->GetZ(), npnip->GetX(), npnip->GetY(), npnip->GetTime());
        }
    }
}

void R3BNeulandMCMon::Finish()
{
    TDirectory* tmp = gDirectory;
    FairRootManager::Instance()->GetOutFile()->cd();

    gDirectory->mkdir("NeulandMCMon");
    gDirectory->cd("NeulandMCMon");

    fhPDG->Write();
    fhEPrimarys->Write();
    fhEPrimarys2->Write();
    fhEPrimaryNeutrons->Write();
    fhErelMC->Write();
    fhEtot->Write();
    fhEtotPrim->Write();
    fhESecondaryNeutrons->Write();
    fhMotherIDs->Write();
    fhPrimaryDaughterIDs->Write();
    fhMCToF->Write();
    fhNPNIPsEToFVSTime->Write();
    fhNPNIPSrvsz->Write();
    fhNPNIPSxy->Write();
    fhnNPNIPs->Write();
    fhThetaLight->Write();
    fhElossVSLight->Write();
    fhElossVSLightLog->Write();

    gDirectory = tmp;
    gDirectory->cd("NeulandMCMon");
    gDirectory->mkdir("LightYield");
    gDirectory->cd("LightYield");

    gDirectory->mkdir("LightYieldByProductPdg");
    gDirectory->cd("LightYieldByProductPdg");
    for (const auto& kv : fhmEPdg)
    {
        kv.second->Write();
    }
    for (const auto& kv : fhmElossVSLightLogPdg)
    {
        kv.second->Write();
    }
    gDirectory->cd("..");

    gDirectory->mkdir("SumLightYieldByProductPdg");
    gDirectory->cd("SumLightYieldByProductPdg");
    for (const auto& kv : fhmEtotPdg)
    {
        kv.second->Write();
    }
    gDirectory->cd("..");

    gDirectory->mkdir("RelSumLightYieldByProductPdg");
    gDirectory->cd("RelSumLightYieldByProductPdg");
    for (const auto& kv : fhmEtotPdgRel)
    {
        kv.second->Write();
    }
    gDirectory->cd("..");

    gDirectory = tmp;
    gDirectory->cd("NeulandMCMon");
    gDirectory->mkdir("Products");
    gDirectory->cd("Products");

    fhnProducts->Write();
    fhnProductsCharged->Write();
    fhSumProductEnergy->Write();
    fhnSecondaryProtons->Write();
    fhnSecondaryNeutrons->Write();

    writeout(fhmEnergyByProductPdg, "EnergyByProduct");
    writeout(fhmEnergyByProductPdgReduced, "EnergyByProductReduced");
    writeout(fhmCountByProductPdg, "CountByProduct", nEvents);
    writeout(fhmCountByProductPdgReduced, "CountByProductReduced", nEvents);
    writeout(fhmEnergyByReaction, "EnergyByReaction");

    if (!fhmmEnergyByReactionByProductPdg.empty())
    {
        gDirectory->mkdir("EnergyByReactionByProductPdg");
        gDirectory->cd("EnergyByReactionByProductPdg");
        for (const auto& kv : fhmmEnergyByReactionByProductPdg)
        {
            gDirectory->mkdir(kv.first);
            gDirectory->cd(kv.first);
            for (const auto& kv2 : kv.second)
            {
                kv2.second->Write();
            }
            gDirectory->cd("..");
        }
        gDirectory->cd("..");
    }

    gDirectory = tmp;
}

ClassImp(R3BNeulandMCMon)
