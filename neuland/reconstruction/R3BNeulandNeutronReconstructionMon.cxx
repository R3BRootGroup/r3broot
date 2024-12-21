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

#include "R3BNeulandNeutronReconstructionMon.h"

#include <algorithm>
#include <cmath>
#include <functional>
#include <iostream>
#include <numeric>
#include <utility>

#include "TClonesArray.h"
#include "TDirectory.h"
#include "TH1D.h"
#include "TH2D.h"
#include <TFile.h>

#include "FairLogger.h"
#include "FairMCPoint.h"
#include "FairRootManager.h"
#include "FairRtdbRun.h"
#include "FairRuntimeDb.h"

#include "R3BMCTrack.h"
#include "R3BNeulandNeutron.h"
#include "R3BNeulandPoint.h"

static const Double_t c2 = 898.75517873681758374; // cm²/ns²
static const Double_t massNeutron = 939.565379;   // MeV/c²

Double_t Distance(const R3BNeulandNeutron& nn, const FairMCPoint& mc)
{
    TVector3 v;
    mc.Position(v);
    v -= nn.GetPosition();
    return std::sqrt(v.Dot(v));
}

Double_t Score(const std::vector<std::pair<R3BNeulandNeutron, FairMCPoint>>& combination)
{
    return std::accumulate(combination.begin(),
                           combination.end(),
                           0.,
                           [](const Double_t sum, const std::pair<R3BNeulandNeutron, FairMCPoint>& pair)
                           { return sum + Distance(pair.first, pair.second); });
}

template <typename T, typename U>
std::vector<std::vector<std::pair<T, U>>> GetAllCombinations(std::vector<T> ts /* Copy! */,
                                                             std::vector<U> us /* Copy! */,
                                                             std::function<bool(const U&, const U&)> comparator)
{
    std::vector<std::vector<std::pair<T, U>>> out;

    // Bring both inputs up to the same length
    if (ts.size() < us.size())
    {
        ts.resize(us.size());
    }
    if (us.size() < ts.size())
    {
        us.resize(ts.size());
    }

    std::sort(us.begin(), us.end(), comparator);
    do
    {
        std::vector<std::pair<T, U>> tmp;
        for (size_t i = 0; i < ts.size(); i++)
        {
            tmp.push_back({ ts.at(i), us.at(i) });
        }
        out.push_back(std::move(tmp));
    } while (std::next_permutation(us.begin(), us.end(), comparator));

    return out;
}

R3BNeulandNeutronReconstructionMon::R3BNeulandNeutronReconstructionMon(const TString input, const TString output)
    : FairTask("R3B Neuland Neutron Reconstruction Evaluation")
    , fInput(input)
    , fOutput(output)
{
    LOG(info) << "Using R3B Neuland Neutron Reconstruction Evaluation";
}

R3BNeulandNeutronReconstructionMon::~R3BNeulandNeutronReconstructionMon() {}

InitStatus R3BNeulandNeutronReconstructionMon::Init()
{
    FairRootManager* ioman = FairRootManager::Instance();
    if (!ioman)
    {
        LOG(fatal) << "R3BNeulandNeutronReconstructionMon::Init: No FairRootManager";
        return kFATAL;
    }

    // Set Input: TClonesArray of R3BNeulandNeutrons
    if (dynamic_cast<TClonesArray*>(ioman->GetObject(fInput)) == nullptr)
    {
        LOG(fatal) << "R3BNeulandNeutronReconstructionMon::Init No " << fInput << "!";
        return kFATAL;
    }
    if (!TString((dynamic_cast<TClonesArray*>(ioman->GetObject(fInput)))->GetClass()->GetName())
             .EqualTo("R3BNeulandNeutron"))
    {
        LOG(fatal) << "R3BNeulandNeutronReconstructionMon::Init Branch " << fInput
                   << " does not contain "
                      "R3BNeulandNeutrons!";
        return kFATAL;
    }
    fReconstructedNeutrons = dynamic_cast<TClonesArray*>(ioman->GetObject(fInput));

    // Set Input: TClonesArray of FairMCPoints
    if (dynamic_cast<TClonesArray*>(ioman->GetObject("NeulandPrimaryPoints")) == nullptr)
    {
        LOG(fatal) << "R3BNeulandNeutronReconstructionMon::Init No NeulandPrimaryPoints!";
        return kFATAL;
    }
    if (!TString((dynamic_cast<TClonesArray*>(ioman->GetObject("NeulandPrimaryPoints")))->GetClass()->GetName())
             .EqualTo("R3BNeulandPoint"))
    {
        LOG(fatal) << "R3BNeulandNeutronReconstructionMon::Init Branch NeulandPrimaryPoints "
                      "does not contain R3BNeulandPoint!";
        return kFATAL;
    }
    fPrimaryNeutronInteractionPoints = dynamic_cast<TClonesArray*>(ioman->GetObject("NeulandPrimaryPoints"));

    // Set Input: TClonesArray of R3BMCTrack
    if (dynamic_cast<TClonesArray*>(ioman->GetObject("MCTrack")) == nullptr)
    {
        LOG(fatal) << "R3BNeulandNeutronReconstructionMon::Init No MCTrack!";
        return kFATAL;
    }
    if (!TString((dynamic_cast<TClonesArray*>(ioman->GetObject("MCTrack")))->GetClass()->GetName())
             .EqualTo("R3BMCTrack"))
    {
        LOG(fatal) << "R3BNeulandNeutronReconstructionMon::Init Branch MCTrack "
                      "does not contain FairMCPoints!";
        return kFATAL;
    }
    fMCTracks = dynamic_cast<TClonesArray*>(ioman->GetObject("MCTrack"));

    TH1::AddDirectory(kFALSE);

    fhScore = new TH1D("fhScore", "Neuland Neutron Reconstruction Score (lower is better)", 5000, 0, 5000);
    fhCountN = new TH1D("fhCountN", "Number of reconstructed Neutrons", 11, -0.5, 10.5);
    fhCountNdiff =
        new TH1D("fhCountNdiff", "Number of reacted primary Neutrons - Number of reconstructed Neutrons", 11, -5, 5);
    fhEdiff = new TH1D("fhEdiff", "Energy of primary Neutron - Energy of reconstructed Neutron", 2001, -1000, 1000);

    fhErel = new TH1D("fhErel", "fhErel", 5000, 0, 5000);
    fhErelMC = new TH1D("fhErelMC", "fhErelMC", 5000, 0, 5000);
    fhErelVSnNreco = new TH2D("fhErelVSnNreco", "fhErelVSnNreco", 5000, 0, 5000, 11, -0.5, 10.5);

    fhErelVSnNrecoNPNIPs = new TH2D("fhErelVSnNrecoNPNIPs", "fhErelVSnNrecoNPNIPs", 5000, 0, 5000, 11, -0.5, 10.5);

    fhNreacNreco = new TH2D("fhNreacNreco", "fhNreacNreco", 11, -0.5, 10.5, 11, -0.5, 10.5);
    fhNreacNreco->GetXaxis()->SetTitle("nReac");
    fhNreacNreco->GetYaxis()->SetTitle("nReco");

    return kSUCCESS;
}

void R3BNeulandNeutronReconstructionMon::Exec(Option_t*)
{
    const Int_t nReconstructedNeutrons = fReconstructedNeutrons->GetEntries();
    std::vector<R3BNeulandNeutron> neutrons;
    neutrons.reserve(nReconstructedNeutrons);
    for (Int_t i = 0; i < nReconstructedNeutrons; i++)
    {
        neutrons.push_back(*(dynamic_cast<R3BNeulandNeutron*>(fReconstructedNeutrons->At(i))));
    }

    const Int_t nNPNIPS = fPrimaryNeutronInteractionPoints->GetEntries();
    std::vector<R3BNeulandPoint> npnips;
    npnips.reserve(nReconstructedNeutrons);
    for (Int_t i = 0; i < nNPNIPS; i++)
    {
        npnips.push_back(*(dynamic_cast<R3BNeulandPoint*>(fPrimaryNeutronInteractionPoints->At(i))));
    }
    fhCountN->Fill(nReconstructedNeutrons);
    fhCountNdiff->Fill((Int_t)nNPNIPS - (Int_t)nReconstructedNeutrons);
    fhNreacNreco->Fill(nNPNIPS, nReconstructedNeutrons);

    // Find the lowes score by brute force creating all combinations of MC-known and reconstructed neutron hits.
    // Score is the sum of distances. Note that detecting a wrong amount of neutrons will create a huge score
    // WiP. Slow because of all the copying
    /*const auto allCombinations =
        GetAllCombinations<R3BNeulandNeutron, FairMCPoint>(neutrons,
                                                           npnips,
                                                           [](const FairMCPoint& a, const FairMCPoint& b)
                                                           {
                                                               return a.GetTime() < b.GetTime();
                                                           });

    const auto bestCombination = *std::min_element(allCombinations.begin(),
                                                   allCombinations.end(),
                                                   [](const std::vector<std::pair<R3BNeulandNeutron, FairMCPoint>>& a,
                                                      const std::vector<std::pair<R3BNeulandNeutron, FairMCPoint>>& b)
                                                   {
                                                       return Score(a) < Score(b);
                                                   });

    const Double_t minScore = Score(bestCombination);
    fhScore->Fill(minScore);

    for (const auto& pair : bestCombination)
    {
        if (pair.second.GetEnergyLoss() > 0 || pair.first.GetEkin() > 0)
        {
            // TODO: Only for massNeutron
            if (pair.second.GetEnergyLoss() > 0)
            {
                fhEdiff->Fill(pair.second.GetEnergyLoss() * 1000. - massNeutron - pair.first.GetEkin());
            }
            else
            {
                fhEdiff->Fill(-pair.first.GetEkin());
            }
        }
    }*/

    // Erel
    // if (nReconstructedNeutrons == 0)
    {
        TLorentzVector p4_reco;
        Double_t m0_reco = 0;

        TLorentzVector p4_mc;
        Double_t m0_mc = 0;

        TLorentzVector p4_npnips;
        Double_t m0_npnips = 0;

        for (const auto& neutron : neutrons)
        {
            p4_reco += TLorentzVector(neutron.GetP(), neutron.GetEtot());
            m0_reco += massNeutron;
        }

        for (const auto& npnip : npnips)
        {
            // Workaround mom = 0
            R3BMCTrack* MCTrack = dynamic_cast<R3BMCTrack*>(fMCTracks->At(npnip.GetTrackID()));
            p4_npnips += TLorentzVector(MCTrack->GetPx() * 1000.,
                                        MCTrack->GetPy() * 1000.,
                                        MCTrack->GetPz() * 1000.,
                                        MCTrack->GetEnergy() * 1000.);
            m0_npnips += MCTrack->GetMass() * 1000.;
        }

        R3BMCTrack* MCTrack;
        const Int_t nMCTracks = fMCTracks->GetEntries();
        for (Int_t i = 0; i < nMCTracks; i++)
        {
            MCTrack = dynamic_cast<R3BMCTrack*>(fMCTracks->At(i));
            // We are only interested in primary particles, these should come first,
            // thus we can stop once the Mother Id != -1
            if (MCTrack->GetMotherId() != -1)
            {
                // break;
                continue;
            }

            // For E_rel with reconstructed neutrons, only sum up non-neutron primary particles
            // This information has to come from the other detectors later.
            if (MCTrack->GetPdgCode() != 2112)
            {
                p4_reco += TLorentzVector(MCTrack->GetPx() * 1000.,
                                          MCTrack->GetPy() * 1000.,
                                          MCTrack->GetPz() * 1000.,
                                          MCTrack->GetEnergy() * 1000.);
                m0_reco += MCTrack->GetMass() * 1000.;

                p4_npnips += TLorentzVector(MCTrack->GetPx() * 1000.,
                                            MCTrack->GetPy() * 1000.,
                                            MCTrack->GetPz() * 1000.,
                                            MCTrack->GetEnergy() * 1000.);
                m0_npnips += MCTrack->GetMass() * 1000.;
            }

            // For comparison, calculate E_rel for all primary particles from MC
            p4_mc += TLorentzVector(MCTrack->GetPx() * 1000.,
                                    MCTrack->GetPy() * 1000.,
                                    MCTrack->GetPz() * 1000.,
                                    MCTrack->GetEnergy() * 1000.);
            m0_mc += MCTrack->GetMass() * 1000.;
        }

        const Double_t minv = p4_reco.Mag() - m0_reco;
        fhErel->Fill(minv * 1000.);
        fhErelVSnNreco->Fill(minv * 1000., nReconstructedNeutrons);

        if (fhmErelnReco[nReconstructedNeutrons] == nullptr)
        {
            fhmErelnReco[nReconstructedNeutrons] = new TH1D("fhErel" + TString::Itoa(nReconstructedNeutrons, 10),
                                                            "fhErel" + TString::Itoa(nReconstructedNeutrons, 10),
                                                            5000,
                                                            0,
                                                            5000);
        }
        fhmErelnReco.at(nReconstructedNeutrons)->Fill(minv * 1000.);

        const Double_t minv_mc = p4_mc.Mag() - m0_mc;
        fhErelMC->Fill(minv_mc * 1000.);

        const Double_t minv_npnips = p4_npnips.Mag() - m0_npnips;
        fhErelVSnNrecoNPNIPs->Fill(minv_npnips * 1000., nNPNIPS);
    }
}

void R3BNeulandNeutronReconstructionMon::Finish()
{
    TDirectory* tmp = gDirectory;
    FairRootManager::Instance()->GetOutFile()->cd();

    gDirectory->mkdir(fOutput);
    gDirectory->cd(fOutput);

    fhCountN->Write();
    fhCountNdiff->Write();
    fhScore->Write();
    fhEdiff->Write();
    fhErel->Write();
    fhErelMC->Write();
    fhErelVSnNreco->Write();
    fhErelVSnNrecoNPNIPs->Write();
    fhNreacNreco->Write();

    for (auto& nh : fhmErelnReco)
    {
        nh.second->Write();
    }

    gDirectory = tmp;
}

ClassImp(R3BNeulandNeutronReconstructionMon)
