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

#include "R3BCosmicGenerator.h"
#include "R3BDistribution1D.h"
#include "R3BDistribution2D.h"
#include "R3BDistribution3D.h"

#include "FairLogger.h"
#include "FairPrimaryGenerator.h"

#include "TDatabasePDG.h"
#include "TMath.h"
#include "TParticlePDG.h"

#include <vector>

// Data from doi.org/10.1016/0029-554X(81)91039-9
const std::vector<Double_t> _Distribution_EnergyValues = { .004, .01, .02, .04, .1,  .2,  .4,
                                                           1.,   2.,  4.,  10., 20., 40., 100. };
const std::vector<Double_t> _Distribution_PropabilityValues = { 5.89E-7, 1.05E-6, 1.63E-6, 2.34E-6, 3.13E-6,
                                                                3.89E-6, 4.18E-6, 2.02E-6, 9.11E-7, 3.81E-7,
                                                                1.11E-7, 3.00E-8, 3.40E-9, 9.01E-11 };

constexpr Double_t Sqr(const Double_t val) { return val * val; }

Double_t Sin2(const Double_t val) { return Sqr(sin(val)); }

R3BCosmicGenerator::R3BCosmicGenerator(const Int_t pdgid, const Int_t mult, const UInt_t seed)
    : fVertexPositionDist_cm(R3BDistribution3D::Delta(0, 0, 0))
    , fAngularDist_Rad(
          { R3BDistribution1D::Function(Sin2, 0, TMath::Pi()), R3BDistribution1D::Flat(0, 2 * TMath::Pi()) })
    , fRngGen(seed)
    , fPDGType(pdgid)
    , fMult(mult)
{
    SetEnergyRange_GeV();
}

void R3BCosmicGenerator::SetEnergyRange_GeV(const Double_t minEnergy, const Double_t maxEnergy)
{
    if (minEnergy < _Distribution_EnergyValues[0] || maxEnergy > _Distribution_EnergyValues.back())
        LOG(fatal) << "R3BCosmicGenerator::SetEnergyRange_GeV : Bad Energy-Range: Maximum Range is "
                   << _Distribution_EnergyValues[0] << " < E < " << _Distribution_EnergyValues.back();

    TGraph g(_Distribution_EnergyValues.size(), &_Distribution_EnergyValues[0], &_Distribution_PropabilityValues[0]);

    // Get Log-Log integrated distribution
    fEnergyDist_GeV = R3BDistribution1D::DataLogLog(g, minEnergy, maxEnergy);
}

Bool_t R3BCosmicGenerator::Init()
{
    // Check for particle type
    TDatabasePDG* pdgBase = TDatabasePDG::Instance();
    TParticlePDG* particle = pdgBase->GetParticle(fPDGType);
    if (!particle)
        LOG(fatal) << "R3BCosmicGenerator: PDG code " << fPDGType << " not defined.";
    fPDGMass = particle->Mass();
    return kTRUE;
}

Bool_t R3BCosmicGenerator::ReadEvent(FairPrimaryGenerator* primGen)
{
    for (Int_t k = 0; k < fMult; k++)
    {
        const auto energy = fEnergyDist_GeV.GetRandomValues({ fRngGen.Rndm() })[0]; // in GeV
        const auto angDist = fAngularDist_Rad.GetRandomValues({ fRngGen.Rndm(), fRngGen.Rndm() });
        const auto vertex = fVertexPositionDist_cm.GetRandomValues({ fRngGen.Rndm(), fRngGen.Rndm(), fRngGen.Rndm() });

        const auto ptot = TMath::Sqrt(energy * energy + 2 * fPDGMass * energy);
        const auto px = ptot * TMath::Sin(angDist[0]) * TMath::Cos(angDist[1]);
        const auto pz = ptot * TMath::Sin(angDist[0]) * TMath::Sin(angDist[1]);
        const auto py = ptot * TMath::Cos(angDist[0]);

        primGen->AddTrack(fPDGType, px, py, pz, vertex[0], vertex[1], vertex[2]);
    }

    return kTRUE;
}

ClassImp(R3BCosmicGenerator)
