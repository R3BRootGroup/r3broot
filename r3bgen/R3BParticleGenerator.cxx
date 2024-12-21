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

#include "R3BParticleGenerator.h"
#include "R3BDistribution2D.h"
#include "R3BDistribution3D.h"

#include "FairIon.h"
#include "FairLogger.h"
#include "FairParticle.h"
#include "FairPrimaryGenerator.h"

#include "TDatabasePDG.h"
#include "TMath.h"
#include "TVector3.h"

R3BParticleGenerator::R3BParticleGenerator(int pdgCode, unsigned int seed)
    : fVertex_cm(R3BDistribution3D::Delta(0, 0, 0))
    , fSpread_mRad({ R3BDistribution1D::Delta(0), R3BDistribution1D::Flat(0, 2 * TMath::Pi() * 1e3) })
    , fEnergy_AMeV(R3BDistribution1D::Delta(500))
    , fPDGCode(0)
    , fRNG(seed)
{
    AddParticle(pdgCode);
}

R3BParticleGenerator::R3BParticleGenerator(unsigned int seed)
    : fVertex_cm(R3BDistribution3D::Delta(0, 0, 0))
    , fSpread_mRad({ R3BDistribution1D::Delta(0), R3BDistribution1D::Flat(0, 2 * TMath::Pi() * 1e3) })
    , fEnergy_AMeV(R3BDistribution1D::Delta(500))
    , fPDGCode(0)
    , fRNG(seed)
{
}

bool R3BParticleGenerator::Init()
{
    if (fPDGCode == 0)
        LOG(fatal) << "R3BParticleGenerator::Init: No Particle added!";

    return true;
}

bool R3BParticleGenerator::ReadEvent(FairPrimaryGenerator* primGen)
{
    const auto vertex_cm = Beam.GetVertexDistribution().GetRandomValues({ fRNG.Rndm(), fRNG.Rndm(), fRNG.Rndm() });
    const auto theta_phi_mRad = Beam.GetSpreadDistribution().GetRandomValues({ fRNG.Rndm(), fRNG.Rndm() });
    const auto beamBeta = Beam.GetBetaDistribution().GetRandomValues({ fRNG.Rndm() })[0];

    const auto beamGamma = 1. / sqrt(1. - beamBeta * beamBeta);

    TVector3 momentum(0., 0., beamBeta * beamGamma * fMass_GeV);
    momentum.RotateX(theta_phi_mRad[0] * 1e-3);
    momentum.RotateZ(theta_phi_mRad[1] * 1e-3);

    const auto totalEnergy = sqrt(momentum.Mag() * momentum.Mag() + fMass_GeV * fMass_GeV);

    primGen->AddTrack(fPDGCode,
                      momentum[0],
                      momentum[1],
                      momentum[2],
                      vertex_cm[0],
                      vertex_cm[1],
                      vertex_cm[2],
                      -1,
                      true,
                      totalEnergy);

    return kTRUE;
}

void R3BParticleGenerator::addParticle(const int pdgCode, const double mass)
{
    if (fPDGCode != 0)
    {
        // one particle was already added
        LOG(fatal) << "R3BParticleGenerator: Trying to add a second Particle. Use a second generator!";
    }

    fPDGCode = pdgCode;
    fMass_GeV = mass;
}

ClassImp(R3BParticleGenerator)
