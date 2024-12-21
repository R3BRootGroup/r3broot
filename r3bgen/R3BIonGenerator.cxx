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

// -------------------------------------------------------------------------
// -----            Based on FairIonGenerator source file              -----
// -----            Created 18/03/11  by M. Labiche                    -----
// -------------------------------------------------------------------------

#include "R3BIonGenerator.h"

#include "FairIon.h"
#include "FairLogger.h"
#include "FairPrimaryGenerator.h"
#include "FairRunSim.h"

#include "TDatabasePDG.h"
#include "TMath.h"
#include "TObjArray.h"
#include "TParticle.h"
#include "TParticlePDG.h"
#include "TRandom.h"
#include "TString.h"

R3BIonGenerator::R3BIonGenerator(UInt_t seed)
    : fMult(0)
    , fIon(nullptr)
    , fRNG(seed)
{
}

R3BIonGenerator::R3BIonGenerator(const Char_t* ionName, Int_t mult, Double_t momentum_AGeV_per_c, UInt_t seed)
    : fMult(mult)
    , fIon(nullptr)
    , fRNG(seed)
{
    fIon = dynamic_cast<FairIon*>(FairRunSim::Instance()->GetUserDefIons()->FindObject(ionName));

    if (!fIon)
        LOG(fatal) << "R3BIonGenerator: Ion is not defined!";

    const auto mass = fIon->GetMass();
    const auto totalMomentum_GeV_per_c = fIon->GetA() * momentum_AGeV_per_c;
    const auto beta = momentum_AGeV_per_c / sqrt(mass * mass + totalMomentum_GeV_per_c * totalMomentum_GeV_per_c);
    Beam.SetBetaDistribution(R3BDistribution1D::Delta(beta));
}

R3BIonGenerator::R3BIonGenerator(Int_t z, Int_t a, Int_t q, Int_t mult, Double_t momentum_AGeV_per_c, UInt_t seed)
    : fMult(mult)
    , fIon(NULL)
    , fRNG(seed)
{
    fIon = new FairIon(TString::Format("FairIon_%d_%d_%d", z, a, q), z, a, q);

    auto run = FairRunSim::Instance();
    if (!run)
        LOG(fatal) << "FairIonGenerator: No FairRun instantised!";

    const auto mass = fIon->GetMass();
    const auto totalMomentum_GeV_per_c = fIon->GetA() * momentum_AGeV_per_c;
    const auto beta = momentum_AGeV_per_c / sqrt(mass * mass + totalMomentum_GeV_per_c * totalMomentum_GeV_per_c);
    Beam.SetBetaDistribution(R3BDistribution1D::Delta(beta));

    run->AddNewIon(fIon);
}

R3BIonGenerator::R3BIonGenerator(const R3BIonGenerator& right)
    : Beam(right.Beam)
    , fMult(right.fMult)
    , fIon(right.fIon)
    , fRNG(right.fRNG)
{
}

R3BIonGenerator::~R3BIonGenerator() {}

void R3BIonGenerator::SetExcitationEnergy(Double_t eExc) { fIon->SetExcEnergy(eExc); }

void R3BIonGenerator::SetMass(Double_t mass) { fIon->SetMass(mass); }

Bool_t R3BIonGenerator::Init()
{
    TParticlePDG* thisPart = TDatabasePDG::Instance()->GetParticle(fIon->GetName());
    if (!thisPart)
        LOG(fatal) << "FairIonGenerator: Ion " << fIon->GetName() << " not found in database!";

    return kTRUE;
}

Bool_t R3BIonGenerator::ReadEvent(FairPrimaryGenerator* primGen)
{
    TParticlePDG* thisPart = TDatabasePDG::Instance()->GetParticle(fIon->GetName());
    const auto pdgCode = thisPart->PdgCode();

    const auto vertex_cm = Beam.GetVertexDistribution().GetRandomValues({ fRNG.Rndm(), fRNG.Rndm(), fRNG.Rndm() });
    const auto theta_phi_mRad = Beam.GetSpreadDistribution().GetRandomValues({ fRNG.Rndm(), fRNG.Rndm() });
    const auto beamBeta = Beam.GetBetaDistribution().GetRandomValues({ fRNG.Rndm() })[0];

    const auto beamGamma = 1. / sqrt(1. - beamBeta * beamBeta);

    TVector3 momentum(0., 0., beamBeta * beamGamma * fIon->GetMass());
    momentum.RotateX(theta_phi_mRad[0] * 1e-3);
    momentum.RotateZ(theta_phi_mRad[1] * 1e-3);

    const auto totalEnergy = sqrt(momentum.Mag() * momentum.Mag() + fIon->GetMass() * fIon->GetMass());

    for (Int_t i = 0; i < fMult; i++)
        primGen->AddTrack(pdgCode,
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

ClassImp(R3BIonGenerator)
