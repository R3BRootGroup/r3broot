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

#include "R3BCALIFATestGenerator.h"

#include "FairLogger.h"
#include "FairPrimaryGenerator.h"

#include "TDatabasePDG.h"
#include "TMath.h"
#include "TParticlePDG.h"
#include "TRandom.h"

R3BCALIFATestGenerator::R3BCALIFATestGenerator()
    : fPDGType(0)
    , fMult(0)
    , fPDGMass(0)
    , fPtMin(0)
    , fPtMax(0)
    , fPhiMin(0)
    , fPhiMax(0)
    , fEtaMin(0)
    , fEtaMax(0)
    , fYMin(0)
    , fYMax(0)
    , fPMin(0)
    , fPMax(0)
    , fThetaMin(0)
    , fThetaMax(0)
    , fX(0)
    , fY(0)
    , fZ(0)
    , fX1(0)
    , fY1(0)
    , fZ1(0)
    , fX2(0)
    , fY2(0)
    , fZ2(0)
    , fEtaRangeIsSet(0)
    , fYRangeIsSet(0)
    , fThetaRangeIsSet(0)
    , fCosThetaIsSet(0)
    , fPtRangeIsSet(0)
    , fPRangeIsSet(0)
    , fPointVtxIsSet(0)
    , fBoxVtxIsSet(0)
    , fDebug(0)
    , fGammasDefinedInNuclearDecay(0)
    , fBetaOfEmittingFragment(0)
    , fGammaFactor(1)
    , fLorentzBoostIsSet(0)
    , fNuclearDecayChainIsSet(0)
{
    // Default constructor
}

R3BCALIFATestGenerator::R3BCALIFATestGenerator(Int_t pdgid, Int_t mult)
    : fPDGType(pdgid)
    , fMult(mult)
    , fPDGMass(0)
    , fPtMin(0)
    , fPtMax(0)
    , fPhiMin(0)
    , fPhiMax(0)
    , fEtaMin(0)
    , fEtaMax(0)
    , fYMin(0)
    , fYMax(0)
    , fPMin(0)
    , fPMax(0)
    , fThetaMin(0)
    , fThetaMax(0)
    , fX(0)
    , fY(0)
    , fZ(0)
    , fX1(0)
    , fY1(0)
    , fZ1(0)
    , fX2(0)
    , fY2(0)
    , fZ2(0)
    , fEtaRangeIsSet(0)
    , fYRangeIsSet(0)
    , fThetaRangeIsSet(0)
    , fCosThetaIsSet(0)
    , fPtRangeIsSet(0)
    , fPRangeIsSet(0)
    , fPointVtxIsSet(0)
    , fBoxVtxIsSet(0)
    , fDebug(0)
    , fGammasDefinedInNuclearDecay(0)
    , fBetaOfEmittingFragment(0)
    , fGammaFactor(1)
    , fLorentzBoostIsSet(0)
    , fNuclearDecayChainIsSet(0)
{
    // Constructor. Set default kinematics limits
    SetPhiRange();
}

Bool_t R3BCALIFATestGenerator::Init()
{
    // Initialize generator

    if (fPhiMax - fPhiMin > 360)
        LOG(fatal) << "Init(): R3BCALIFATestGenerator: phi range is too wide: " << fPhiMin << "<phi<" << fPhiMax;
    if (fPRangeIsSet && fPtRangeIsSet)
        LOG(fatal) << "Init(): R3BCALIFATestGenerator: Cannot set P and Pt ranges simultaneously";
    if (fPRangeIsSet && fYRangeIsSet)
        LOG(fatal) << "Init(): R3BCALIFATestGenerator: Cannot set P and Y ranges simultaneously";
    if ((fThetaRangeIsSet && fYRangeIsSet) || (fThetaRangeIsSet && fEtaRangeIsSet) || (fYRangeIsSet && fEtaRangeIsSet))
        LOG(fatal) << "Init(): R3BCALIFATestGenerator: Cannot set Y, Theta or Eta ranges simultaneously";
    if (fPointVtxIsSet && fBoxVtxIsSet)
        LOG(fatal) << "Init(): R3BCALIFATestGenerator: Cannot set point and box vertices simultaneously";

    // CALIFA specifics
    if (fBetaOfEmittingFragment > 1)
        LOG(fatal) << "Init(): R3BCALIFATestGenerator: beta of fragment larger than 1!";

    Double32_t sumBranchingRatios = 0;
    for (Int_t i = 0; i < fGammasDefinedInNuclearDecay; i++)
    {
        if (fGammaBranchingRatios[i] > 1)
            LOG(fatal) << "Init(): R3BCALIFATestGenerator: gamma branching ratio in position " << i
                       << " larger than 1!";
        sumBranchingRatios += fGammaBranchingRatios[i];
    }
    if (sumBranchingRatios > 1)
        LOG(fatal) << "Init(): R3BCALIFATestGenerator: gamma branching ratio sum larger than 1!";

    // Check for particle type
    TDatabasePDG* pdgBase = TDatabasePDG::Instance();
    TParticlePDG* particle = pdgBase->GetParticle(fPDGType);
    if (!particle)
        LOG(fatal) << "R3BCALIFATestGenerator: PDG code " << fPDGType << " not defined.";
    fPDGMass = particle->Mass();
    return kTRUE;
}

Bool_t R3BCALIFATestGenerator::ReadEvent(FairPrimaryGenerator* primGen)
{
    // Generate one event: produce primary particles emitted from one vertex.
    // Primary particles are distributed uniformly along
    // those kinematics variables which were limitted by setters.
    // if SetCosTheta() function is used, the distribution will be uniform in
    // cos(theta)

    Double32_t pabs = 0, phi, pt = 0, theta = 0, eta, y, mt, px, py, pz = 0;
    Double32_t br = 0;
    Bool_t doNotBoost = 0;

    // Generate particles
    for (Int_t k = 0; k < fMult; k++)
    {
        phi = gRandom->Uniform(fPhiMin, fPhiMax) * TMath::DegToRad();

        if (fPRangeIsSet)
            pabs = gRandom->Uniform(fPMin, fPMax);
        else if (fPtRangeIsSet)
            pt = gRandom->Uniform(fPtMin, fPtMax);

        if (fThetaRangeIsSet)
        {
            if (fCosThetaIsSet)
                theta = acos(gRandom->Uniform(cos(fThetaMin * TMath::DegToRad()), cos(fThetaMax * TMath::DegToRad())));
            else
                theta = gRandom->Uniform(fThetaMin, fThetaMax) * TMath::DegToRad();
        }
        else if (fEtaRangeIsSet)
        {
            eta = gRandom->Uniform(fEtaMin, fEtaMax);
            theta = 2 * TMath::ATan(TMath::Exp(-eta));
        }
        else if (fYRangeIsSet)
        {
            y = gRandom->Uniform(fYMin, fYMax);
            mt = TMath::Sqrt(fPDGMass * fPDGMass + pt * pt);
            pz = mt * TMath::SinH(y);
        }

        if (fThetaRangeIsSet || fEtaRangeIsSet)
        {
            if (fPRangeIsSet)
            {
                pz = pabs * TMath::Cos(theta);
                pt = pabs * TMath::Sin(theta);
            }
            else if (fPtRangeIsSet)
                pz = pt / TMath::Tan(theta);
        }

        px = pt * TMath::Cos(phi);
        py = pt * TMath::Sin(phi);

        if (fBoxVtxIsSet)
        {
            fX = gRandom->Uniform(fX1, fX2);
            fY = gRandom->Uniform(fY1, fY2);
            fZ = gRandom->Uniform(fZ1, fZ2);
        }

        if (fNuclearDecayChainIsSet)
        {
            if (fPDGType != 22)
                LOG(fatal) << "R3BCALIFATestGenerator: PDG code " << fPDGType << " is not a gamma!";
            br = gRandom->Uniform();
            for (Int_t i = 0; i < fGammasDefinedInNuclearDecay; i++)
            {
                if (br < fGammaBranchingRatios[i])
                {
                    Double32_t gammaMomentum = TMath::Sqrt(px * px + py * py + pz * pz);
                    px = px * fGammaEnergies[i] / gammaMomentum;
                    py = py * fGammaEnergies[i] / gammaMomentum;
                    pz = pz * fGammaEnergies[i] / gammaMomentum;
                    break;
                }
            }
            // if Sum(branchingRatios)<1, the leftover probability (up to 1) is defined as environmental noise
            doNotBoost = 1;
        }
        /*
      if (fLorentzBoostIsSet && !doNotBoost){

        //Lorentz transformation Pz(lab) = gamma * Pz(cm) + gamma * beta * E
        //As each Lorentz transformation can be performed sequencially,
        //we can separate the gamma factor corresponding to each direction
        Double32_t gammaMomentum=TMath::Sqrt(px*px+py*py+pz*pz);
        pz = (pz + fBetaOfEmittingFragment * gammaMomentum) / fGammaFactor;
          */

        if (fPDGType == 22 && fLorentzBoostIsSet && !doNotBoost)
        { /// for gamma-rays
            // Lorentz transformation Pz(lab) = gamma * Pz(cm) + gamma * beta * E
            // As each Lorentz transformation can be performed sequencially,
            // we can separate the gamma factor corresponding to each direction
            Double32_t gammaMomentum = TMath::Sqrt(px * px + py * py + pz * pz);
            pz = (pz + fBetaOfEmittingFragment * gammaMomentum) / fGammaFactor;
        }
        else if (fLorentzBoostIsSet && !doNotBoost)
        { /// for any massive particle
            // Lorentz transformation Pz(lab) = gamma * Pz(cm) + gamma * beta * E
            // As each Lorentz transformation can be performed sequencially,
            // we can separate the gamma factor corresponding to each direction
            Double32_t particleEnergy = TMath::Sqrt(px * px + py * py + pz * pz + fPDGMass * fPDGMass);
            pz = (pz + fBetaOfEmittingFragment * particleEnergy) / fGammaFactor;
        }

        if (fDebug)
            printf("CALIFATestGen: kf=%d, p=(%.2f, %.2f, %.2f) GeV, x=(%.1f, %.1f, %.1f) cm\n",
                   fPDGType,
                   px,
                   py,
                   pz,
                   fX,
                   fY,
                   fZ);

        primGen->AddTrack(fPDGType, px, py, pz, fX, fY, fZ);
    }
    return kTRUE;
}

void R3BCALIFATestGenerator::SetFragmentVelocity(Double32_t beta)
{
    // Sets the velocity and gamma factor of the fragment emitting the gammas
    fBetaOfEmittingFragment = beta;
    fGammaFactor = TMath::Sqrt(1 - fBetaOfEmittingFragment * fBetaOfEmittingFragment);
}

void R3BCALIFATestGenerator::SetDecayChainPoint(Double32_t gammaEnergy, Double32_t branchingRatio)
{
    //
    //
    //
    if (fGammasDefinedInNuclearDecay > 7)
        printf("CALIFATestGen: Maximum number (8) of gammas defined in the chain\n");
    else
    {
        fGammaEnergies[fGammasDefinedInNuclearDecay] = gammaEnergy;
        fGammaBranchingRatios[fGammasDefinedInNuclearDecay] = branchingRatio;
        fGammasDefinedInNuclearDecay++;
    }
}

ClassImp(R3BCALIFATestGenerator)
