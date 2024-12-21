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

/**  R3BCALIFATestGenerator.h
 *@author H. Alvarez Pol <hector.alvarez@usc.es>
 *
 * The R3BCALIFATestGenerator generates gammas with different options
 * for testing CALIFA. Copies from FairBoxGenerator (I tried to derive from it,
 * but requires a modification of their data members to protected).
 * Derived from FairGenerator.
 **/

#ifndef FAIR_CALIFATESTGENERATOR_H
#define FAIR_CALIFATESTGENERATOR_H

#include "FairGenerator.h"

#include <iostream>

class FairPrimaryGenerator;

class R3BCALIFATestGenerator : public FairGenerator
{
  public:
    /** Default constructor. **/
    R3BCALIFATestGenerator();

    /** Constructor with PDG-ID, multiplicity
     **@param pdgid Particle type (PDG encoding)
     **@param mult  Multiplicity (default is 1)
     **/
    R3BCALIFATestGenerator(Int_t pdgid, Int_t mult = 1);

    /** Destructor **/
    virtual ~R3BCALIFATestGenerator(){};

    /** Modifiers **/
    void SetPDGType(Int_t pdg) { fPDGType = pdg; };

    void SetMultiplicity(Int_t mult) { fMult = mult; };

    void SetPRange(Double32_t pmin = 0, Double32_t pmax = 10)
    {
        fPMin = pmin;
        fPMax = pmax;
        fPRangeIsSet = kTRUE;
    }

    void SetPtRange(Double32_t ptmin = 0, Double32_t ptmax = 10)
    {
        fPtMin = ptmin;
        fPtMax = ptmax;
        fPtRangeIsSet = kTRUE;
    };

    void SetPhiRange(Double32_t phimin = 0, Double32_t phimax = 360)
    {
        fPhiMin = phimin;
        fPhiMax = phimax;
    };

    void SetEtaRange(Double32_t etamin = -5, Double32_t etamax = 7)
    {
        fEtaMin = etamin;
        fEtaMax = etamax;
        fEtaRangeIsSet = kTRUE;
    };

    void SetYRange(Double32_t ymin = -5, Double32_t ymax = 7)
    {
        fYMin = ymin;
        fYMax = ymax;
        fYRangeIsSet = kTRUE;
    };

    void SetThetaRange(Double32_t thetamin = 0, Double32_t thetamax = 90)
    {
        fThetaMin = thetamin;
        fThetaMax = thetamax;
        fThetaRangeIsSet = kTRUE;
    };

    void SetCosTheta() { fCosThetaIsSet = kTRUE; };

    void SetXYZ(Double32_t x = 0, Double32_t y = 0, Double32_t z = 0)
    {
        fX = x;
        fY = y;
        fZ = z;
        fPointVtxIsSet = kTRUE;
    }

    void SetBoxXYZ(Double32_t x1 = 0,
                   Double32_t y1 = 0,
                   Double32_t z1 = 0,
                   Double32_t x2 = 0,
                   Double32_t y2 = 0,
                   Double32_t z2 = 0)
    {
        fX1 = x1;
        fY1 = y1;
        fZ1 = z1;
        fX2 = x2;
        fY2 = y2;
        fZ2 = z2;
        fBoxVtxIsSet = kTRUE;
    }

    void SetDebug(Bool_t debug = 0) { fDebug = debug; }

    void SetLorentzBoost(Double32_t beta = 0)
    {
        SetFragmentVelocity(beta);
        fLorentzBoostIsSet = kTRUE;
    };

    void SetFragmentVelocity(Double32_t beta = 0);

    void SetNuclearDecayChain() { fNuclearDecayChainIsSet = kTRUE; };

    void SetDecayChainPoint(Double32_t gammaEnergy = 0, Double32_t branchingRatio = 0);

    /** Initializer **/
    virtual Bool_t Init();

    /** Creates an event with given type and multiplicity.
     **@param primGen  pointer to the FairPrimaryGenerator
     **/
    virtual Bool_t ReadEvent(FairPrimaryGenerator* primGen);

  private:
    Int_t fPDGType; // Particle type (PDG encoding)
    Int_t fMult;    // Multiplicity

    Double32_t fPDGMass;                     // Particle mass [GeV]
    Double32_t fPtMin, fPtMax;               // Transverse momentum range [GeV]
    Double32_t fPhiMin, fPhiMax;             // Azimuth angle range [degree]
    Double32_t fEtaMin, fEtaMax;             // Pseudorapidity range in lab system
    Double32_t fYMin, fYMax;                 // Rapidity range in lab system
    Double32_t fPMin, fPMax;                 // Momentum range in lab system
    Double32_t fThetaMin, fThetaMax;         // Polar angle range in lab system [degree]
    Double32_t fX, fY, fZ;                   // Point vertex coordinates [cm]
    Double32_t fX1, fY1, fZ1, fX2, fY2, fZ2; // Box vertex coords (x1,y1,z1)->(x2,y2,z2)

    Bool_t fEtaRangeIsSet;   // True if eta range is set
    Bool_t fYRangeIsSet;     // True if rapidity range is set
    Bool_t fThetaRangeIsSet; // True if theta range is set
    Bool_t fCosThetaIsSet;   // True if uniform distribution in
                             // cos(theta) is set (default -> not set)
    Bool_t fPtRangeIsSet;    // True if transverse momentum range is set
    Bool_t fPRangeIsSet;     // True if abs.momentum range is set
    Bool_t fPointVtxIsSet;   // True if point vertex is set
    Bool_t fBoxVtxIsSet;     // True if box vertex is set
    Bool_t fDebug;           // Debug switch

    // SPECIFIC OF CALIFA
    Double32_t fGammaEnergies[8];        // Gamma energies for the nuclear decay chain [GeV] TEST UNITS!
    Double32_t fGammaBranchingRatios[8]; // Gamma branching ratios for the nuclear decay chain
    Int_t fGammasDefinedInNuclearDecay;  // Number of gammas defined (max 8)
    Double32_t fBetaOfEmittingFragment;  // Velocity of the fragment emitting the gammas
    Double32_t fGammaFactor;             // Velocity of the fragment emitting the gammas

    Bool_t fLorentzBoostIsSet;      // True if Lorentz Boost is set
    Bool_t fNuclearDecayChainIsSet; // True if a nuclear decay chain is set

    ClassDef(R3BCALIFATestGenerator, 1);
};

#endif
