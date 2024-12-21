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

/////////////////////////////////////////////////////////////////
//*-- AUTHOR : Marc Labiche		labi-ph0@paisley.ac.uk
//	       Oleg Kiselev		O.Kiselev@gsi.de
//             Hector Alvarez-Pol       hapol@fpddux.usc.es
//             Saúl Beceiro             saul.beceiro@rai.usc.es
//             Denis Bertini          <D.Bertini@gsi.de>
//*-- Date: 11/2005
//*-- Last Update: 3/08/09  <D.Bertini@gsi.de>
// --------------------------------------------------------------
// Description:
//   Event generator interface for R3B
//
// --------------------------------------------------------------
// Comments:
//               3/08/09 - adapted for R3BRoot
//                         <D.Bertini@gsi.de>
// --------------------------------------------------------------
/////////////////////////////////////////////////////////////////

#ifndef R3BSpecificGenerator_H
#define R3BSpecificGenerator_H 1

#include "TString.h"
#include "TVector3.h"

#include "TDatabasePDG.h"
#include "TMath.h"
#include "TParticlePDG.h"
#include "TRandom.h"

#include "FairGenerator.h"

#include "R3BBackTracking.h"
#include "R3BCDGenerator.h"
#include "R3BReadKinematics.h"

class TVector3;
class TString;

class R3BSpecificGenerator : public FairGenerator
{

  private:
    R3BReadKinematics* pReadKinematics;
    R3BCDGenerator* pCDGenerator;
    R3BBackTracking* pBackTrackingGenerator;

    TString gammasFlag;       // flag to select an isotropic gamma emitter
    TString decaySchemeFlag;  //
    TString reactionFlag;     // flag to select a reaction
    TString reactionType;     // reaction type
    TString dissociationFlag; // Coulomb dissociation generator flag
    TString backTrackingFlag; // BackTracking generator flag

    TString targetType; // target elements
    Double_t targetHalfThicknessPara;
    Double_t targetThicknessLiH;
    Double_t targetRadius;

    TString beamInteractionFlag; // flag to select target/beam interaction-like position for gammas origin
    TString rndmFlag;            // flag for a rndm (angle) primary emission
    TString rndmEneFlag;         // flag for a rndm energy primary emission
    TString boostFlag;           // flag for a boosted primary emission
    Int_t fPDGType;              // Particle type (PDG encoding)
    Double_t kinEnergyPrim;      // kinetic energy of the primary
    Double_t meanKinEnergyBeam;  // kinetic energy mean of the beam (per nucleon)
    Double_t sigmaKinEnergyBeam; // kinetic energy sigma of the beam

    TString simEmittanceFlag; // flag to select an emittance spectra for vertex emission in CoulDiss
    Double_t sigmaXInEmittance;
    Double_t sigmaXPrimeInEmittance;

    Double_t fPDGMass; // Particle Mass from internal PDG table

    Int_t fMult;    // Multiplicity
    Double_t fP;    // Momentum
    TVector3 fPdir; // Momentum Direction
    Int_t fCharge;  // Charge
    TVector3 fPol;  // Polarization
    TVector3 fPos;  // Vertex position
    Double_t fTime; // Decay Time
    TString particlePrim;
    Bool_t isDumped;

  public:
    R3BSpecificGenerator();
    R3BSpecificGenerator(Int_t pdg, Float_t beamEnergy);
    R3BSpecificGenerator(const R3BSpecificGenerator&);
    ~R3BSpecificGenerator();

    R3BSpecificGenerator& operator=(const R3BSpecificGenerator&) { return *this; }

    /** Initializer **/
    virtual Bool_t Init();

    /** Creates an event with given type and multiplicity.
    **@param primGen  pointer to the FairSpecificGenerator
    **/
    virtual Bool_t ReadEvent(FairPrimaryGenerator* primGen);

    void SetBeamInteractionFlag(TString val) { beamInteractionFlag = val; }
    void SetRndmFlag(TString val) { rndmFlag = val; }
    void SetRndmEneFlag(TString val) { rndmEneFlag = val; }
    void SetBoostFlag(TString val) { boostFlag = val; }
    void SetBeamEnergy(Double_t val) { meanKinEnergyBeam = val; }
    void SetSigmaBeamEnergy(Double_t val) { sigmaKinEnergyBeam = val; }

    void SetParticlePrim(TString val) { particlePrim = val; }

    void SetEnergyPrim(Double_t val)
    {
        if (reactionFlag.CompareTo("on") == 0)
        {
            cout << "-I- R3BSpecificGenerator::SetEnergyPrim()  \
                 FLAG(Reaction) is on "
                 << endl;
            cout << "-I  R3BSpecificGenerator::SetEnergyPrim()  \
                 Sampling between [50.,400.] MeV activated ..."
                 << endl;
        }

        kinEnergyPrim = val;
    }

    void SetTargetType(TString ans) { targetType = ans; }
    void SetReactionFlag(TString val) { reactionFlag = val; }
    void SetGammasFlag(TString val) { gammasFlag = val; }
    void SetDecaySchemeFlag(TString val) { decaySchemeFlag = val; }
    void SetReactionType(TString val) { reactionType = val; }
    void SetTargetHalfThicknessPara(Double_t para) { targetHalfThicknessPara = para; }
    void SetTargetThicknessLiH(Double_t para) { targetThicknessLiH = para; }
    void SetTargetRadius(Double_t para) { targetRadius = para; }

    void SetDissociationFlag(TString val) { dissociationFlag = val; }
    void SetBackTrackingFlag(TString val) { backTrackingFlag = val; }
    void SetSimEmittanceFlag(TString val) { simEmittanceFlag = val; }
    void SetSigmaXInEmittance(Double_t val) { sigmaXInEmittance = val; }
    void SetSigmaXPrimeInEmittance(Double_t val) { sigmaXPrimeInEmittance = val; }

    //
    TString GetBeamInteractionFlag() { return beamInteractionFlag; }
    TString GetRndmFlag() { return rndmFlag; }
    TString GetRndmEneFlag() { return rndmEneFlag; }
    TString GetBoostFlag() { return boostFlag; }
    Double_t GetBeamEnergy() { return meanKinEnergyBeam; }
    TString GetParticlePrim() { return particlePrim; }
    Double_t GetEnergyPrim() { return kinEnergyPrim; }
    TString GetTargetType() { return targetType; }
    TString GetReactionFlag() { return reactionFlag; }
    TString GetGammasFlag() { return gammasFlag; }
    TString GetDecaySchemeFlag() { return decaySchemeFlag; }
    TString GetReactionType() { return reactionType; }

    Double_t GetTargetHalfThicknessPara() { return targetHalfThicknessPara; }
    Double_t GetTargetThicknessLiH() { return targetThicknessLiH; }
    Double_t GetTargetRadius() { return targetRadius; }

    TString GetDissociationFlag() { return dissociationFlag; }
    TString GetBackTrackingFlag() { return backTrackingFlag; }
    TString GetSimEmittanceFlag() { return simEmittanceFlag; }
    Double_t GetSigmaXInEmittance() { return sigmaXInEmittance; }
    Double_t GetSigmaXPrimeInEmittance() { return sigmaXPrimeInEmittance; }

    void SetParticleDefinition(Int_t pdg_id) { fPDGType = pdg_id; }

    inline void SetParticleMomentum(Double32_t aMomentum) { fP = aMomentum; }

    inline void SetParticleMomentumDirection(TVector3 aMomentumDirection) { fPdir = aMomentumDirection; }

    inline void SetParticleCharge(Int_t aCharge) { fCharge = aCharge; }

    inline void SetParticlePolarization(TVector3 aVal) { fPol = aVal; }

    inline void SetParticlePosition(TVector3 aPos) { fPos = aPos; }

    inline void SetParticleTime(Double32_t aTime) { fTime = aTime; }

    inline void SetNumberOfParticles(Int_t i) { fMult = i; }

    inline Int_t GetPrimPDGCode() { return fPDGType; }

    inline TParticlePDG* GetParticleDefinition()
    {
        TDatabasePDG* pdgBase = TDatabasePDG::Instance();
        TParticlePDG* particle = pdgBase->GetParticle(fPDGType);

        return particle;
    }

    inline TVector3 GetParticleMomentumDirection() { return TVector3(fP, 0, 0); }

    inline Int_t GetParticleCharge() { return fCharge; }

    inline TVector3 GetParticlePolarization() { return fPol; }

    inline TVector3 GetParticlePosition() { return fPos; }

    inline Double_t GetParticleTime() { return fTime; }

    inline Int_t GetNumberOfParticles() { return fMult; }

    void PrintParameters();

    ClassDef(R3BSpecificGenerator, 1);
};

#endif
