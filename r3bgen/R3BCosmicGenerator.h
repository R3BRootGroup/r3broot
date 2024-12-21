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
// -----                R3BCosmicGenerator header file                 -----
// -----               Created 09/09/04  by Marc Labiche               -----
// -----            Derived from R3BBoxGenerator header file           -----
// -------------------------------------------------------------------------

/**  R3BCosmicGenerator.h
 *@author Marc Labiche <marc.labiche@stfc.ac.uk>
 *
 * The R3BBoxGenerator generates particles within given range
 * over phi,theta and a fixed multiplicity per event.
 * Derived from R3BGenerator.
 **/

#ifndef R3B_COSMICGENERATOR_H
#define R3B_COSMICGENERATOR_H

#include "R3BDistribution.h"

#include "FairGenerator.h"

#include "TRandom3.h"

#include <iostream>

class FairPrimaryGenerator;

class R3BCosmicGenerator : public FairGenerator
{
  public:
    R3BCosmicGenerator(const Int_t pdgid = 13, const Int_t mult = 1, const UInt_t seed = 0U);

    /** Destructor **/
    virtual ~R3BCosmicGenerator(){};

    /** Modifiers **/
    void SetPDGType(const Int_t pdg) { fPDGType = pdg; };
    void SetMultiplicity(const Int_t mult) { fMult = mult; };
    void SetVertexDistribution_cm(const R3BDistribution<3> vertexDist_cm) { fVertexPositionDist_cm = vertexDist_cm; }
    void SetAngularDistribution_Rad(const R3BDistribution<2> angularDist_Rad) { fAngularDist_Rad = angularDist_Rad; }
    void SetEnergyRange_GeV(const Double_t minEnergy = 0.004, const Double_t maxEnergy = 100.);

    R3BDistribution<3>& GetVertexDistribution_cm() { return fVertexPositionDist_cm; }
    R3BDistribution<2>& GetAngularDistribution_Rad() { return fAngularDist_Rad; }
    R3BDistribution<1>& GetEnergyDistribution_GeV() { return fEnergyDist_GeV; };

    /** Initializer **/
    Bool_t Init() override;

    /** Creates an event with given type and multiplicity.
     **@param primGen  pointer to the R3BPrimaryGenerator
     **/
    Bool_t ReadEvent(FairPrimaryGenerator* primGen) override;

  private:
    R3BDistribution<3> fVertexPositionDist_cm; //!
    R3BDistribution<2> fAngularDist_Rad;       //!
    R3BDistribution<1> fEnergyDist_GeV;        //!

    TRandom3 fRngGen;

    Int_t fPDGType; // Particle type (PDG encoding)
    Int_t fMult;    // Multiplicity

    Double_t fPDGMass; // Particle mass [GeV]

    ClassDefOverride(R3BCosmicGenerator, 2)
};

#endif
