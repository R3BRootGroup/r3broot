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
// -----                   R3BIonGenerator header file                 -----
// -----                    Based on FairIonGenerator                  -----
// -----                 Created 18/03/11  by M. Labiche               -----
// -------------------------------------------------------------------------

/** Include beam spot size **/

#pragma once

#include "R3BBeamProperties.h"

#include <FairGenerator.h>
#include <FairIon.h>

#include <TRandom3.h>

#include <iostream>

class FairPrimaryGenerator;

class R3BIonGenerator : public FairGenerator
{

  public:
    /** Default constructor **/
    R3BIonGenerator(UInt_t seed = 0U);

    /** Constructor with ion name
     ** For the generation of ions with pre-defined FairIon
     ** By default, the  excitation energy is zero. This can be changed with the
     ** respective modifiers.
     **@param ionName  Ion name
     **@param mult      Number of ions per event
     **@param px,py,pz  Momentum components [GeV] per nucleon!
     **@param vx,vy,vz  Vertex coordinates [cm]
     **/
    R3BIonGenerator(const Char_t* ionName, Int_t mult, Double_t momentum_AGeV_per_c, UInt_t seed = 0U);

    /** Default constructor
     ** For the generation of ions with atomic number z and mass number a.
     ** By default, the mass equals a times the proton mass and the
     ** excitation energy is zero. This can be changed with the
     ** respective modifiers.
     **@param z         Atomic number
     **@param a         Atomic mass
     **@param q         Electric charge [e]
     **@param mult      Number of ions per event
     **@param px,py,pz  Momentum components [GeV] per nucleon!
     **@param vx,vy,vz  Vertex coordinates [cm]
     **/
    R3BIonGenerator(Int_t z, Int_t a, Int_t q, Int_t mult, Double_t momentum_AGeV_per_c, UInt_t seed = 0U);

    R3BIonGenerator(const R3BIonGenerator&);

    R3BIonGenerator& operator=(const R3BIonGenerator&) { return *this; }

    /** Destructor **/
    virtual ~R3BIonGenerator();

    /** Modifiers **/
    void SetExcitationEnergy(Double_t eExc);
    void SetMass(Double_t mass);

    Bool_t Init() override;

    /** Method ReadEvent
     ** Generates <mult> of the specified ions and hands hem to the
     ** FairPrimaryGenerator.
     **/
    virtual Bool_t ReadEvent(FairPrimaryGenerator* primGen) override;

    R3BBeamProperties Beam; //!

  private:
    Int_t fMult;   // Multiplicity per event
    FairIon* fIon; // Pointer to the FairIon to be generated

    TRandom3 fRNG; // the RNG
    ClassDefOverride(R3BIonGenerator, 1)
};
