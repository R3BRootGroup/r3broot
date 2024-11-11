/******************************************************************************
 *   Copyright (C) 2019 GSI Helmholtzzentrum für Schwerionenforschung GmbH    *
 *   Copyright (C) 2019 Members of R3B Collaboration                          *
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

#ifndef R3BIONGENERATOR_H
#define R3BIONGENERATOR_H

#include "FairGenerator.h"
#include "FairIon.h"

#include <iostream>

class FairPrimaryGenerator;

class R3BIonGenerator : public FairGenerator
{

  public:
    /** Default constructor **/
    R3BIonGenerator();

    /** Constructor with ion name
     ** For the generation of ions with pre-defined FairIon
     ** By default, the  excitation energy is zero. This can be changed with the
     ** respective modifiers.
     **@param ionName  Ion name
     **@param mult      Number of ions per event
     **@param px,py,pz  Momentum components [GeV] per nucleon!
     **@param vx,vy,vz  Vertex coordinates [cm]
     **/
    R3BIonGenerator(const Char_t* ionName, Int_t mult, Double_t px, Double_t py, Double_t pz);

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
    R3BIonGenerator(Int_t z, Int_t a, Int_t q, Int_t mult, Double_t px, Double_t py, Double_t pz);

    R3BIonGenerator(const R3BIonGenerator&);

    R3BIonGenerator& operator=(const R3BIonGenerator&) { return *this; }

    /** Destructor **/
    virtual ~R3BIonGenerator();

    /** Modifiers **/
    void SetCharge(Int_t charge) { fQ = charge; }
    void SetExcitationEnergy(Double_t eExc);
    void SetMass(Double_t mass);

    void SetSpotRadius(Double32_t r = 0, Double32_t z = 0, Double32_t offX = 0, Double32_t offY = 0)
    {
        fR = r;
        fz = z;
        fOffsetX = offX;
        fOffsetY = offY;
        fBeamSpotIsSet = kTRUE;
    }

    void SetBeamParameter(Double32_t sigmaP = 0, Double32_t angleX = 0, Double32_t angleY = 0)
    {
        fSigmaP = sigmaP;
        fAngleX = angleX;
        fAngleY = angleY;
    }

    /** Method ReadEvent
     ** Generates <mult> of the specified ions and hands hem to the
     ** FairPrimaryGenerator.
     **/
    virtual Bool_t ReadEvent(FairPrimaryGenerator* primGen);

  private:
    static Int_t fgNIon;        //! Number of the instance of this class
    Int_t fMult;                // Multiplicity per event
    Double_t fPx, fPy, fPz;     // Momentum components [GeV] per nucleon
    Double32_t fR, fz, fOffsetX, fOffsetY; // beam Spot radius [cm], z source, x, y source
    Double_t fVx, fVy, fVz;     // Vertex coordinates [cm]
    FairIon* fIon;              // Pointer to the FairIon to be generated
    Int_t fQ;                   // Electric charge [e]
    Double32_t fSigmaP, fAngleX, fAngleY; // beam variation in momentum, angle of emittance

    Bool_t fBeamSpotIsSet; // True if point beamspot is set

    ClassDef(R3BIonGenerator, 1)
};

#endif
