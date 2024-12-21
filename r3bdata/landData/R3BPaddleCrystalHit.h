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
// -----                      R3BLandPoint header file                  -----
// -----                  Created 26/03/09  by D.Bertini               -----
// -------------------------------------------------------------------------
/**  R3BPaddleCrystalHit.h
 **/

#ifndef R3BPADDLECRYSTALHIT_H
#define R3BPADDLECRYSTALHIT_H

#include "FairHit.h"
#include "TObject.h"
#include "TVector3.h"

class R3BPaddleCrystalHit : public FairHit
{

  public:
    /** Default constructor **/
    R3BPaddleCrystalHit();
    /** Copy constructor **/
    R3BPaddleCrystalHit(const R3BPaddleCrystalHit&);
    /** Destructor **/
    virtual ~R3BPaddleCrystalHit();

    R3BPaddleCrystalHit& operator=(const R3BPaddleCrystalHit&) { return *this; }

    /** Output to screen **/
    virtual void Print(const Option_t* opt) const;

    void SetType(UInt_t ty) { type = ty; }
    void SetCopy(UInt_t co) { copy = co; }
    void SetEnergy(Double_t en) { fEnergy = en; }
    void SetNbOfSteps(UInt_t nb) { nbOfSteps = nb; }

    void SetTimeFirstStep(Double_t ti) { fTimeFirstStep = ti; }
    void SetTimeLastStep(Double_t ti) { fTimeLastStep = ti; }
    void SetPositionX(Double_t posX) { fX = posX; }
    void SetPositionY(Double_t posY) { fY = posY; }
    void SetPositionZ(Double_t posZ) { fZ = posZ; }

    void SetNbOfPrimaries(UInt_t nb) { nbOfPrimaries = nb; }
    void SetEnergyPrimary(Double_t ene) { fEnergyPrimary = ene; }
    void SetThetaPrimary(Double_t the) { fThetaPrimary = the; }
    void SetPhiPrimary(Double_t phi) { fPhiPrimary = phi; }

    void SetEventID(UInt_t ev) { eventID = ev; }
    void SetRunID(UInt_t run) { runID = run; }

    Int_t GetType() { return type; }
    Int_t GetCopy() { return copy; }
    Double_t GetEnergy() { return fEnergy; }
    UInt_t GetNbOfSteps() { return nbOfSteps; }

    Double_t GetTimeFirstStep() { return fTimeFirstStep; }
    Double_t GetTimeLastStep() { return fTimeLastStep; }

    UInt_t GetNbOfPrimaries() { return nbOfPrimaries; }
    Double_t GetEnergyPrimary() { return fEnergyPrimary; }
    Double_t GetThetaPrimary() { return fThetaPrimary; }
    Double_t GetPhiPrimary() { return fPhiPrimary; }

    UInt_t GetEventID() { return eventID; }
    UInt_t GetRunID() { return runID; }

  protected:
    // Crystal identification
    UInt_t type; // crystal type
    UInt_t copy; // crystal copy

    // Basic Hit information
    Double32_t fEnergy; // total energy deposited on the paddle
    UInt_t nbOfSteps;   // number of step depositing energy in the crystal

    Double32_t fTimeFirstStep; // time of the first interaction
    Double32_t fTimeLastStep;  // time of the last interaction

    // Primary identification
    UInt_t nbOfPrimaries; // for a gamma source
    Double32_t fEnergyPrimary;
    Double32_t fThetaPrimary;
    Double32_t fPhiPrimary;

    // Event identification
    UInt_t eventID; // event ID
    UInt_t runID;   // run ID

    ClassDef(R3BPaddleCrystalHit, 1)
};
#endif
