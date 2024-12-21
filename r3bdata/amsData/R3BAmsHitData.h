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
// -----              R3BAmsHitData header file                        -----
// -----             Created 01/06/18  by J.L. Rodriguez-Sanchez       -----
// -------------------------------------------------------------------------

/**  R3BAmsHitData.h
 **  A HitData is the representation of the
 **  information obtained from each ams detector.
 **/

#ifndef R3BAmsHitData_H
#define R3BAmsHitData_H

#include "FairMultiLinkedData.h"
#include "TObject.h"
#include "TVector3.h"

class R3BAmsHitData : public FairMultiLinkedData
{
  public:
    /** Default constructor **/
    R3BAmsHitData();

    /** Constructor with arguments
     *@param fDetId     Detector unique identifier
     *@param fNumHit    Number of hits
     *@param fX         Position X (S-Side) [mm] local coordinate (detector frame)
     *@param fY         Position Y (K-Side) [mm] local coordinate (detector frame)
     *@param fTheta     Angle theta [rad] (lab frame)
     *@param fPhi       Angle Phi [rad] (lab frame)
     *@param fEnergyS   Total energy deposited on S-side by hit ([GeV] in sim)
     *@param fEnergyK   Total energy deposited on K-side by hit ([GeV] in sim)
     **/
    R3BAmsHitData(Int_t detid,
                  Int_t numhit,
                  Double_t x,
                  Double_t y,
                  TVector3 master,
                  Double_t energy_x,
                  Double_t energy_y,
                  Int_t mulS = 0,
                  Int_t mulK = 0);

    /** Copy constructor **/
    R3BAmsHitData(const R3BAmsHitData&);

    R3BAmsHitData& operator=(const R3BAmsHitData&) { return *this; }

    /** Destructor **/
    virtual ~R3BAmsHitData() {}

    /** Accessors **/
    inline const Int_t& GetDetId() const { return fDetId; }
    inline const Int_t& GetNumHit() const { return fNumHit; }
    inline const Int_t& GetMulS() const { return fMulS; }
    inline const Int_t& GetMulK() const { return fMulK; }
    inline const Double_t& GetX() const { return fX; }
    inline const Double_t& GetY() const { return fY; }
    inline const Double_t& GetPos_S() const { return fX; }
    inline const Double_t& GetPos_K() const { return fY; }
    inline const Double_t& GetTheta() const { return fTheta; }
    inline const Double_t& GetPhi() const { return fPhi; }
    inline const TVector3 GetPosLab() const { return fmaster; }
    inline const Double_t& GetEnergyS() const { return fEnergyS; }
    inline const Double_t& GetEnergyK() const { return fEnergyK; }

  protected:
    Int_t fDetId;
    Int_t fNumHit;
    Int_t fMulS, fMulK;
    Double_t fX, fY, fTheta, fPhi;
    TVector3 fmaster;
    Double_t fEnergyS, fEnergyK;

    ClassDef(R3BAmsHitData, 1)
};

#endif
