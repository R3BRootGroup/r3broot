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
// -----                      R3BVetoPoint header file                  -----
// -----                  Created 26/03/09  by D.Bertini               -----
// -------------------------------------------------------------------------

/**  R3BVetoPoint.h
 **/

#ifndef R3BVETOPOINT_H
#define R3BVETOPOINT_H

#include "TObject.h"
#include "TVector3.h"

#include "FairMCPoint.h"

class R3BVetoPoint : public FairMCPoint
{

  public:
    /** Default constructor **/
    R3BVetoPoint();

    /** Constructor with arguments
     *@param trackID  Index of MCTrack
     *@param detID    Detector ID
     *@param posIn    Ccoordinates at entrance to active volume [cm]
     *@param posOut   Coordinates at exit of active volume [cm]
     *@param momIn    Momentum of track at entrance [GeV]
     *@param momOut   Momentum of track at exit [GeV]
     *@param tof      Time since event start [ns]
     *@param length   Track length since creation [cm]
     *@param eLoss    Energy deposit [GeV]
     **/
    R3BVetoPoint(Int_t trackID,
                 Int_t detID,
                 Int_t box,
                 Int_t sector,
                 Int_t paddle,
                 TVector3 posIn,
                 TVector3 posOut,
                 TVector3 momIn,
                 TVector3 momOut,
                 Double_t tof,
                 Double_t length,
                 Double_t eLoss);

    /** Copy constructor **/
    R3BVetoPoint(const R3BVetoPoint& point) { *this = point; };

    /** Destructor **/
    virtual ~R3BVetoPoint();

    /** Accessors **/
    Double_t GetXIn() const { return fX; }
    Double_t GetYIn() const { return fY; }
    Double_t GetZIn() const { return fZ; }
    Double_t GetXOut() const { return fX_out; }
    Double_t GetYOut() const { return fY_out; }
    Double_t GetZOut() const { return fZ_out; }
    Double_t GetPxOut() const { return fPx_out; }
    Double_t GetPyOut() const { return fPy_out; }
    Double_t GetPzOut() const { return fPz_out; }
    Int_t GetSector() const { return fSector; }
    Int_t GetPaddleType() const { return fPaddleTyp; }
    Int_t GetPaddleNb() const { return fPaddleNb; }
    void PositionIn(TVector3& pos) { pos.SetXYZ(fX, fY, fZ); }
    void PositionOut(TVector3& pos) { pos.SetXYZ(fX_out, fY_out, fZ_out); }
    void MomentumOut(TVector3& mom) { mom.SetXYZ(fPx_out, fPy_out, fPz_out); }

    /** Point coordinates at given z from linear extrapolation **/
    Double_t GetX(Double_t z) const;
    Double_t GetY(Double_t z) const;

    /** Check for distance between in and out **/
    Bool_t IsUsable() const;

    /** Modifiers **/
    void SetPositionOut(TVector3 pos);
    void SetMomentumOut(TVector3 mom);

    /** Output to screen **/
    virtual void Print(const Option_t* opt) const;

  protected:
    Double32_t fX_out, fY_out, fZ_out;
    Double32_t fPx_out, fPy_out, fPz_out;
    Int_t fSector;
    Int_t fPaddleNb;
    Int_t fPaddleTyp;

    ClassDef(R3BVetoPoint, 1)
};

inline void R3BVetoPoint::SetPositionOut(TVector3 pos)
{
    fX_out = pos.X();
    fY_out = pos.Y();
    fZ_out = pos.Z();
}

inline void R3BVetoPoint::SetMomentumOut(TVector3 mom)
{
    fPx_out = mom.Px();
    fPy_out = mom.Py();
    fPz_out = mom.Pz();
}

#endif
