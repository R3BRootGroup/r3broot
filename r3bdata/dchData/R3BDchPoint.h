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
// -----                      R3BDchPoint header file                  -----
// -----                  Created 26/03/09  by D.Bertini               -----
// -------------------------------------------------------------------------

/**  R3BDchPoint.h
 **/

#ifndef R3BDCHPOINT_H
#define R3BDCHPOINT_H

#include "TObject.h"
#include "TVector3.h"

#include "FairMCPoint.h"

class R3BDchPoint : public FairMCPoint
{

  public:
    /** Default constructor **/
    R3BDchPoint();

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
    R3BDchPoint(Int_t trackId,
                Int_t mod,
                Int_t layer,
                Int_t cell,
                TVector3 posIn,
                TVector3 posOut,
                TVector3 momIn,
                TVector3 momOut,
                TVector3 lpos1,
                TVector3 lmom1,
                TVector3 lpos2,
                TVector3 lmom2,
                Double_t tof,
                Double_t length,
                Double_t eLoss);

    /** Copy constructor **/
    R3BDchPoint(const R3BDchPoint&);

    R3BDchPoint& operator=(const R3BDchPoint&) { return *this; }

    /** Destructor **/
    virtual ~R3BDchPoint();

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

    Double_t GetXLocalIn() const { return fLocalX_in; }
    Double_t GetYLocalIn() const { return fLocalY_in; }

    Double_t GetXLocalOut() const { return fLocalX_out; }
    Double_t GetYLocalOut() const { return fLocalY_out; }

    void PositionIn(TVector3& pos) { pos.SetXYZ(fX, fY, fZ); }
    void MomentumIn(TVector3& mom) { mom.SetXYZ(fPx, fPy, fPz); }

    void PositionOut(TVector3& pos) { pos.SetXYZ(fX_out, fY_out, fZ_out); }
    void MomentumOut(TVector3& mom) { mom.SetXYZ(fPx_out, fPy_out, fPz_out); }

    void LocalPositionIn(TVector3& pos) { pos.SetXYZ(fLocalX_in, fLocalY_in, fLocalZ_in); }
    void LocalMomentumIn(TVector3& mom) { mom.SetXYZ(fLocalPx_in, fLocalPy_in, fLocalPz_in); }

    void LocalPositionOut(TVector3& pos) { pos.SetXYZ(fLocalX_out, fLocalY_out, fLocalZ_out); }
    void LocalMomentumOut(TVector3& mom) { mom.SetXYZ(fLocalPx_out, fLocalPy_out, fLocalPz_out); }

    /** Point coordinates at given z from linear extrapolation **/
    Double_t GetX(Double_t z) const;
    Double_t GetY(Double_t z) const;
    Double_t GetLocalX(Double_t z) const;
    Double_t GetLocalY(Double_t z) const;

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
    Double32_t fLocalX_in, fLocalY_in, fLocalZ_in;
    Double32_t fLocalPx_in, fLocalPy_in, fLocalPz_in;
    Double32_t fLocalX_out, fLocalY_out, fLocalZ_out;
    Double32_t fLocalPx_out, fLocalPy_out, fLocalPz_out;

    Int_t fModule;
    Int_t fLayer;
    Int_t fCell;

    ClassDef(R3BDchPoint, 1)
};

inline void R3BDchPoint::SetPositionOut(TVector3 pos)
{
    fX_out = pos.X();
    fY_out = pos.Y();
    fZ_out = pos.Z();
}

inline void R3BDchPoint::SetMomentumOut(TVector3 mom)
{
    fPx_out = mom.Px();
    fPy_out = mom.Py();
    fPz_out = mom.Pz();
}

#endif
