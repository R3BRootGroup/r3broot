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
// -----                      R3BTraPoint header file                  -----
// -----                  Created 26/03/09  by D.Bertini               -----
// -------------------------------------------------------------------------

/**  R3BTraPoint.h
 **/

#ifndef R3BTRAPOINT_H
#define R3BTRAPOINT_H 1

#include "TObject.h"
#include "TVector3.h"

#include "FairMCPoint.h"

class R3BTraPoint : public FairMCPoint
{
  public:
    /** Default constructor **/
    R3BTraPoint();

    /** Constructor with arguments
     *@param trackID  Index of MCTrack
     *@param detID    Detector ID
     *@param detVolID Detector Copy ID  // added by Marc
     *@param posIn    Ccoordinates at entrance to active volume [cm]
     *@param posOut   Coordinates at exit of active volume [cm]
     *@param momIn    Momentum of track at entrance [GeV]
     *@param momOut   Momentum of track at exit [GeV]
     *@param tof      Time since event start [ns]
     *@param length   Track length since creation [cm]
     *@param eLoss    Energy deposit [GeV]
     *@param pid      Particle ID
     **/
    R3BTraPoint(Int_t trackID,
                Int_t detID,
                Int_t detCopyID,
                TVector3 posIn,
                TVector3 posOut,
                TVector3 momIn,
                TVector3 momOut,
                Double_t tof,
                Double_t length,
                Double_t eLoss,
                Int_t pid);

    /** Copy constructor **/
    R3BTraPoint(const R3BTraPoint& point) { *this = point; };

    /** Destructor **/
    virtual ~R3BTraPoint();

    /** Accessors **/
    Int_t GetDetCopyID() const { return fDetCopyID; }
    Double_t GetXIn() const { return fX; }
    Double_t GetYIn() const { return fY; }
    Int_t GetPid() const { return fPid; }
    Double_t GetZIn() const { return fZ; }
    Double_t GetXOut() const { return fX_out; }
    Double_t GetYOut() const { return fY_out; }
    Double_t GetZOut() const { return fZ_out; }
    Double_t GetPxOut() const { return fPx_out; }
    Double_t GetPyOut() const { return fPy_out; }
    Double_t GetPzOut() const { return fPz_out; }
    inline Double_t GetTheta() const
    {
        TVector3 v;
        v.SetXYZ(fX, fY, fZ);
        return v.Theta();
    }
    inline Double_t GetPhi() const
    {
        TVector3 v;
        v.SetXYZ(fX, fY, fZ);
        return v.Phi();
    }

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
    void SetDetCopyID(Int_t id) { fDetCopyID = id; }; // added by Marc

    /** Output to screen **/
    virtual void Print(const Option_t* opt) const;

  protected:
    Double32_t fX_out, fY_out, fZ_out, fEloss;
    Double32_t fPx_out, fPy_out, fPz_out;
    Int_t fDetCopyID, fPid;

    ClassDef(R3BTraPoint, 1)
};

inline void R3BTraPoint::SetPositionOut(TVector3 pos)
{
    fX_out = pos.X();
    fY_out = pos.Y();
    fZ_out = pos.Z();
}

inline void R3BTraPoint::SetMomentumOut(TVector3 mom)
{
    fPx_out = mom.Px();
    fPy_out = mom.Py();
    fPz_out = mom.Pz();
}

#endif
