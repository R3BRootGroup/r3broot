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

#ifndef R3BNEUTRONTRACK_H
#define R3BNEUTRONTRACK_H

#include "TObject.h"
#include "TVector3.h"

class R3BNeutronTrack : public TObject
{

  public:
    /** Default constructor **/
    R3BNeutronTrack();

    R3BNeutronTrack(TVector3 posIn, TVector3 posOut, TVector3 momOut, Double_t tof);

    /** Copy constructor **/
    R3BNeutronTrack(const R3BNeutronTrack&);

    /** Destructor **/
    virtual ~R3BNeutronTrack();

    R3BNeutronTrack& operator=(const R3BNeutronTrack&) { return *this; }

    /** Output to screen **/
    virtual void Print(const Option_t* opt) const;

    Double_t GetXIn() const { return fX; }
    Double_t GetYIn() const { return fY; }
    Double_t GetZIn() const { return fZ; }
    Double_t GetXOut() const { return fX_out; }
    Double_t GetYOut() const { return fY_out; }
    Double_t GetZOut() const { return fZ_out; }
    Double_t GetPxOut() const { return fPx_out; }
    Double_t GetPyOut() const { return fPy_out; }
    Double_t GetPzOut() const { return fPz_out; }

    void PositionIn(TVector3& pos) { pos.SetXYZ(fX, fY, fZ); }
    void PositionOut(TVector3& pos) { pos.SetXYZ(fX_out, fY_out, fZ_out); }
    void MomentumOut(TVector3& mom) { mom.SetXYZ(fPx_out, fPy_out, fPz_out); }

    /** Modifiers **/
    void SetPositionIn(TVector3 pos);
    void SetPositionOut(TVector3 pos);
    void SetMomentumOut(TVector3 mom);

  protected:
    Double32_t fX, fY, fZ;
    Double32_t fX_in, fY_in, fZ_in;
    Double32_t fX_out, fY_out, fZ_out;
    Double32_t fPx_out, fPy_out, fPz_out;
    Int_t fPaddleNb;

    ClassDef(R3BNeutronTrack, 1)
};

inline void R3BNeutronTrack::SetPositionIn(TVector3 pos)
{
    fX = pos.X();
    fY = pos.Y();
    fZ = pos.Z();
}

inline void R3BNeutronTrack::SetPositionOut(TVector3 pos)
{
    fX_out = pos.X();
    fY_out = pos.Y();
    fZ_out = pos.Z();
}

inline void R3BNeutronTrack::SetMomentumOut(TVector3 mom)
{
    fPx_out = mom.Px();
    fPy_out = mom.Py();
    fPz_out = mom.Pz();
}

#endif
