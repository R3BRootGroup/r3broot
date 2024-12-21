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
// -----                      R3BDchPoint source file                  -----
// -------------------------------------------------------------------------

#include "R3BDchPoint.h"

#include <iostream>

using std::cout;
using std::endl;
using std::flush;

// -----   Default constructor   -------------------------------------------
R3BDchPoint::R3BDchPoint()
    : FairMCPoint()
    , fX_out(-99999.)
    , fY_out(-99999.)
    , fZ_out(-99999.)
    , fPx_out(-99999.)
    , fPy_out(-99999.)
    , fPz_out(-99999.)
    , fLocalX_in(-99999.)
    , fLocalY_in(-99999.)
    , fLocalZ_in(-99999.)
    , fLocalPx_in(-99999.)
    , fLocalPy_in(-99999.)
    , fLocalPz_in(-99999.)
    , fLocalX_out(-99999.)
    , fLocalY_out(-99999.)
    , fLocalZ_out(-99999.)
    , fLocalPx_out(-99999.)
    , fLocalPy_out(-99999.)
    , fLocalPz_out(-99999.)
    , fModule(-99999)
    , fLayer(-99999)
    , fCell(-99999)
{
}
// -------------------------------------------------------------------------

// -----   Standard constructor   ------------------------------------------
R3BDchPoint::R3BDchPoint(Int_t trackId,
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
                         Double_t eLoss)
    : FairMCPoint(trackId, mod, posIn, momIn, tof, length, eLoss)
    , fX_out(posOut.X())
    , fY_out(posOut.Y())
    , fZ_out(posOut.Z())
    , fPx_out(momOut.Px())
    , fPy_out(momOut.Py())
    , fPz_out(momOut.Pz())
    , fLocalX_in(lpos1.X())
    , fLocalY_in(lpos1.Y())
    , fLocalZ_in(lpos1.Z())
    , fLocalPx_in(lmom1.Px())
    , fLocalPy_in(lmom1.Py())
    , fLocalPz_in(lmom1.Pz())
    , fLocalX_out(lpos2.X())
    , fLocalY_out(lpos2.Y())
    , fLocalZ_out(lpos2.Z())
    , fLocalPx_out(lmom2.Px())
    , fLocalPy_out(lmom2.Py())
    , fLocalPz_out(lmom2.Pz())
    , fModule(mod)
    , fLayer(layer)
    , fCell(cell)
{
}

R3BDchPoint::R3BDchPoint(const R3BDchPoint& right)
    : FairMCPoint(right)
    , fX_out(right.fX_out)
    , fY_out(right.fY_out)
    , fZ_out(right.fZ_out)
    , fPx_out(right.fPx_out)
    , fPy_out(right.fPy_out)
    , fPz_out(right.fPz_out)
    , fLocalX_in(right.fLocalX_in)
    , fLocalY_in(right.fLocalY_in)
    , fLocalZ_in(right.fLocalZ_in)
    , fLocalPx_in(right.fLocalPx_in)
    , fLocalPy_in(right.fLocalPy_in)
    , fLocalPz_in(right.fLocalPz_in)
    , fLocalX_out(right.fLocalX_out)
    , fLocalY_out(right.fLocalY_out)
    , fLocalZ_out(right.fLocalZ_out)
    , fLocalPx_out(right.fLocalPx_out)
    , fLocalPy_out(right.fLocalPy_out)
    , fLocalPz_out(right.fLocalPz_out)
    , fModule(right.fModule)
    , fLayer(right.fLayer)
    , fCell(right.fCell)
{
}

// -----   Destructor   ----------------------------------------------------
R3BDchPoint::~R3BDchPoint() {}
// -------------------------------------------------------------------------

// -----   Public method Print   -------------------------------------------
void R3BDchPoint::Print(const Option_t* opt) const
{
    cout << "-I- R3BDchPoint: STS Point for track " << fTrackID << " in detector " << fDetectorID << endl;
    cout << "    Position (" << fX << ", " << fY << ", " << fZ << ") cm" << endl;
    cout << "    Momentum (" << fPx << ", " << fPy << ", " << fPz << ") GeV" << endl;
    cout << "    Time " << fTime << " ns,  Length " << fLength << " cm,  Energy loss " << fELoss * 1.0e06 << " keV"
         << endl;
}
// -------------------------------------------------------------------------

// -----   Local Point x coordinate from linear extrapolation   ------------------
Double_t R3BDchPoint::GetLocalX(Double_t z) const
{
    //  cout << fZ << " " << z << " " << fZ_out << endl;
    if ((fLocalZ_out - z) * (fLocalZ_in - z) >= 0.)
        return (fLocalX_out + fLocalX_in) / 2.;
    Double_t dz = fLocalZ_out - fLocalZ_in;
    return (fLocalX_in + (z - fLocalZ_in) / dz * (fLocalX_out - fLocalX_in));
}
// -------------------------------------------------------------------------

// -----   Local Point y coordinate from linear extrapolation   ------------------
Double_t R3BDchPoint::GetLocalY(Double_t z) const
{
    if ((fLocalZ_out - z) * (fLocalZ_in - z) >= 0.)
        return (fLocalY_out + fLocalY_in) / 2.;
    Double_t dz = fLocalZ_out - fLocalZ_in;
    //  if ( TMath::Abs(dz) < 1.e-3 ) return (fY_out+fY)/2.;
    return (fLocalY_in + (z - fLocalZ_in) / dz * (fLocalY_out - fLocalY_in));
}
// -------------------------------------------------------------------------

// -----   Point x coordinate from linear extrapolation   ------------------
Double_t R3BDchPoint::GetX(Double_t z) const
{
    //  cout << fZ << " " << z << " " << fZ_out << endl;
    if ((fZ_out - z) * (fZ - z) >= 0.)
        return (fX_out + fX) / 2.;
    Double_t dz = fZ_out - fZ;
    return (fX + (z - fZ) / dz * (fX_out - fX));
}
// -------------------------------------------------------------------------

// -----   Point y coordinate from linear extrapolation   ------------------
Double_t R3BDchPoint::GetY(Double_t z) const
{
    if ((fZ_out - z) * (fZ - z) >= 0.)
        return (fY_out + fY) / 2.;
    Double_t dz = fZ_out - fZ;
    //  if ( TMath::Abs(dz) < 1.e-3 ) return (fY_out+fY)/2.;
    return (fY + (z - fZ) / dz * (fY_out - fY));
}
// -------------------------------------------------------------------------

// -----   Public method IsUsable   ----------------------------------------
Bool_t R3BDchPoint::IsUsable() const
{
    Double_t dz = fZ_out - fZ;
    if (TMath::Abs(dz) < 1.e-4)
        return kFALSE;
    return kTRUE;
}
// -------------------------------------------------------------------------

ClassImp(R3BDchPoint)
