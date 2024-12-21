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
// -----                      R3BTofdPoint source file                 -----
// -------------------------------------------------------------------------

#include "R3BTofdPoint.h"
#include "R3BLogger.h"

// -----   Default constructor   -------------------------------------------
R3BTofdPoint::R3BTofdPoint()
    : FairMCPoint()
{
    fPlane = 0;
    fPaddle = 0;
    fZfrag = 0.;
    fAfrag = 0.;
    fX_out = fY_out = fZ_out = 0.;
    fPx_out = fPy_out = fPz_out = 0.;
}
// -------------------------------------------------------------------------

// -----   Standard constructor   ------------------------------------------
R3BTofdPoint::R3BTofdPoint(Int_t trackID,
                           Int_t detID,
                           TVector3 posIn,
                           TVector3 posOut,
                           TVector3 momIn,
                           TVector3 momOut,
                           Double_t tof,
                           Double_t length,
                           Double_t eLoss)
    : R3BTofdPoint(trackID, detID, 0, 0, posIn, posOut, momIn, momOut, tof, length, eLoss, 0, 0)
{
}

// -----   Standard constructor   ------------------------------------------
R3BTofdPoint::R3BTofdPoint(Int_t trackID,
                           Int_t detID,
                           Int_t planeID,
                           Int_t paddleID,
                           TVector3 posIn,
                           TVector3 posOut,
                           TVector3 momIn,
                           TVector3 momOut,
                           Double_t tof,
                           Double_t length,
                           Double_t eLoss,
                           Double_t Z,
                           Double_t A)
    : FairMCPoint(trackID, detID, posIn, momIn, tof, length, eLoss)
{
    fPlane = planeID;
    fPaddle = paddleID;
    fX_out = posOut.X();
    fY_out = posOut.Y();
    fZ_out = posOut.Z();
    fPx_out = momOut.Px();
    fPy_out = momOut.Py();
    fPz_out = momOut.Pz();
    fZfrag = Z;
    fAfrag = A;
}
// -------------------------------------------------------------------------

// -----   Destructor   ----------------------------------------------------
R3BTofdPoint::~R3BTofdPoint() {}
// -------------------------------------------------------------------------

// -----   Public method Print   -------------------------------------------
void R3BTofdPoint::Print(const Option_t* opt) const
{
    R3BLOG(info,
           "Point for track " << fTrackID << " in plane " << fPlane << " and paddle " << fPaddle << "    Position ("
                              << fX << ", " << fY << ", " << fZ << ") cm"
                              << "    Momentum (" << fPx << ", " << fPy << ", " << fPz << ") GeV"
                              << "    Time " << fTime << " ns,  Length " << fLength << " cm,  Energy loss "
                              << fELoss * 1.0e06 << " keV");
}
// -------------------------------------------------------------------------

// -----   Point x coordinate from linear extrapolation   ------------------
Double_t R3BTofdPoint::GetX(Double_t z) const
{
    //  cout << fZ << " " << z << " " << fZ_out << endl;
    if ((fZ_out - z) * (fZ - z) >= 0.)
        return (fX_out + fX) / 2.;
    Double_t dz = fZ_out - fZ;
    return (fX + (z - fZ) / dz * (fX_out - fX));
}
// -------------------------------------------------------------------------

// -----   Point y coordinate from linear extrapolation   ------------------
Double_t R3BTofdPoint::GetY(Double_t z) const
{
    if ((fZ_out - z) * (fZ - z) >= 0.)
        return (fY_out + fY) / 2.;
    Double_t dz = fZ_out - fZ;
    //  if ( TMath::Abs(dz) < 1.e-3 ) return (fY_out+fY)/2.;
    return (fY + (z - fZ) / dz * (fY_out - fY));
}
// -------------------------------------------------------------------------

// -----   Public method IsUsable   ----------------------------------------
Bool_t R3BTofdPoint::IsUsable() const
{
    Double_t dz = fZ_out - fZ;
    if (TMath::Abs(dz) < 1.e-4)
        return kFALSE;
    return kTRUE;
}
// -------------------------------------------------------------------------

ClassImp(R3BTofdPoint);
