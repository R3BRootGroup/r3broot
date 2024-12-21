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
// -----                      R3BDchFullPoint source file                  -----
// -------------------------------------------------------------------------

#include "R3BDchFullPoint.h"

#include <iostream>

// -----   Default constructor   -------------------------------------------
R3BDchFullPoint::R3BDchFullPoint()
    : FairMCPoint()
    , fLocalX(0.)
    , fLocalY(0.)
    , fLocalZ(0.)
    , fLocalPx(0.)
    , fLocalPy(0.)
    , fLocalPz(0.)
    , fModule(-1)
    , fLayer(-1)
    , fCell(-1)
{
}
// -------------------------------------------------------------------------

// -----   Standard constructor   ------------------------------------------
R3BDchFullPoint::R3BDchFullPoint(Int_t trackId,
                                 Int_t mod,
                                 Int_t layer,
                                 Int_t cell,
                                 TVector3 pos,
                                 TVector3 lpos,
                                 TVector3 mom,
                                 TVector3 lmom,
                                 Double_t tof,
                                 Double_t length,
                                 Double_t eLoss)
    : FairMCPoint(trackId, mod, pos, mom, tof, length, eLoss)
    , fLocalX(lpos.X())
    , fLocalY(lpos.Y())
    , fLocalZ(lpos.Z())
    , fLocalPx(lmom.X())
    , fLocalPy(lmom.Y())
    , fLocalPz(lmom.Z())
    , fModule(mod)
    , fLayer(layer)
    , fCell(cell)
{
}

// -----   Destructor   ----------------------------------------------------
R3BDchFullPoint::~R3BDchFullPoint() {}
// -------------------------------------------------------------------------

ClassImp(R3BDchFullPoint)
