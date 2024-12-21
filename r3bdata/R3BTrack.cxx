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

// -----------------------------------------------------------------------------
// -----                              R3BTrack                             -----
// -----                 Created on 09.03.2020 by M.Heil                   -----
// -----------------------------------------------------------------------------

#include "R3BTrack.h"

R3BTrack::R3BTrack()
    : fX(0.)
    , fY(0.)
    , fZ(0.)
    , fPx(0.)
    , fPy(0.)
    , fPz(0.)
    , fQ(0)
    , fAoZ(0.)
    , fChix(0.)
    , fChiy(0.)
    , fQuality(0.)
{
}

R3BTrack::R3BTrack(Double_t x,
                   Double_t y,
                   Double_t z,
                   Double_t px,
                   Double_t py,
                   Double_t pz,
                   Double_t q,
                   Double_t AoZ,
                   Double_t chix,
                   Double_t chiy,
                   Int_t quality)
    : fX(x)
    , fY(y)
    , fZ(z)
    , fPx(px)
    , fPy(py)
    , fPz(pz)
    , fQ(q)
    , fAoZ(AoZ)
    , fChix(chix)
    , fChiy(chiy)
    , fQuality(quality)
{
}

R3BTrack::~R3BTrack() {}

ClassImp(R3BTrack)
