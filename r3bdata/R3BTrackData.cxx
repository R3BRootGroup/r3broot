/******************************************************************************
 *   Copyright (C) 2019 GSI Helmholtzzentrum für Schwerionenforschung GmbH    *
 *   Copyright (C) 2019 Members of R3B Collaboration                          *
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
// -----                 Created on 03.03.2023 by A.Kelic-Heil             -----
// -----------------------------------------------------------------------------

#include "R3BTrackData.h"

R3BTrackData::R3BTrackData()
    : fDetPos()
    , fTime(0.)
{
}

R3BTrackData::R3BTrackData(Double_t x,
                           Double_t y,
                           Double_t z,
                           Double_t px,
                           Double_t py,
                           Double_t pz,
                           Int_t q,
                           Double_t AoZ,
                           Double_t chix,
                           Double_t chiy,
                           Int_t quality,
                           std::vector<Double_t> detPos,
                           Double_t t)
    : R3BTrack(x, y, z, px, py, pz, q, AoZ, chix, chiy, quality)
    , fDetPos(detPos)
    , fTime(t)
{
}

R3BTrackData::~R3BTrackData() {}

ClassImp(R3BTrackData)
