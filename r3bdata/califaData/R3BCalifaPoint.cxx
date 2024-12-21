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

#include "R3BCalifaPoint.h"

R3BCalifaPoint::R3BCalifaPoint(int trackID,
                               int detID,
                               int trackPID,
                               int ident,
                               TVector3 posIn,
                               TVector3 momIn,
                               double tof,
                               double length,
                               double eLoss)
    : FairMCPoint(trackID, detID, posIn, momIn, tof, length, eLoss)
    , fCrystalId(ident)
    , fTrackPID(trackPID)
{
}

R3BCalifaPoint::R3BCalifaPoint(const R3BCalifaPoint& right)
    : FairMCPoint(right)
    , fCrystalId(right.fCrystalId)
    , fTrackPID(right.fTrackPID)
{
}

ClassImp(R3BCalifaPoint)
