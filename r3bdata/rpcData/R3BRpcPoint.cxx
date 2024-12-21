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

#include "R3BRpcPoint.h"

#include <iostream>

using std::cout;
using std::endl;
using std::flush;

R3BRpcPoint::R3BRpcPoint()
    : FairMCPoint()
    , fChannelId(-1)
{
}

R3BRpcPoint::R3BRpcPoint(Int_t trackID,
                         Int_t detID,
                         Int_t ident,
                         TVector3 posIn,
                         TVector3 momIn,
                         Double_t tof,
                         Double_t length,
                         Double_t eLoss,
                         UInt_t EventId)
    : FairMCPoint(trackID, detID, posIn, momIn, tof, length, eLoss, EventId)
    , fChannelId(ident)
{
}

R3BRpcPoint::R3BRpcPoint(const R3BRpcPoint& right)
    : FairMCPoint(right)
    , fChannelId(right.fChannelId)
{
}

R3BRpcPoint::~R3BRpcPoint() {}

void R3BRpcPoint::Print(const Option_t* opt) const
{
    cout << "-I- R3BRpcPoint: RPC Point for track " << fTrackID << " in channel " << fChannelId << endl;
    cout << "    Position (" << fX << ", " << fY << ", " << fZ << ") cm" << endl;
    cout << "    Momentum (" << fPx << ", " << fPy << ", " << fPz << ") GeV" << endl;
    cout << "    Time " << fTime << " ns,  Length " << fLength << " cm,  Energy loss " << fELoss * 1.0e06 << " keV"
         << endl;
}

ClassImp(R3BRpcPoint)
