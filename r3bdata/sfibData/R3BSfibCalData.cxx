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

#include "R3BSfibCalData.h"

R3BSfibCalData::R3BSfibCalData()
    : fSide(-1)
    , fChannel(-1)
    , fIsLeading()
    , fTime_ns(-1)
{
}

R3BSfibCalData::R3BSfibCalData(Int_t a_side, Int_t a_channel, Bool_t a_is_leading, Double_t a_time_ns)
    : fSide(a_side)
    , fChannel(a_channel)
    , fIsLeading(a_is_leading)
    , fTime_ns(a_time_ns)
{
}

R3BSfibCalData::~R3BSfibCalData() {}

Int_t R3BSfibCalData::GetChannel() const { return fChannel; }

Double_t R3BSfibCalData::GetTime_ns() const { return fTime_ns; }

Bool_t R3BSfibCalData::IsTop() const { return 2 == fSide; }

Bool_t R3BSfibCalData::IsBot() const { return !IsTop(); }

Bool_t R3BSfibCalData::IsLeading() const { return fIsLeading; }

Bool_t R3BSfibCalData::IsTrailing() const { return !fIsLeading; }

ClassImp(R3BSfibCalData)
