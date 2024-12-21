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

#include "R3BSfibHitData.h"

R3BSfibHitData::R3BSfibHitData()
    : fFiberId(-1)
    , fTopTime_ns(-1)
    , fBotTime_ns(-1)
    , fTime_ns(-1)
    , fTopToT_ns(-1)
    , fBotToT_ns(-1)
{
}

R3BSfibHitData::R3BSfibHitData(Int_t a_fiber_id,
                               Double_t a_top_time_ns,
                               Double_t a_bot_time_ns,
                               Double_t a_top_tot_ns,
                               Double_t a_bot_tot_ns)
    : R3BHit(0, 0, 0, 0, 0)
    , fFiberId(a_fiber_id)
    , fTopTime_ns(a_top_time_ns)
    , fBotTime_ns(a_bot_time_ns)
    , fTime_ns((a_top_time_ns + a_bot_time_ns) / 2)
    , fTopToT_ns(a_top_tot_ns)
    , fBotToT_ns(a_bot_tot_ns)
{
}

R3BSfibHitData::~R3BSfibHitData() {}

Int_t R3BSfibHitData::GetFiberId() const { return fFiberId; }

Double_t R3BSfibHitData::GetTopTime_ns() const { return fTopTime_ns; }

Double_t R3BSfibHitData::GetBotTime_ns() const { return fBotTime_ns; }

Double_t R3BSfibHitData::GetTime_ns() const { return fTime_ns; }

Double_t R3BSfibHitData::GetTopToT_ns() const { return fTopToT_ns; }

Double_t R3BSfibHitData::GetBotToT_ns() const { return fBotToT_ns; }

ClassImp(R3BSfibHitData)
