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

#include "R3BFiberMAPMTHitData.h"

R3BFiberMAPMTHitData::R3BFiberMAPMTHitData()
    : fFiberId(0)
    , fBottomTime_ns(-1)
    , fTopTime_ns(-1)
    , fTime_ns(-1)
    , fBottomToT_ns(-1)
    , fTopToT_ns(-1)
{
}

R3BFiberMAPMTHitData::R3BFiberMAPMTHitData(Int_t ID,
                                           Double_t x,
                                           Double_t y,
                                           Double_t eloss,
                                           Double_t time,
                                           Int_t a_fiber_id,
                                           Double_t a_bottom_time_ns,
                                           Double_t a_top_time_ns,
                                           Double_t a_bottom_tot_ns,
                                           Double_t a_top_tot_ns)
    : R3BHit(ID, x, y, eloss, time)
    , fFiberId(a_fiber_id)
    , fBottomTime_ns(a_bottom_time_ns)
    , fTopTime_ns(a_top_time_ns)
    , fTime_ns((a_bottom_time_ns + a_top_time_ns) / 2.)
    , fBottomToT_ns(a_bottom_tot_ns)
    , fTopToT_ns(a_top_tot_ns)
{
}

R3BFiberMAPMTHitData::~R3BFiberMAPMTHitData() {}

Int_t R3BFiberMAPMTHitData::GetFiberId() const { return fFiberId; }

Double_t R3BFiberMAPMTHitData::GetBottomTime_ns() const { return fBottomTime_ns; }

Double_t R3BFiberMAPMTHitData::GetTopTime_ns() const { return fTopTime_ns; }

Double_t R3BFiberMAPMTHitData::GetTime_ns() const { return fTime_ns; }

Double_t R3BFiberMAPMTHitData::GetBottomToT_ns() const { return fBottomToT_ns; }

Double_t R3BFiberMAPMTHitData::GetTopToT_ns() const { return fTopToT_ns; }

ClassImp(R3BFiberMAPMTHitData);
