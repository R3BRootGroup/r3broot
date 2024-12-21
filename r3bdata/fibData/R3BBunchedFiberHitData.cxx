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

#include "R3BBunchedFiberHitData.h"

R3BBunchedFiberHitData::R3BBunchedFiberHitData()
    : fFiberId(-1)
    , fMAPMTTime_ns(-1)
    , fSPMTTime_ns(-1)
    , fTime_ns(-1)
    , fMAPMTToT_ns(-1)
    , fSPMTToT_ns(-1)
{
}

R3BBunchedFiberHitData::R3BBunchedFiberHitData(Int_t ID,
                                               Double_t x,
                                               Double_t y,
                                               Double_t eloss,
                                               Double_t time,
                                               Int_t a_fiber_id,
                                               Double_t a_mapmt_time_ns,
                                               Double_t a_spmt_time_ns,
                                               Double_t a_mapmt_tot_ns,
                                               Double_t a_spmt_tot_ns)
    : R3BHit(ID, x, y, eloss, time)
    , fFiberId(a_fiber_id)
    , fMAPMTTime_ns(a_mapmt_time_ns)
    , fSPMTTime_ns(a_spmt_time_ns)
    , fTime_ns((a_mapmt_time_ns + a_spmt_time_ns) / 2)
    , fMAPMTToT_ns(a_mapmt_tot_ns)
    , fSPMTToT_ns(a_spmt_tot_ns)
{
}

R3BBunchedFiberHitData::~R3BBunchedFiberHitData() {}

Int_t R3BBunchedFiberHitData::GetFiberId() const { return fFiberId; }

Double_t R3BBunchedFiberHitData::GetMAPMTTime_ns() const { return fMAPMTTime_ns; }

Double_t R3BBunchedFiberHitData::GetSPMTTime_ns() const { return fSPMTTime_ns; }

Double_t R3BBunchedFiberHitData::GetTime_ns() const { return fTime_ns; }

Double_t R3BBunchedFiberHitData::GetMAPMTToT_ns() const { return fMAPMTToT_ns; }

Double_t R3BBunchedFiberHitData::GetSPMTToT_ns() const { return fSPMTToT_ns; }

ClassImp(R3BBunchedFiberHitData)
