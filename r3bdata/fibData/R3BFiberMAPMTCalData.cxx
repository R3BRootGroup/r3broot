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

#include "R3BFiberMAPMTCalData.h"

R3BFiberMAPMTCalData::R3BFiberMAPMTCalData()
    : fSide(-1)
    , fChannel(-1)
    , fIsLeading()
    , fTime_ns(-1)
{
}

R3BFiberMAPMTCalData::R3BFiberMAPMTCalData(Int_t a_side, Int_t a_channel, Bool_t a_is_leading, Double_t a_time_ns)
    : fSide(a_side)
    , fChannel(a_channel)
    , fIsLeading(a_is_leading)
    , fTime_ns(a_time_ns)
{
}

R3BFiberMAPMTCalData::~R3BFiberMAPMTCalData() {}

Int_t R3BFiberMAPMTCalData::GetChannel() const { return fChannel; }

Double_t R3BFiberMAPMTCalData::GetTime_ns() const { return fTime_ns; }

Int_t R3BFiberMAPMTCalData::GetSide() const { return fSide; }

Double_t R3BFiberMAPMTCalData::GetTimeLeading_ns() const { return fIsLeading ? fTime_ns : 0; }

Double_t R3BFiberMAPMTCalData::GetTimeTrailing_ns() const { return fIsLeading ? 0 : fTime_ns; }

Bool_t R3BFiberMAPMTCalData::IsBottom() const { return 0 == fSide; }

Bool_t R3BFiberMAPMTCalData::IsTop() const { return 1 == fSide; }

Bool_t R3BFiberMAPMTCalData::IsTrigger() const { return 2 == fSide; }

Bool_t R3BFiberMAPMTCalData::IsLeading() const { return fIsLeading; }

Bool_t R3BFiberMAPMTCalData::IsTrailing() const { return !fIsLeading; }

Bool_t R3BFiberMAPMTCalData::IsSortable() const { return kTRUE; }

ClassImp(R3BFiberMAPMTCalData)
