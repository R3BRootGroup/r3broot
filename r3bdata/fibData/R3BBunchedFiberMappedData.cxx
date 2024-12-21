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

// ----------------------------------------------------------------
// -----              R3BBunchedFiberMappedData                -----
// -----             Created Jan 2018 by M.Heil        -----
// ----------------------------------------------------------------

#include "R3BBunchedFiberMappedData.h"

R3BBunchedFiberMappedData::R3BBunchedFiberMappedData()
    : fSide(-1)
    , fChannel(-1)
    , fIsLeading()
    , fCoarse(-1)
    , fFine(-1)
{
}

R3BBunchedFiberMappedData::R3BBunchedFiberMappedData(Int_t a_side,
                                                     Int_t a_channel,
                                                     Bool_t a_is_leading,
                                                     Int_t a_coarse,
                                                     Int_t a_fine)
    : fSide(a_side)
    , fChannel(a_channel)
    , fIsLeading(a_is_leading)
    , fCoarse(a_coarse)
    , fFine(a_fine)
{
}

R3BBunchedFiberMappedData::~R3BBunchedFiberMappedData() {}

Int_t R3BBunchedFiberMappedData::GetChannel() const { return fChannel; }

Int_t R3BBunchedFiberMappedData::GetCoarse() const { return fCoarse; }

Int_t R3BBunchedFiberMappedData::GetFine() const { return fFine; }

Int_t R3BBunchedFiberMappedData::GetSide() const { return fSide; }

Bool_t R3BBunchedFiberMappedData::IsMAPMT() const { return 0 == fSide; }

Bool_t R3BBunchedFiberMappedData::IsSPMT() const { return 1 == fSide; }

Bool_t R3BBunchedFiberMappedData::IsMAPMTTrigger() const { return 2 == fSide; }

Bool_t R3BBunchedFiberMappedData::IsSPMTTrigger() const { return 3 == fSide; }

Bool_t R3BBunchedFiberMappedData::IsLeading() const { return fIsLeading; }

Bool_t R3BBunchedFiberMappedData::IsTrailing() const { return !fIsLeading; }

ClassImp(R3BBunchedFiberMappedData)
