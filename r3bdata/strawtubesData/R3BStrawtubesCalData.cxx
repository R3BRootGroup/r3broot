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

#include "R3BStrawtubesCalData.h"
#include <cmath>

R3BStrawtubesCalData::R3BStrawtubesCalData(UChar_t plane, UChar_t straw)
    : fPlane(plane)
    , fStraw(straw)
    , fTime_ns()
{
    fTime_ns[0] = NAN;
    fTime_ns[1] = NAN;
}

R3BStrawtubesCalData::~R3BStrawtubesCalData() {}

UChar_t R3BStrawtubesCalData::GetPlane() const { return fPlane; }

UChar_t R3BStrawtubesCalData::GetStraw() const { return fStraw; }

Double_t R3BStrawtubesCalData::GetTime(UChar_t side) const { return fTime_ns[side]; }

void R3BStrawtubesCalData::SetTime(UChar_t side, Double_t time) { fTime_ns[side] = time; }

ClassImp(R3BStrawtubesCalData)
