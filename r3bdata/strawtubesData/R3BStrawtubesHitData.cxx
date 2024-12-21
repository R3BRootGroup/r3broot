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

#include "R3BStrawtubesHitData.h"
#include "R3BStrawtubesCalData.h"

R3BStrawtubesHitData::R3BStrawtubesHitData(R3BStrawtubesCalData const& cal)
    : fPlane(cal.GetPlane())
    , fStraw(cal.GetStraw())
    , fAverage_ns(0.5 * (cal.GetTime(1) + cal.GetTime(0)))
    , fDiff_ns(cal.GetTime(1) - cal.GetTime(0))
{
}

R3BStrawtubesHitData::~R3BStrawtubesHitData() {}

ClassImp(R3BStrawtubesHitData)
