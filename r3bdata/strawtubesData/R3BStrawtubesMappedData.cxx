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

#include "R3BStrawtubesMappedData.h"

R3BStrawtubesMappedData::R3BStrawtubesMappedData(UChar_t plane,
                                                 UChar_t straw,
                                                 UChar_t side,
                                                 UInt_t timeCoarse,
                                                 UInt_t timeFine)
    : fPlane(plane)
    , fStraw(straw)
    , fSide(side)
    , fTimeCoarse(timeCoarse)
    , fTimeFine(timeFine)
{
}

R3BStrawtubesMappedData::~R3BStrawtubesMappedData() {}

ClassImp(R3BStrawtubesMappedData)
