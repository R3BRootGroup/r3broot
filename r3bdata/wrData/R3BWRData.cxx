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

// ---------------------------------------------------------------------
// -----                      R3BWRData                            -----
// -----          Created 28/02/2019 by J.L. Rodriguez             -----
// ---------------------------------------------------------------------

#include "R3BWRData.h"

R3BWRData::R3BWRData()
    : fTimeStamp(0)
    , fId(1)
{
}

//------------------------------

R3BWRData::R3BWRData(uint64_t timestamp, UInt_t id)
    : fTimeStamp(timestamp)
    , fId(id)
{
}

ClassImp(R3BWRData);
