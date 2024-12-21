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

// --------------------------------------------------------------
// -----            R3BTrloiiData source file               -----
// --------------------------------------------------------------

#include "R3BTrloiiData.h"

R3BTrloiiData::R3BTrloiiData()
    : fType(0)
    , fCh(0)
    , fCounts(0)
{
}

R3BTrloiiData::R3BTrloiiData(UInt_t type, Int_t ch, uint32_t counts)
    : fType(type)
    , fCh(ch)
    , fCounts(counts)
{
}

ClassImp(R3BTrloiiData);
