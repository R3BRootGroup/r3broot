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

// -----------------------------------------------------------------------
// -----            R3BFootMappedData source file                    -----
// -----------------------------------------------------------------------

#include "R3BFootMappedData.h"

R3BFootMappedData::R3BFootMappedData(uint8_t detId, uint16_t stripId, int energy)
    : fDetId(detId)
    , fStripId(stripId)
    , fEnergy(energy)
{
}

ClassImp(R3BFootMappedData)
