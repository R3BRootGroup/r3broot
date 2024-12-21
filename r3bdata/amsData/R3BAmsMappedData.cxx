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

// ---------------------------------------------------------------------------
// -----                                                                 -----
// -----                      R3BAmsMappedData                           -----
// -----                  Created 12/05/2018 by J.L. Rodriguez           -----
// -----                                                                 -----
// ---------------------------------------------------------------------------

#include "R3BAmsMappedData.h"

R3BAmsMappedData::R3BAmsMappedData()
    : fDetId(0)
    , fStripId(0)
    , fEnergy(0)
{
}

//------------------------------

R3BAmsMappedData::R3BAmsMappedData(Int_t detId, Int_t stripId, Int_t energy)
    : fDetId(detId)
    , fStripId(stripId)
    , fEnergy(energy)
{
}

ClassImp(R3BAmsMappedData);
