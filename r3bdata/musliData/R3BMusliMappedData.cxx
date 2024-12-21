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

#include "R3BMusliMappedData.h"

// -----   Default constructor   -------------------------------------------
R3BMusliMappedData::R3BMusliMappedData()
    : fType(0)
    , fIndex(0)
    , fSignal(0)
    , fTime(0)
    , fEnergy(0)
    , fPileup(kFALSE)
    , fOverflow(kFALSE)
{
}
// -------------------------------------------------------------------------

// -----   Standard constructor   ------------------------------------------
R3BMusliMappedData::R3BMusliMappedData(UInt_t type,
                                       UInt_t index,
                                       UInt_t signal,
                                       UInt_t time,
                                       UInt_t energy,
                                       Bool_t pu,
                                       Bool_t ov)
    : fType(type)
    , fIndex(index)
    , fSignal(signal)
    , fTime(time)
    , fEnergy(energy)
    , fPileup(pu)
    , fOverflow(ov)
{
}
// -------------------------------------------------------------------------

ClassImp(R3BMusliMappedData);
