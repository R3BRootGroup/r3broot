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

// ----------------------------------------------------------------------
// -----              R3BMwpcMappedData source file                 -----
// ----------------------------------------------------------------------

#include "R3BMwpcMappedData.h"

// -----   Default constructor   -------------------------------------------
R3BMwpcMappedData::R3BMwpcMappedData()
    : fPlane(0)
    , fPad(0)
    , fCharge(0)
{
}
// -------------------------------------------------------------------------

// -----   Standard constructor   ------------------------------------------
R3BMwpcMappedData::R3BMwpcMappedData(UInt_t plane, UInt_t pad, UInt_t charge)
    : fPlane(plane)
    , fPad(pad)
    , fCharge(charge)
{
}
// -------------------------------------------------------------------------

ClassImp(R3BMwpcMappedData);
