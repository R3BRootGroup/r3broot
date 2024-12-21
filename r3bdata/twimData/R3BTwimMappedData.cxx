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

// -------------------------------------------------------------------------
// -----                      R3BTwimMappedData source file            -----
// -------------------------------------------------------------------------

#include "R3BTwimMappedData.h"

// -----   Default constructor   -------------------------------------------
R3BTwimMappedData::R3BTwimMappedData()
    : fSecID(0)
    , fAnodeID(0)
    , fTime(0)
    , fEnergy(0)
    , fPileup(kFALSE)
    , fOverflow(kFALSE)
{
}
// -------------------------------------------------------------------------

// -----   Standard constructor   ------------------------------------------
R3BTwimMappedData::R3BTwimMappedData(UInt_t secID, UInt_t anodeID, Int_t time, Int_t energy, Bool_t pu, Bool_t ov)
    : fSecID(secID)
    , fAnodeID(anodeID)
    , fTime(time)
    , fEnergy(energy)
    , fPileup(pu)
    , fOverflow(ov)
{
}
// -------------------------------------------------------------------------

ClassImp(R3BTwimMappedData);
