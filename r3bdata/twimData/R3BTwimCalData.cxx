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
// -----                      R3BTwimCalData source file              -----
// -------------------------------------------------------------------------

#include "R3BTwimCalData.h"

// -----   Default constructor   -------------------------------------------
R3BTwimCalData::R3BTwimCalData()
    : fSecID(0)
    , fAnodeID(0)
    , fDT(0.)
    , fEnergy(0.)
{
}
// -------------------------------------------------------------------------

// -----   Standard constructor   ------------------------------------------
R3BTwimCalData::R3BTwimCalData(UInt_t secID, UInt_t anodeID, Double_t dt, Double_t energy)
    : fSecID(secID)
    , fAnodeID(anodeID)
    , fDT(dt)
    , fEnergy(energy)
{
}
// -------------------------------------------------------------------------

ClassImp(R3BTwimCalData);
