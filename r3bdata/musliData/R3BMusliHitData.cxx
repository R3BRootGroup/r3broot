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

#include "R3BMusliHitData.h"

// -----   Default constructor   -------------------------------------------
R3BMusliHitData::R3BMusliHitData()
    : fType(0)
    , fE(-1.)
    , fZ(-1.)
    , fX(0.)
    , fTheta(-1.)
{
}
// -------------------------------------------------------------------------

// -----   Standard constructor   ------------------------------------------
R3BMusliHitData::R3BMusliHitData(UInt_t n, Double_t e, Double_t z, Double_t x, Double_t theta)
    : fType(n)
    , fE(e)
    , fZ(z)
    , fX(x)
    , fTheta(theta)
{
}
// -------------------------------------------------------------------------

ClassImp(R3BMusliHitData);
