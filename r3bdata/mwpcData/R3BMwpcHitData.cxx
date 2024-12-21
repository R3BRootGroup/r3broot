/******************************************************************************
 **   Copyright (C) 2019 GSI Helmholtzzentrum für Schwerionenforschung GmbH    *
 **   Copyright (C) 2019-2025 Members of R3B Collaboration                     *
 **                                                                            *
 **             This software is distributed under the terms of the            *
 **                 GNU General Public Licence (GPL) version 3,                *
 **                    copied verbatim in the file "LICENSE".                  *
 **                                                                            *
 ** In applying this license GSI does not waive the privileges and immunities  *
 ** granted to it by virtue of its status as an Intergovernmental Organization *
 ** or submit itself to any jurisdiction.                                      *
 *******************************************************************************/
// -------------------------------------------------------------------------
// -----                   R3BMwpcHitData source file                  -----
// -------------------------------------------------------------------------

#include "R3BMwpcHitData.h"

// -----   Default constructor   -------------------------------------------
R3BMwpcHitData::R3BMwpcHitData()
    : fX(0.)
    , fY(0.)
    , fPlane(0)
{
}
// -------------------------------------------------------------------------

// -----   Standard constructor   ------------------------------------------
R3BMwpcHitData::R3BMwpcHitData(Double_t x, Double_t y, Int_t plane)
    : fX(x)
    , fY(y)
    , fPlane(plane)
{
}
// -------------------------------------------------------------------------

ClassImp(R3BMwpcHitData);
