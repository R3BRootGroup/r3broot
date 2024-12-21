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

#include "R3BLosHitData.h"

// for the data analysis of the Los detector.
// Introduced by Ralf, Feb 2016

R3BLosHitData::R3BLosHitData()
    : fDetector(0)
{
}

R3BLosHitData::R3BLosHitData(UChar_t detector, Double_t t, Double_t x, Double_t y, Double_t Z)

    : fDetector(detector)
    , fTime_ns(t)
    , fX_cm(x)
    , fY_cm(y)
    , fZ(Z)
//    , fTimeDiff2Ref_ns( diff )
{
}

ClassImp(R3BLosHitData)
