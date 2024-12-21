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

#include "R3BTofiHitData.h"

// for the data analysis of the Tofi detector.
// Introduced by M.Heil, May 2016

R3BTofiHitData::R3BTofiHitData()
    : fBarId(-1)
{
}

R3BTofiHitData::R3BTofiHitData(Double_t t,
                               Double_t x,
                               Double_t y,
                               Double_t Z,
                               Double_t tdiff,
                               Double_t ELoss,
                               Double_t ID,
                               Int_t iBar)
    : R3BHit(ID, x, y, ELoss, t)
    , fBarId(iBar)
{
}

Int_t R3BTofiHitData::GetBarId() const { return fBarId; }

ClassImp(R3BTofiHitData)
