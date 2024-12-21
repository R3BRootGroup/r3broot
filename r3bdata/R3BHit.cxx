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

// -----------------------------------------------------------------------------
// -----                              R3BHit                               -----
// -----                 Created on 01.03.2018 by D.Kresan                 -----
// -----------------------------------------------------------------------------

#include "R3BHit.h"

R3BHit::R3BHit()
    : fDetId(0)
    , fX(0.)
    , fY(0.)
    , fEloss(0.)
    , fTime(0.)
    , fHitId(-1)
{
}

R3BHit::R3BHit(Int_t detId, Double_t x, Double_t y, Double_t eloss, Double_t time, Int_t hitId)
    : fDetId(detId)
    , fX(x)
    , fY(y)
    , fEloss(eloss)
    , fTime(time)
    , fHitId(hitId)
{
}

R3BHit::~R3BHit() {}

ClassImp(R3BHit)
