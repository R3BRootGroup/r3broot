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

// ------------------------------------------------------------------------
// -----                        R3BNeutHit                            -----
// -----              Created 21.05.2012 by D.Kresan                  -----
// ------------------------------------------------------------------------

#include "R3BNeutHit.h"

// ------------------------------------------------------------------------
R3BNeutHit::R3BNeutHit()
    : fX(0.)
    , fY(0.)
    , fZ(0.)
    , fT(0.)
{
}
// ------------------------------------------------------------------------

// ------------------------------------------------------------------------
R3BNeutHit::R3BNeutHit(Double_t x, Double_t y, Double_t z, Double_t t)
{
    fX = x;
    fY = y;
    fZ = z;
    fT = t;
}
// ------------------------------------------------------------------------

// ------------------------------------------------------------------------
R3BNeutHit::~R3BNeutHit() {}
// ------------------------------------------------------------------------

// ------------------------------------------------------------------------
Double_t R3BNeutHit::GetP() const
{
    Double_t beta = GetD() / GetT() / 29.9792458;
    Double_t gamma = 1. / sqrt(1. - beta * beta);
    Double_t p = beta * gamma * 1.0086649156 * 931.494028;
    return p;
}
// ------------------------------------------------------------------------

ClassImp(R3BNeutHit)
