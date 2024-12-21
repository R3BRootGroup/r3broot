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
// -----                      R3BPspPoint source file                  -----
// -------------------------------------------------------------------------

#include "R3BPspDigi.h"

#include <iostream>

using std::cout;
using std::endl;
using std::flush;

// -----   Default constructor   -------------------------------------------
R3BPspDigi::R3BPspDigi()
    : Ps03mul(0)
{
}

R3BPspDigi::R3BPspDigi(Int_t psp3mul, Double_t psp3x, Double_t psp3y, Double_t psp3e)
    : R3BHit(0, psp3x, psp3y, psp3e, 0.)
    , Ps03mul(psp3mul)
{
}

// -----   Destructor   ----------------------------------------------------
R3BPspDigi::~R3BPspDigi() {}

// -----   Public method Print   -------------------------------------------
void R3BPspDigi::Print(const Option_t* opt) const {}
// -------------------------------------------------------------------------

ClassImp(R3BPspDigi)
