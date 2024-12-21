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
// -----                      R3BLandPoint source file                  -----
// -------------------------------------------------------------------------

#include "R3BLandFirstHits.h"

#include <iostream>

using std::cout;
using std::endl;
using std::flush;

// -----   Default constructor   -------------------------------------------
R3BLandFirstHits::R3BLandFirstHits()
    : fx0(0.)
    , fy0(0.)
    , fz0(0.)
    , fT0(0.)
    , fx1(0.)
    , fy1(0.)
    , fz1(0.)
    , fT1(0.)
    , fx2(0.)
    , fy2(0.)
    , fz2(0.)
    , fT2(0.)
    , fx3(0.)
    , fy3(0.)
    , fz3(0.)
    , fT3(0.)
    , fx4(0.)
    , fy4(0.)
    , fz4(0.)
    , fT4(0.)
    , fx5(0.)
    , fy5(0.)
    , fz5(0.)
    , fT5(0.)
{
}

R3BLandFirstHits::R3BLandFirstHits(Double_t x0,
                                   Double_t y0,
                                   Double_t z0,
                                   Double_t T0,
                                   Double_t x1,
                                   Double_t y1,
                                   Double_t z1,
                                   Double_t T1,
                                   Double_t x2,
                                   Double_t y2,
                                   Double_t z2,
                                   Double_t T2,
                                   Double_t x3,
                                   Double_t y3,
                                   Double_t z3,
                                   Double_t T3,
                                   Double_t x4,
                                   Double_t y4,
                                   Double_t z4,
                                   Double_t T4,
                                   Double_t x5,
                                   Double_t y5,
                                   Double_t z5,
                                   Double_t T5)
    : fx0(x0)
    , fy0(y0)
    , fz0(z0)
    , fT0(T0)
    , fx1(x1)
    , fy1(y1)
    , fz1(z1)
    , fT1(T1)
    , fx2(x2)
    , fy2(y2)
    , fz2(z2)
    , fT2(T2)
    , fx3(x3)
    , fy3(y3)
    , fz3(z3)
    , fT3(T3)
    , fx4(x4)
    , fy4(y4)
    , fz4(z4)
    , fT4(T4)
    , fx5(x5)
    , fy5(y5)
    , fz5(z5)
    , fT5(T5)
{
}

R3BLandFirstHits::R3BLandFirstHits(const R3BLandFirstHits& right)
    : fx0(right.fx0)
    , fy0(right.fy0)
    , fz0(right.fz0)
    , fT0(right.fT0)
    , fx1(right.fx1)
    , fy1(right.fy1)
    , fz1(right.fz1)
    , fT1(right.fT1)
    , fx2(right.fx2)
    , fy2(right.fy2)
    , fz2(right.fz2)
    , fT2(right.fT2)
    , fx3(right.fx3)
    , fy3(right.fy3)
    , fz3(right.fz3)
    , fT3(right.fT3)
    , fx4(right.fx4)
    , fy4(right.fy4)
    , fz4(right.fz4)
    , fT4(right.fT4)
    , fx5(right.fx5)
    , fy5(right.fy5)
    , fz5(right.fz5)
    , fT5(right.fT5)
{
}

// -----   Destructor   ----------------------------------------------------
R3BLandFirstHits::~R3BLandFirstHits() {}

// -----   Public method Print   -------------------------------------------
void R3BLandFirstHits::Print(const Option_t* opt) const
{
    cout << "-I- R3BLandFirstHits: " << endl;
    cout << "x0: " << fx0 << " y0: " << fy0 << " z0: " << fz0 << " T0: " << fT0 << endl;
    cout << "x1: " << fx1 << " y1: " << fy1 << " z1: " << fz1 << " T1: " << fT1 << endl;
    cout << "x2: " << fx2 << " y2: " << fy2 << " z2: " << fz2 << " T2: " << fT2 << endl;
    cout << "x3: " << fx3 << " y3: " << fy3 << " z3: " << fz3 << " T3: " << fT3 << endl;
    cout << "x4: " << fx4 << " y4: " << fy4 << " z4: " << fz4 << " T4: " << fT4 << endl;
    cout << "x5: " << fx5 << " y5: " << fy5 << " z5: " << fz5 << " T5: " << fT5 << endl;
}
// -------------------------------------------------------------------------

ClassImp(R3BLandFirstHits)
