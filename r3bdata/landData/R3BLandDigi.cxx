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

#include "R3BLandDigi.h"

#include <iostream>

using std::cout;
using std::endl;
using std::flush;

// -----   Default constructor   -------------------------------------------
R3BLandDigi::R3BLandDigi()
    : fPaddleNb(-1)
    , fTdcL(0.)
    , fTdcR(0.)
    , fTdc(0.)
    , fQdcL(0.)
    , fQdcR(0.)
    , fQdc(0.)
    , fxx(0.)
    , fyy(0.)
    , fzz(0.)
{
}

R3BLandDigi::R3BLandDigi(Int_t paddle,
                         Double_t tdcL,
                         Double_t tdcR,
                         Double_t tdc,
                         Double_t qdcL,
                         Double_t qdcR,
                         Double_t qdc,
                         Double_t xx,
                         Double_t yy,
                         Double_t zz)
    : fPaddleNb(paddle)
    , fTdcL(tdcL)
    , fTdcR(tdcR)
    , fTdc(tdc)
    , fQdcL(qdcL)
    , fQdcR(qdcR)
    , fQdc(qdc)
    , fxx(xx)
    , fyy(yy)
    , fzz(zz)
{
}

R3BLandDigi::R3BLandDigi(Int_t paddle, Double_t tdcR, Double_t tdcL, Double_t qdcR, Double_t qdcL)
    : fPaddleNb(paddle)
    , fTdcL(tdcL)
    , fTdcR(tdcR)
    , fTdc(0.)
    , fQdcL(qdcL)
    , fQdcR(qdcR)
    , fQdc(0.)
    , fxx(0.)
    , fyy(0.)
    , fzz(0.)
{
}

R3BLandDigi::R3BLandDigi(const R3BLandDigi& right)
    : fPaddleNb(right.fPaddleNb)
    , fTdcL(right.fTdcL)
    , fTdcR(right.fTdcR)
    , fTdc(right.fTdc)
    , fQdcL(right.fQdcL)
    , fQdcR(right.fQdcR)
    , fQdc(right.fQdc)
    , fxx(right.fxx)
    , fyy(right.fyy)
    , fzz(right.fzz)
{
}

// -----   Destructor   ----------------------------------------------------
R3BLandDigi::~R3BLandDigi() {}

// -----   Public method Print   -------------------------------------------
void R3BLandDigi::Print(const Option_t* opt) const
{
    cout << "-I- R3BLandDigi: LAND Digi in detector " << endl;
    cout << "    PaddleNr: " << fPaddleNb << endl;
    cout << "    TdcR: " << fTdcR << "    TdcL " << fTdcL << endl;
    cout << "    QdcR: " << fQdcR << "    QdcL " << fQdcL << endl;
}
// -------------------------------------------------------------------------

ClassImp(R3BLandDigi)
