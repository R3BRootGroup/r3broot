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

#include <iostream>

#include "R3BPtofHit.h"

using std::cout;
using std::endl;

R3BPtofHit::R3BPtofHit()
    : fPaddleNb(-1)
    , fCharge(0)
    , fTdc(0.)
    , fxx(0.)
    , fyy(0.)
    , fzz(0.)
{
}

R3BPtofHit::R3BPtofHit(Int_t paddle, Double_t charge, Double_t tdc, Double_t xx, Double_t yy, Double_t zz)
    : fPaddleNb(paddle)
    , fCharge(charge)
    , fTdc(tdc)
    , fxx(xx)
    , fyy(yy)
    , fzz(zz)
{
}

R3BPtofHit::R3BPtofHit(const R3BPtofHit& digi)
    : fPaddleNb(digi.fPaddleNb)
    , fCharge(digi.fCharge)
    , fTdc(digi.fTdc)
    , fxx(digi.fxx)
    , fyy(digi.fyy)
    , fzz(digi.fzz)
{
}

R3BPtofHit::~R3BPtofHit() {}

void R3BPtofHit::Print(const Option_t* opt) const
{
    cout << "-I- R3BPtofHit: Ptof Digi in detector " << endl;
    cout << "    PaddleNr: " << fPaddleNb << endl;
    cout << "    Charge: " << fCharge << endl;
    cout << "    Time: " << fTdc << endl;
    cout << "    Position:   x: " << fxx << "   y: " << fyy << "   z: " << fzz << endl;
}

ClassImp(R3BPtofHit)
