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
// -----                      R3BTrackerHit source file                  -----
// -------------------------------------------------------------------------

#include "R3BTrackerHit.h"

#include <iostream>

using std::cout;
using std::endl;
using std::flush;

// -----   Default constructor   -------------------------------------------
R3BTrackerHit::R3BTrackerHit()
    : FairMultiLinkedData()
{
    // fNbOfTrackerHits = 0;
    fEnergy = -1.;
    fDetector = 0;
    fXHit = 0.;
    fYHit = 0.;
    fZHit = 0.;
    fPx = 0.;
    fPy = 0.;
    fPz = 0.;
    fThetaZero = 0.;
    fPhiZero = 0.;
}
// -------------------------------------------------------------------------

// -----   Standard constructor   ------------------------------------------
R3BTrackerHit::R3BTrackerHit(Double_t ene,
                             Int_t det,
                             Double_t x,
                             Double_t y,
                             Double_t z,
                             Double_t px,
                             Double_t py,
                             Double_t pz,
                             Double_t th,
                             Double_t phi)
    // R3BTrackerHit::R3BTrackerHit(Double_t ene, Int_t det)
    : FairMultiLinkedData()
{

    // fNbOfTrackerHits = Nb;
    fEnergy = ene;
    fDetector = det;
    fXHit = x;
    fYHit = y;
    fZHit = z;
    fPx = px;
    fPy = py;
    fPz = pz;
    fThetaZero = th;
    fPhiZero = phi;
}
// -------------------------------------------------------------------------

// -----   Destructor   ----------------------------------------------------
R3BTrackerHit::~R3BTrackerHit() {}
// -------------------------------------------------------------------------

// -----   Public method Print   -------------------------------------------
// void R3BTrackerHit::Print(const Option_t* opt) const {
// cout << "-I- R3BTrackerHit: a Si tracker Hit made of " << fNbOfTrackerHits
//<< " Energy = " << fEnergy << " GeV" << endl;
//}
// -------------------------------------------------------------------------
