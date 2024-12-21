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
// -----                      R3BTraHit source file                  -----
// -------------------------------------------------------------------------

#include "R3BTraHit.h"

#include <iostream>

using std::cout;
using std::endl;
using std::flush;

// -----   Default constructor   -------------------------------------------
R3BTraHit::R3BTraHit()
    : FairHit()
{
    fEnergy = fTof = 0.;
    fDetectorID = -1;
}

// -----   Destructor   ----------------------------------------------------
R3BTraHit::~R3BTraHit() {}

// -----   Public method Print   -------------------------------------------
void R3BTraHit::Print(const Option_t* opt) const
{
    cout << "-I- R3BTraHit: Tracker Hit in detector " << fDetectorID << endl;
    cout << "    Tof: " << fTof << " cm,  Energy loss " << fEnergy * 1.0e06 << " keV" << endl;
}
// -------------------------------------------------------------------------

ClassImp(R3BTraHit)
