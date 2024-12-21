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
// -----                      R3BSitrackerHit source file                  -----
// -------------------------------------------------------------------------

#include "R3BSitrackerHit.h"

#include <iostream>

using std::cout;
using std::endl;
using std::flush;

// -----   Default constructor   -------------------------------------------
R3BSitrackerHit::R3BSitrackerHit()
    : FairHit()
{
    fX = fY = fZ = 0.;
    type = 0;
    copy = 0;
    fEnergy = 0.;
    nbOfSteps = 0;
    eventID = 0;
    runID = 0;

    fTimeFirstStep = 0.;
    fTimeLastStep = 0.;

    nbOfPrimaries = 0;
    fEnergyPrimary = 0.;
    fThetaPrimary = 0.;
    fPhiPrimary = 0.;
}

// -----   Destructor   ----------------------------------------------------
R3BSitrackerHit::~R3BSitrackerHit() {}
// -------------------------------------------------------------------------

// -----   Public method Print   -------------------------------------------
void R3BSitrackerHit::Print(const Option_t* opt) const
{
    cout << "==============================================" << endl;
    cout << "-I- R3BSitrackerHit:  Hit in detector: " << fDetectorID << endl;
    cout << "    Position (" << fX << ", " << fY << ", " << fZ << ") cm" << endl;

    cout << " type: " << type << ", copy: " << copy << endl;
    cout << " energy: " << fEnergy * 1000. << " MeV"
         << ", nbOfSteps: " << nbOfSteps << endl;
    cout << " nbOfPrimaries: " << nbOfPrimaries << ", energyPrimary: " << fEnergyPrimary * 1000. << " MeV"
         << ", thetaPrimary: " << fThetaPrimary << " deg"
         << ", phiPrimary: " << fPhiPrimary << " deg" << endl;

    cout << "----------------------------------------------" << endl;
}
// -------------------------------------------------------------------------

ClassImp(R3BSitrackerHit)
