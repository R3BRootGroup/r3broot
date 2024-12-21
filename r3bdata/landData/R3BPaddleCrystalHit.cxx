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
// -----                      R3BPaddleCrystalHit source file                  -----
// -------------------------------------------------------------------------

#include "R3BPaddleCrystalHit.h"

#include <iostream>

using std::cout;
using std::endl;
using std::flush;

// -----   Default constructor   -------------------------------------------
R3BPaddleCrystalHit::R3BPaddleCrystalHit()
    : FairHit()
    , type(0)
    , copy(0)
    , fEnergy(0.)
    , nbOfSteps(0)
    , fTimeFirstStep(0.)
    , fTimeLastStep(0.)
    , nbOfPrimaries(0)
    , fEnergyPrimary(0.)
    , fThetaPrimary(0.)
    , fPhiPrimary(0.)
    , eventID(0)
    , runID(0)
{
}

R3BPaddleCrystalHit::R3BPaddleCrystalHit(const R3BPaddleCrystalHit& right)
    : FairHit(right)
    , type(right.type)
    , copy(right.copy)
    , fEnergy(right.fEnergy)
    , nbOfSteps(right.nbOfSteps)
    , fTimeFirstStep(right.fTimeFirstStep)
    , fTimeLastStep(right.fTimeLastStep)
    , nbOfPrimaries(right.nbOfPrimaries)
    , fEnergyPrimary(right.fEnergyPrimary)
    , fThetaPrimary(right.fThetaPrimary)
    , fPhiPrimary(right.fPhiPrimary)
    , eventID(right.eventID)
    , runID(right.runID)
{
}

// -----   Destructor   ----------------------------------------------------
R3BPaddleCrystalHit::~R3BPaddleCrystalHit() {}
// -------------------------------------------------------------------------

// -----   Public method Print   -------------------------------------------
void R3BPaddleCrystalHit::Print(const Option_t* opt) const
{
    cout << "==============================================" << endl;
    cout << "-I- R3BPaddleCrystalHit:  Hit in detector: " << fDetectorID << endl;
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

ClassImp(R3BPaddleCrystalHit)
