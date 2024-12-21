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

#ifndef GARFIELD_MESSENGER_H
#define GARFIELD_MESSENGER_H

//------------------------------------------------
// The Virtual Monte Carlo examples
// Copyright (C) 2007 - 2016 Ivana Hrivnacova
// All rights reserved.
//
// For the licensing terms see geant4_vmc/LICENSE.
// Contact: root-vmc@cern.ch
//-------------------------------------------------

/// \file  Garfield/geant4/include/GarfieldMessenger.h
/// \brief Definition of the GarfieldMessengerclass
///
/// Garfield garfieldpp example adapted to Virtual Monte Carlo.
/// This class is imported from garfieldpp example.
/// Kept only commands relevant to VMC version of the example.
///
/// \date 28/10/2015
/// \author D. Pheiffer, CERN

#include "G4UImessenger.hh"
#include "globals.hh"

class G4UIdirectory;
class G4UIcommand;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class GarfieldMessenger : public G4UImessenger
{
  public:
    GarfieldMessenger();
    ~GarfieldMessenger();

    virtual void SetNewValue(G4UIcommand*, G4String);

  private:
    G4UIdirectory* fExampleDir;
    G4UIdirectory* fGarfieldPhysicsDir;
    G4UIcommand* fIonizationModelCmd;
    G4UIcommand* fGarfieldParticleCmd;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif // GARFIELDMESSENGER_H
