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

//------------------------------------------------
// The Virtual Monte Carlo examples
// Copyright (C) 2007 - 2016 Ivana Hrivnacova
// All rights reserved.
//
// For the licensing terms see geant4_vmc/LICENSE.
// Contact: root-vmc@cern.ch
//-------------------------------------------------

/// \file  Garfield/geant4/src/GarfieldG4FastSimulationModel.cxx
/// \brief Implementation of the GarfieldG4FastSimulationModel class
///
/// Garfield garfieldpp example adapted to Virtual Monte Carlo.
/// This class is imported from garfieldpp example.
/// There is disabled WriteGeometryToGDML function in order to avoid
/// requiring Geant4 optional library.
///
/// \date 28/10/2015
/// \author D. Pheiffer, CERN

#include "GarfieldMessenger.h"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithoutParameter.hh"
#include "G4UIcommand.hh"
#include "G4UIdirectory.hh"
#include "G4UIparameter.hh"
#include "GarfieldPhysics.h"

GarfieldMessenger::GarfieldMessenger()
    : G4UImessenger()
    , fExampleDir(0)
    , fGarfieldPhysicsDir(0)
    , fIonizationModelCmd(0)
    , fGarfieldParticleCmd(0)
{
    fExampleDir = new G4UIdirectory("/exampleGarfield/");
    fExampleDir->SetGuidance("Commands specific to this example");

    G4bool broadcast = false;
    fGarfieldPhysicsDir = new G4UIdirectory("/exampleGarfield/physics/", broadcast);
    fGarfieldPhysicsDir->SetGuidance("Particle and energy ranges for Garfield++ physics model");

    fIonizationModelCmd = new G4UIcommand("/exampleGarfield/physics/setIonizationModel", this);
    fIonizationModelCmd->SetGuidance("Select ionization model for Garfield++");
    fIonizationModelCmd->SetGuidance("  and choose whether to use default particles");
    fIonizationModelCmd->SetGuidance("  and energy ranges for the chosen model");
    //
    G4UIparameter* ionizationModelPrm = new G4UIparameter("ionizationModel", 's', false);
    ionizationModelPrm->SetGuidance("ionization model (1. PAIPhot, 2. Heed)");
    ionizationModelPrm->SetGuidance("  1. Geant4 model, delta electrons transported by Heed");
    ionizationModelPrm->SetGuidance("  2. Use directly Heed");
    fIonizationModelCmd->SetParameter(ionizationModelPrm);
    //
    G4UIparameter* useDefaultsPrm = new G4UIparameter("useDefaults", 'b', false);
    useDefaultsPrm->SetGuidance("true to use default, false to manually choose particles and energies");
    fIonizationModelCmd->SetParameter(useDefaultsPrm);
    //
    fIonizationModelCmd->AvailableForStates(G4State_PreInit);

    fGarfieldParticleCmd = new G4UIcommand("/exampleGarfield/physics/setGarfieldParticleTypeAndEnergy", this);
    fGarfieldParticleCmd->SetGuidance("Select particle types and energies for Heed model.");
    fGarfieldParticleCmd->SetGuidance(" For PAI and PAIPhot model choose at which energy electrons are");
    fGarfieldParticleCmd->SetGuidance(" transported as delta electrons by Heed, and treatment of gammas");
    //
    G4UIparameter* particleGarfieldPrm = new G4UIparameter("particleName", 's', false);
    particleGarfieldPrm->SetGuidance(
        "Particle name (gamma, e-, e+, mu-, mu+, proton, anti_proton, pi-, pi+, kaon, kaon+, alpha, deuteron)");
    fGarfieldParticleCmd->SetParameter(particleGarfieldPrm);
    //
    G4UIparameter* minEnergyGarfieldPrm = new G4UIparameter("minimumEnergyGarfield", 'd', false);
    minEnergyGarfieldPrm->SetGuidance("minimum energy");
    minEnergyGarfieldPrm->SetParameterRange("minimumEnergyGarfield>=0");
    fGarfieldParticleCmd->SetParameter(minEnergyGarfieldPrm);
    //
    G4UIparameter* maxEnergyGarfieldPrm = new G4UIparameter("maximumEnergyGarfield", 'd', false);
    maxEnergyGarfieldPrm->SetGuidance("maximum energy");
    maxEnergyGarfieldPrm->SetParameterRange("maximumEnergyGarfield>=0");
    fGarfieldParticleCmd->SetParameter(maxEnergyGarfieldPrm);
    //
    G4UIparameter* unitGarfieldPrm = new G4UIparameter("unit", 's', false);
    unitGarfieldPrm->SetGuidance("unit of energy");
    G4String unitListGarfield = G4UIcommand::UnitsList(G4UIcommand::CategoryOf("MeV"));
    unitGarfieldPrm->SetParameterCandidates(unitListGarfield);
    fGarfieldParticleCmd->SetParameter(unitGarfieldPrm);
    //
    fGarfieldParticleCmd->AvailableForStates(G4State_PreInit);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

GarfieldMessenger::~GarfieldMessenger()
{
    delete fExampleDir;
    delete fGarfieldPhysicsDir;
    delete fIonizationModelCmd;
    delete fGarfieldParticleCmd;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void GarfieldMessenger::SetNewValue(G4UIcommand* command, G4String newValue)
{
    if (command == fIonizationModelCmd)
    {
        GarfieldPhysics* garfieldPhysics = GarfieldPhysics::GetInstance();
        G4String modelName;
        G4bool useDefaults;
        std::istringstream is(newValue);
        is >> modelName >> std::boolalpha >> useDefaults;
        garfieldPhysics->SetIonizationModel(modelName, useDefaults);
    }
    else if (command == fGarfieldParticleCmd)
    {
        GarfieldPhysics* garfieldPhysics = GarfieldPhysics::GetInstance();
        G4String particleName, unit, programName;
        G4double minEnergy;
        G4double maxEnergy;
        std::istringstream is(newValue);
        is >> particleName >> minEnergy >> maxEnergy >> unit;
        minEnergy *= G4UIcommand::ValueOf(unit);
        maxEnergy *= G4UIcommand::ValueOf(unit);
        garfieldPhysics->AddParticleName(particleName, minEnergy, maxEnergy);
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
