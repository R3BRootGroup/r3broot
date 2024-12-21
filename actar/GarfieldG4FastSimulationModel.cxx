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

/*
 * GarfieldModel.cpp
 *
 *  Created on: Apr 9, 2014
 *      Author: dpfeiffe
 */
#include "GarfieldG4FastSimulationModel.h"
#include "G4Electron.hh"
#include "G4Gamma.hh"
#include "G4VPhysicalVolume.hh"
#include <iostream>

#include "G4SystemOfUnits.hh"

// I.H. make this optional
// #include "G4GDMLParser.hh"

GarfieldG4FastSimulationModel::GarfieldG4FastSimulationModel(G4String modelName, G4Region* envelope)
    : G4VFastSimulationModel(modelName, envelope)
{
    fGarfieldPhysics = GarfieldPhysics::GetInstance();
    fGarfieldPhysics->InitializePhysics();
}

GarfieldG4FastSimulationModel::GarfieldG4FastSimulationModel(G4String modelName)
    : G4VFastSimulationModel(modelName)
{
    // I.H.
    fGarfieldPhysics = GarfieldPhysics::GetInstance();
    fGarfieldPhysics->InitializePhysics();
}

GarfieldG4FastSimulationModel::~GarfieldG4FastSimulationModel() {}

// I.H. make this optional
// void GarfieldG4FastSimulationModel::WriteGeometryToGDML(
// 		G4VPhysicalVolume* physicalVolume) {

// 	G4GDMLParser* parser = new G4GDMLParser();
// 	remove("garfieldGeometry.gdml");
// 	parser->Write("garfieldGeometry.gdml", physicalVolume, false);
// 	delete parser;
// }

G4bool GarfieldG4FastSimulationModel::IsApplicable(const G4ParticleDefinition& particleType)
{

    G4String particleName = particleType.GetParticleName();

    if (fGarfieldPhysics->FindParticleName(particleName))
    {
        return true;
    }
    return false;
}

G4bool GarfieldG4FastSimulationModel::ModelTrigger(const G4FastTrack& fastTrack)
{
    double ekin_MeV = fastTrack.GetPrimaryTrack()->GetKineticEnergy() / MeV;
    G4String particleName = fastTrack.GetPrimaryTrack()->GetParticleDefinition()->GetParticleName();
    if (fGarfieldPhysics->FindParticleNameEnergy(particleName, ekin_MeV))
    {
        return true;
    }
    return false;
}

void GarfieldG4FastSimulationModel::DoIt(const G4FastTrack& fastTrack, G4FastStep& fastStep)
{

    G4TouchableHandle theTouchable = fastTrack.GetPrimaryTrack()->GetTouchableHandle();
    G4String name = theTouchable->GetVolume()->GetName();

    G4ThreeVector pdirection = fastTrack.GetPrimaryTrack()->GetMomentum().unit();
    G4ThreeVector localdir = fastTrack.GetPrimaryTrackLocalDirection();

    G4ThreeVector worldPosition = fastTrack.GetPrimaryTrack()->GetPosition();
    G4ThreeVector localPosition = fastTrack.GetPrimaryTrackLocalPosition();

    double ekin_MeV = fastTrack.GetPrimaryTrack()->GetKineticEnergy() / MeV;
    G4double globalTime = fastTrack.GetPrimaryTrack()->GetGlobalTime();

    G4String particleName = fastTrack.GetPrimaryTrack()->GetParticleDefinition()->GetParticleName();

    fastStep.KillPrimaryTrack();
    fastStep.SetPrimaryTrackPathLength(0.0);

    if (particleName == "kaon+")
    {
        particleName = "K+";
    }
    else if (particleName == "kaon-")
    {
        particleName = "K-";
    }
    else if (particleName == "anti_proton")
    {
        particleName = "anti-proton";
    }

    fGarfieldPhysics->DoIt(particleName,
                           ekin_MeV,
                           globalTime,
                           localPosition.x() / CLHEP::cm,
                           localPosition.y() / CLHEP::cm,
                           localPosition.z() / CLHEP::cm,
                           localdir.x(),
                           localdir.y(),
                           localdir.z());

    fastStep.SetTotalEnergyDeposited(fGarfieldPhysics->GetEnergyDeposit_MeV());

    if (fGarfieldPhysics->GetCreateSecondariesInGeant4())
    {
        std::vector<GarfieldParticle*>* secondaryParticles = fGarfieldPhysics->GetSecondaryParticles();

        if (secondaryParticles->size() > 0)
        {
            fastStep.SetNumberOfSecondaryTracks(secondaryParticles->size());

            G4double totalEnergySecondaries_MeV = 0;

            for (std::vector<GarfieldParticle*>::iterator it = secondaryParticles->begin();
                 it != secondaryParticles->end();
                 ++it)
            {
                G4double x = (*it)->getX_mm();
                G4double y = (*it)->getY_mm();
                G4double z = (*it)->getZ_mm();
                G4double eKin_MeV = (*it)->getEkin_MeV();
                G4double dx = (*it)->getDX();
                G4double dy = (*it)->getDY();
                G4double dz = (*it)->getDZ();
                G4double time = (*it)->getTime();
                G4ThreeVector momentumDirection(dx, dy, dz);
                G4ThreeVector position(x, y, z);
                if ((*it)->getParticleName() == "e-")
                {
                    G4DynamicParticle particle(G4Electron::ElectronDefinition(), momentumDirection, eKin_MeV);
                    fastStep.CreateSecondaryTrack(particle, position, time, true);
                    totalEnergySecondaries_MeV += eKin_MeV;
                }
                else if ((*it)->getParticleName() == "gamma")
                {
                    G4DynamicParticle particle(G4Gamma::GammaDefinition(), momentumDirection, eKin_MeV);
                    fastStep.CreateSecondaryTrack(particle, position, time, true);
                    totalEnergySecondaries_MeV += eKin_MeV;
                }
            }
        }
    }
}
