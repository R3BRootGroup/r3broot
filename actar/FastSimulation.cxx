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

/// \file Garfield/geant4/src/FastSimulation.cxx
/// \brief Implementation of the Garfield::FastSimulation class
///
/// Garfield garfieldpp example adapted to Virtual Monte Carlo.
///
/// \date 28/10/2015
/// \author I. Hrivnacova; IPN, Orsay

#include "FastSimulation.h"
#include "GarfieldG4FastSimulationModel.h"
#include "GarfieldMessenger.h"
// #include "GarfieldPhysics.h"

#include <Random.hh>

#include <G4FastSimulationManager.hh>
#include <G4Material.hh>
#include <G4NistManager.hh>
#include <G4RegionStore.hh>
#include <GFlashHitMaker.hh>
#include <GFlashHomoShowerParameterisation.hh>
#include <GFlashParticleBounds.hh>
#include <GFlashShowerModel.hh>
#include <Randomize.hh>

#include <Riostream.h>

using namespace std;

namespace VMC
{
    namespace Garfield
    {

        //_____________________________________________________________________________
        FastSimulation::FastSimulation()
            : TG4VUserFastSimulation()
            , fMessenger(0)
        {
            /// Standard constructor

            cout << "FastSimulation::FastSimulation" << endl;

            // Choose the Random engine
            G4Random::setTheEngine(new CLHEP::RanecuEngine);
            G4Random::setTheSeed(1);
            // Set seed to Garfield random engine
            ::Garfield::randomEngine.Seed(1);

            // Construct the Garfield messenger which defines the Garfield physics specific
            // command
            fMessenger = new GarfieldMessenger();

            /*
              // In the following calls users can select the particles and regions
              // which the fast simulation model(s) will be applied to.
              // The setting is an alternative to the setting via UI commands
              // in physics.in macro.

              // Create fast simulation model configuration.
              // This will generate UI commands which can be used to set particles
              // and regions where the model will be applied
              SetModel("garfieldModel");

              // In the following calls users can select the particles and regions
              // which the fast simulation model(s) will be applied to.
              // The setting can be done also interactively via UI commands.
              SetModelParticles("garfieldModel", "all");
              SetModelRegions("garfieldModel", "AirB");

              // Enable GarfieldModel for different particle types and energy ranges
              double minEnergy_keV = 100;
              double maxEnergy_keV = 1e+12;

              //GarfieldPhysics* garfieldPhysics = GarfieldPhysics::GetInstance();
              //garfieldPhysics->AddParticleName("e-", minEnergy_keV, maxEnergy_keV);
              //garfieldPhysics->AddParticleName("e+", minEnergy_keV, maxEnergy_keV);

              //garfieldPhysics->AddParticleName("mu-", minEnergy_keV, maxEnergy_keV);
              garfieldPhysics->AddParticleName("mu+", minEnergy_keV, maxEnergy_keV);

              //garfieldPhysics->AddParticleName("pi-", minEnergy_keV, maxEnergy_keV);
              //garfieldPhysics->AddParticleName("pi+", minEnergy_keV, maxEnergy_keV);
              //garfieldPhysics->AddParticleName("kaon-", minEnergy_keV, maxEnergy_keV);
              //garfieldPhysics->AddParticleName("kaon+", minEnergy_keV, maxEnergy_keV);
              //garfieldPhysics->AddParticleName("proton", minEnergy_keV, maxEnergy_keV);
              //garfieldPhysics->AddParticleName("anti_proton", minEnergy_keV, maxEnergy_keV);
              //garfieldPhysics->AddParticleName("deuteron", minEnergy_keV, maxEnergy_keV);
              //garfieldPhysics->AddParticleName("alpha", minEnergy_keV, maxEnergy_keV);

              //garfieldPhysics->EnableCreateSecondariesInGeant4(false);
            */
        }

        //_____________________________________________________________________________
        FastSimulation::~FastSimulation()
        {
            /// Destructor

            // can't we just delete the instance ??
            GarfieldPhysics::Dispose();

            delete fMessenger;
        }

        //
        // protected methods
        //

        //_____________________________________________________________________________
        void FastSimulation::Construct()
        {
            /// This function must be overriden in user class and users should create
            /// the simulation models and register them to VMC framework

            G4cout << "Construct Garfield model." << G4endl;

            // Create the fast simulation model
            GarfieldG4FastSimulationModel* garfieldModel = new GarfieldG4FastSimulationModel("garfieldModel");

            // Register the model in the VMC framework
            Register(garfieldModel);

            G4cout << "end construct Garfield model." << G4endl;
            //
            // end Initializing shower model
        }

    } // namespace Garfield
} // namespace VMC
