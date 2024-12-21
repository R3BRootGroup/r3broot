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

/// \ingroup Garfield
/// \file Garfield/g4Config.C
/// \brief Configuration macro for Geant4 VirtualMC for Garfield example
///
/// For geometry defined with Root and Geant4 native navigation

void Config()
{
    /// The configuration function for Geant4 VMC for Garfield example
    /// called during MC application initialization.
    /// For geometry defined with Root and Geant4 native navigation

    // cout << "Constructing Garfield::RunConfiguration " << endl;

    // Run configuration
    VMC::Garfield::RunConfiguration* runConfiguration = new VMC::Garfield::RunConfiguration("geomRoot", "FTFP_BERT");

    cout << "Constructing TG4RunConfiguration " << endl;

    // Run configuration with special cuts activated
    // VMC::Garfield::RunConfiguration* runConfiguration
    //   = new VMC::Garfield::RunConfiguration("geomRootToGeant4", "FTFP_BERT", "specialCuts");

    // TGeant4
    TGeant4* geant4 = new TGeant4("TGeant4", "The Geant4 Monte Carlo", runConfiguration);

    cout << "Geant4 has been created." << endl;

    // Customise Geant4 setting
    // Fast simulation model configuration
    // + verbose level, global range cuts, etc.
    // geant4->ProcessGeantMacro("g4garfieldconfig.in");

    /// create the Specific stack
    R3BStack* stack = new R3BStack(1000);
    stack->SetDebug(kFALSE);
    stack->StoreSecondaries(kTRUE);
    stack->SetMinPoints(0);
    geant4->SetStack(stack);

    if (FairRunSim::Instance()->IsExtDecayer())
    {
        TVirtualMCDecayer* decayer = TPythia6Decayer::Instance();
        geant4->SetExternalDecayer(decayer);
    }

    /// Customise Geant4 setting
    /// (verbose level, global range cut, ..)

    TString configm(gSystem->Getenv("VMCWORKDIR"));
    configm1 = configm + "/gconfig/g4garfieldconfig.in";
    cout << " -I g4GarfieldConfig() using g4garfieldconf  macro: " << configm1 << endl;

    // set geant4 specific stuff
    geant4->SetMaxNStep(10000); // default is 30000
    geant4->ProcessGeantMacro(configm1.Data());

    cout << "Processing Config() done." << endl;
}
