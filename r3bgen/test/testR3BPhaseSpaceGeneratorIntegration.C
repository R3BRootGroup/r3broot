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

#include <TStopwatch.h>
#include <TString.h>
#include <TSystem.h>
#include <iostream>
#include <memory>

void testR3BPhaseSpaceGeneratorIntegration()
{
    // Timer
    TStopwatch timer;
    timer.Start();

    // Logging
    auto logger = FairLogger::GetLogger();
    logger->SetLogVerbosityLevel("low");
    logger->SetLogScreenLevel("warn");
    logger->SetColoredLog(true);

    // System paths
    const TString workDirectory = getenv("VMCWORKDIR");
    gSystem->Setenv("GEOMPATH", workDirectory + "/geometry");
    gSystem->Setenv("CONFIG_DIR", workDirectory + "/gconfig");

    // Simulation Base
    FairRunSim run;
    run.SetName("TGeant4");
    run.SetOutputFile("testR3BPhaseSpaceGeneratorIntegration.root");
    run.SetMaterials("media_r3b.geo");

    // World
    auto cave = new R3BCave("Cave");
    cave->SetGeometryFileName("r3b_cave_vacuum.geo");
    run.AddModule(cave);

    // Magnet
    // run.AddModule(new R3BGladMagnet("glad_v2023.1.geo.root"));
    auto magField = new R3BGladFieldMap("R3BGladMap");
    magField->SetScale(-0.6);
    run.SetField(magField);

    // Primaries
    auto primGen = new FairPrimaryGenerator();
    auto gen = new R3BPhaseSpaceGenerator();
    gen->GetBeam().SetEnergyDistribution(R3BDistribution1D::Delta(600));
    gen->SetErelDistribution(R3BDistribution1D::Delta(100));
    gen->AddParticle(5, 17);
    gen->AddNeutron();
    gen->AddProton();
    primGen->AddGenerator(gen);
    run.SetGenerator(primGen);

    // Logging
    run.SetStoreTraj(false);
    FairLogger::GetLogger()->SetLogVerbosityLevel("low");
    FairLogger::GetLogger()->SetLogScreenLevel("info");

    // Init & Special MC Settings
    run.Init();

    // Database
    auto rtdb = run.GetRuntimeDb();
    auto fieldPar = (R3BFieldPar*)rtdb->getContainer("R3BFieldPar");
    fieldPar->SetParameters(magField);
    fieldPar->setChanged();
    auto parOut = new FairParRootFileIo(true);
    parOut->open("testR3BPhaseSpaceGeneratorIntegration.para.root");
    rtdb->setOutput(parOut);
    rtdb->saveOutput();
    rtdb->print();

    // Go
    run.Run(10);

    // Test Output
    auto file = TFile::Open("testR3BPhaseSpaceGeneratorIntegration.root");
    auto tree = (TTree*)file->Get("evt");
    auto mctc = new TClonesArray("R3BMCTrack");
    tree->SetBranchAddress("MCTrack", &mctc);

    tree->GetEvent(0);
    if (mctc->GetEntries() < 3)
    {
        std::cout << "Not enough particles produced" << std::endl;
        return;
    }

    auto track = (R3BMCTrack*)mctc->At(1);
    if (track->GetPdgCode() != 2112 || track->GetMotherId() != -1)
    {
        std::cout << "Not the correct primary particle" << std::endl;
        return;
    }

    const Double_t ekin = track->GetEnergy() - track->GetMass();
    std::cout << "Ekin of primary neutron:" << ekin << std::endl;
    if (abs(ekin - 0.6) > 0.02)
    {
        std::cout << "Primary neutron doesn't have the correct energy!" << std::endl;
        return;
    }

    // Report
    timer.Stop();
    std::cout << "Real time: " << timer.RealTime() << "s, CPU time: " << timer.CpuTime() << "s" << std::endl;
    std::cout << "Macro finished successfully." << std::endl;
}
