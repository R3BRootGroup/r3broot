#include "FairBoxGenerator.h"
#include "FairParRootFileIo.h"
#include "FairPrimaryGenerator.h"
#include "FairRootFileSink.h"
#include "FairRunSim.h"
#include "R3BCave.h"
#include "R3BNeuland.h"
#include "TStopwatch.h"
#include <G4RunManager.hh>
#include <G4UserEventAction.hh>
#include <R3BFieldConst.h>
#include <R3BPhaseSpaceGenerator.h>
#include <R3BProgramOptions.h>
#include <TG4EventAction.h>
#include <boost/exception/diagnostic_information.hpp>
#include <boost/program_options.hpp>
#include <iostream>
#include <string>

constexpr int DEFAULT_RUNID = 999;

int main(int argc, const char** argv)
{
    auto timer = TStopwatch{};
    auto const defaultEventNum = 10;
    timer.Start();

    auto programOptions = R3B::ProgramOptions("options for neuland simulation");

    auto help = programOptions.Create_Option<bool>("help,h", "help message", false);
    auto eventNum = programOptions.Create_Option<int>("eventNum", "set total event number", defaultEventNum);
    auto eventPrintNum = programOptions.Create_Option<int>("eventPrint", "set event print number", 1);
    auto runID = programOptions.Create_Option<int>("runID", "set runID", DEFAULT_RUNID);
    auto multi = programOptions.Create_Option<int>("multiplicity", "set particle multiplicity", 1);
    auto pEnergy = programOptions.Create_Option<double>("energy", "set energy value (GeV) of the particle", 1);
    auto simuFileName =
        programOptions.Create_Option<std::string>("simuFile", "set the base filename of simulation ouput", "simu.root");
    auto paraFileName =
        programOptions.Create_Option<std::string>("paraFile", "set the base filename of parameter sink", "para.root");
    auto logLevel = programOptions.Create_Option<std::string>("logLevel,v", "set log level of fairlog", "error");

    if (!programOptions.Verify(argc, argv))
    {
        return EXIT_FAILURE;
    }

    if (help->value())
    {
        std::cout << programOptions.Get_DescRef() << std::endl;
        return 0;
    }

    // Logging
    // FairLogger::GetLogger()->SetLogVerbosityLevel("LOW");
    FairLogger::GetLogger()->SetLogScreenLevel(logLevel->value().c_str());

    // System paths
    const TString workDirectory = getenv("VMCWORKDIR");
    gSystem->Setenv("GEOMPATH", workDirectory + "/geometry");
    gSystem->Setenv("CONFIG_DIR", workDirectory + "/gconfig");

    // Basic simulation setup
    auto run = std::make_unique<FairRunSim>();
    run->SetName("TGeant4");
    run->SetRunId(runID->value());
    run->SetStoreTraj(false);
    run->SetMaterials("media_r3b.geo");
    run->SetSink(std::make_unique<FairRootFileSink>(simuFileName->value().c_str()));
    auto fairField = std::make_unique<R3BFieldConst>();
    run->SetField(fairField.release());

    // Primary particle generator
    auto gen = std::make_unique<R3BPhaseSpaceGenerator>();

    constexpr auto beam_energy = 883.;      // MeV
    constexpr auto rel_energy_max = 10000.; // keV
    constexpr auto Sn_atomic_number = 50;
    constexpr auto Sn_mass = 123;
    constexpr auto neutron_pdg = 2112;
    gen->GetBeam().SetEnergyDistribution(R3BDistribution1D::Delta(beam_energy));
    gen->SetErelDistribution(R3BDistribution1D::Flat(0., rel_energy_max));
    gen->AddParticle(Sn_atomic_number, Sn_mass);
    gen->AddParticle(neutron_pdg);
    gen->EnableWhitelist();
    gen->EnableWrite();
    gen->AddParticleToWhitelist(neutron_pdg);
    auto primGen = std::make_unique<FairPrimaryGenerator>();
    primGen->AddGenerator(gen.release());

    // auto boxGen = std::make_unique<FairBoxGenerator>(PID, multi->value());
    // boxGen->SetXYZ(0, 0, 0.);
    // boxGen->SetThetaRange(0., 3.);
    // boxGen->SetPhiRange(0., 360.);
    // boxGen->SetEkinRange(pEnergy->value(), pEnergy->value());
    // primGen->AddGenerator(boxGen.release());
    run->SetGenerator(primGen.release());

    // Geometry: Cave
    auto cave = std::make_unique<R3BCave>("CAVE");
    cave->SetGeometryFileName("r3b_cave.geo");
    run->AddModule(cave.release());

    // Geometry: Neuland
    auto const nDP = 13;
    auto const neulandGeoTrans = TGeoTranslation{ 0., 0., 1650. };
    auto neuland = std::make_unique<R3BNeuland>(nDP, neulandGeoTrans);
    run->AddModule(neuland.release());

    // Init
    run->Init();

    // event print out:
    auto* grun = G4RunManager::GetRunManager();
    grun->SetPrintProgress(eventPrintNum->value());
    auto* event = dynamic_cast<TG4EventAction*>(const_cast<G4UserEventAction*>(grun->GetUserEventAction())); // NOLINT
    event->VerboseLevel(0);

    // Connect runtime parameter file
    auto parFileIO = std::make_unique<FairParRootFileIo>(true);
    parFileIO->open(paraFileName->value().c_str());
    auto* rtdb = run->GetRuntimeDb();
    rtdb->setOutput(parFileIO.release());
    rtdb->saveOutput();

    // Simulate
    run->Run(eventNum->value());

    // Report
    timer.Stop();
    std::cout << "Macro finished successfully." << std::endl;
    std::cout << "Real time: " << timer.RealTime() << "s, CPU time: " << timer.CpuTime() << "s" << std::endl;
    return 0;
}
