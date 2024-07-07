/******************************************************************************
 *   Copyright (C) 2019 GSI Helmholtzzentrum für Schwerionenforschung GmbH    *
 *   Copyright (C) 2019-2023 Members of R3B Collaboration                     *
 *                                                                            *
 *             This software is distributed under the terms of the            *
 *                 GNU General Public Licence (GPL) version 3,                *
 *                    copied verbatim in the file "LICENSE".                  *
 *                                                                            *
 * In applying this license GSI does not waive the privileges and immunities  *
 * granted to it by virtue of its status as an Intergovernmental Organization *
 * or submit itself to any jurisdiction.                                      *
 ******************************************************************************/

#include "R3BEventHeader.h"
#include "R3BEventHeaderPropagator.h"
#include "R3BFileSource2.h"
#include "R3BNeulandHitPar.h"
#include "R3BNeulandMapToCalParTask.h"
#include "R3BNeulandMapToCalTask.h"
#include "R3BTCalContFact.h"
#include <FairParRootFileIo.h>
#include <FairRootFileSink.h>
#include <FairRunAna.h>
#include <FairRuntimeDb.h>
#include <R3BLogger.h>
#include <R3BNeulandCalToHitParTask.h>
#include <R3BProgramOptions.h>
#include <R3BTCalPar.h>
#include <TObjString.h>
#include <TROOT.h>
#include <TStopwatch.h>
#include <filesystem>
#include <fmt/core.h>
#include <fmt/std.h>
#include <iostream>

namespace fs = std::filesystem;
constexpr int DEFAULT_EVENT_NUM = 10;
constexpr int DEFAULT_RUN_ID = 999;
constexpr int DEFAULT_MIN_STAT = 20;
constexpr auto ERROR_SCALE = 10.F;

using namespace std::string_literals;
auto main(int argc, const char** argv) -> int
{
    auto timer = TStopwatch{};
    timer.Start();

    auto programOptions = R3B::ProgramOptions("options for neuland exp data analysis");
    auto help = programOptions.create_option<bool>("help,h", "help message", false);
    auto enable_mille = programOptions.create_option<bool>("mille", "enable millepede", false);
    auto logLevel = programOptions.create_option<std::string>("logLevel,v", "set log level of fairlog", "info");
    auto input_file = programOptions.create_option<std::string>("in,i", "set the input files (regex)");
    auto input_par = programOptions.create_option<std::string>("in-par,p", "set the input parameter");
    auto output_file = programOptions.create_option<std::string>("out,o", "set the output file");
    auto event_num = programOptions.create_option<int>("eventNum,n", "set the event number", DEFAULT_EVENT_NUM);
    auto error_scale = programOptions.create_option<float>(
        "event_scale,e", "set the error scale for global parameter measurements", ERROR_SCALE);
    auto min_stat =
        programOptions.create_option<int>("min-stat,m", "set minimun statistics for calibration", DEFAULT_MIN_STAT);
    auto run_num = programOptions.create_option<int>("runNum,r", "set the number of runs", 1);

    if (!programOptions.verify(argc, argv))
    {
        return EXIT_FAILURE;
    }

    FairLogger::GetLogger()->SetLogScreenLevel(logLevel().c_str());
    const auto outputDir = R3B::GetParentDir(output_file());

    const auto outputfile_path = fs::path{ output_file() };
    const auto outputParFileName =
        outputDir / fmt::format("{}.par{}", outputfile_path.stem().string(), outputfile_path.extension().string());
    const auto input_filenames = R3B::GetFilesFromRegex(input_file());

    R3BLOG(debug, fmt::format("input data file: {}", fmt::join(input_filenames, ";")).c_str());
    R3BLOG(debug, fmt::format("input data par file: {}", input_par()).c_str());
    R3BLOG(debug, fmt::format("output data file: {}", output_file()).c_str());
    R3BLOG(debug, fmt::format("output data par file: {}", outputParFileName).c_str());

    auto source = std::make_unique<R3BFileSource2>();
    for (const auto& filename : input_filenames)
    {
        source->AddFile(filename);
    }

    auto sinkFile = R3B::make_rootfile(outputfile_path.c_str(), "RECREATE");
    auto sink = std::make_unique<FairRootFileSink>(sinkFile.release());

    auto run = std::make_unique<FairRunAna>();
    source->SetRunId(DEFAULT_RUN_ID);
    FairRootManager::SetTreeName("evt");

    try
    {
        auto EvntHeader = std::make_unique<R3BEventHeader>();
        run->SetEventHeader(std::move(EvntHeader));
        run->SetSource(source.release());
        run->SetSink(sink.release());

        // Add analysis task --------------------------------------------------------
        // auto runIdTask = std::make_unique<R3BEventHeaderPropagator>();
        // run->AddTask(runIdTask.release());

        auto map2Cal = std::make_unique<R3B::Neuland::Map2CalTask>();
        map2Cal->SetTrigger(R3B::Neuland::CalTrigger::all);
        run->AddTask(map2Cal.release());

        auto cal2hit_method =
            (enable_mille.value()) ? R3B::Neuland::Cal2HitParMethod::Millipede : R3B::Neuland::Cal2HitParMethod::LSQT;
        auto cal2hitParTask = std::make_unique<R3B::Neuland::Cal2HitParTask>(cal2hit_method);
        cal2hitParTask->SetTrigger(R3B::Neuland::CalTrigger::all);
        cal2hitParTask->SetMinStat(min_stat.value());
        cal2hitParTask->SetErrorScale(error_scale.value());
        cal2hitParTask->SetMaxModuleNum(1300);
        run->AddTask(cal2hitParTask.release());

        // set par input/output--------------------------------------------------------
        auto* rtdb = run->GetRuntimeDb();

        auto fileio = std::make_unique<FairParRootFileIo>();
        fileio->open(input_par().c_str());
        auto parOut = std::make_unique<FairParRootFileIo>(true);
        parOut->open(outputParFileName.c_str());
        rtdb->setFirstInput(fileio.release());
        rtdb->setOutput(parOut.release());
        rtdb->saveOutput();

        run->Init();

        // histogram calibration
        run->Run(0, event_num.value() <= 0 ? 0 : event_num.value());

        if (enable_mille.value())
        {
            rtdb->saveOutput();

            // millepede calibration
            run->Run(0, event_num.value() <= 0 ? 0 : event_num.value());
        }

        timer.Stop();
        std::cout << "Cal level data written to file " << outputfile_path << "\n";
        std::cout << "Cal level parameter written to file " << outputParFileName << "\n";
        std::cout << "Macro finished successfully." << std::endl;
        std::cout << "Real time: " << timer.RealTime() << "s, CPU time: " << timer.CpuTime() << "s" << std::endl;
    }
    catch (R3B::runtime_error& ex)
    {
        std::cout << "A runtime error has occured: \n";
        std::cerr << ex.what();
        std::cout << "\n\n";
    }
    catch (R3B::logic_error& ex)
    {
        std::cout << "A logic error has occured: \n";
        std::cerr << ex.what();
        std::cout << "\n\n";
    }
    catch (std::exception& ex)
    {
        std::cout << "A error has occured: \n";
        std::cerr << ex.what();
        std::cout << "\n\n";
    }
    catch (...)
    {
        std::cout << "An unrecognized error has occured: \n";
        std::cout << "\n\n";
    }

    auto* sinkfile = run->GetSink();
    if (sinkfile != nullptr)
    {
        sinkfile->Close();
    }

    return 0;
}
