/******************************************************************************
 *   Copyright (C) 2019 GSI Helmholtzzentrum für Schwerionenforschung GmbH    *
 *   Copyright (C) 2019-2024 Members of R3B Collaboration                     *
 *                                                                            *
 *             This software is distributed under the terms of the            *
 *                 GNU General Public Licence (GPL) version 3,                *
 *                    copied verbatim in the file "LICENSE".                  *
 *                                                                            *
 * In applying this license GSI does not waive the privileges and immunities  *
 * granted to it by virtue of its status as an Intergovernmental Organization *
 * or submit itself to any jurisdiction.                                      *
 ******************************************************************************/

void testNeulandDigitizer()
{
    TStopwatch timer;
    timer.Start();

    FairLogger::GetLogger()->SetLogScreenLevel("debug");

    FairRunAna run;
    run.SetSource(new FairFileSource("test.simu.root"));
    run.SetSink(new FairRootFileSink("test.digi.root"));

    auto io = new FairParRootFileIo();
    io->open("test.para.root");
    run.GetRuntimeDb()->setFirstInput(io);

    run.AddTask(new R3BNeulandDigitizer(R3BNeulandDigitizer::Options::neulandTamex));
    // run.AddTask(new R3BNeulandClusterFinder());
    // run.AddTask(new R3BNeulandPrimaryInteractionFinder());
    // run.AddTask(new R3BNeulandPrimaryClusterFinder());

    run.Init();
    run.Run(0, 0);

    timer.Stop();
    cout << "Macro finished successfully." << endl;
    cout << "Real time: " << timer.RealTime() << "s, CPU time: " << timer.CpuTime() << "s" << endl;
}
