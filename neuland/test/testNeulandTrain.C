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

// Note: Just to check if tasks are working. Will not produce reasonably trained models with few events.
void testNeulandTrain()
{
    TStopwatch timer;
    timer.Start();

    FairLogger::GetLogger()->SetLogScreenLevel("debug");

    FairRunAna run;
    run.SetSource(new FairFileSource("test.digi.root"));
    run.SetSink(new FairRootFileSink("test.train.root"));

    auto rtdb = run.GetRuntimeDb();
    auto io = new FairParRootFileIo(true);
    io->open("test.neul.root");
    rtdb->setOutput(io);

    run.AddTask(new R3BNeulandMultiplicityBayesTrain());
    run.AddTask(new R3BNeulandMultiplicityCalorimetricTrain());

    run.Init();
    run.Run(0, 0);
    rtdb->writeContainers();
    rtdb->writeVersions();
    rtdb->saveOutput();
    rtdb->print();
    rtdb->closeOutput();

    timer.Stop();
    cout << "Macro finished successfully." << endl;
    cout << "Real time: " << timer.RealTime() << "s, CPU time: " << timer.CpuTime() << "s" << endl;
}
