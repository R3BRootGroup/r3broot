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

void testNeulandReconstruction()
{
    TStopwatch timer;
    timer.Start();

    FairLogger::GetLogger()->SetLogScreenLevel("debug");

    FairRunAna run;
    run.SetSource(new FairFileSource("test.digi.root"));
    run.SetSink(new FairRootFileSink("test.reco.root"));

    auto io = new FairParRootFileIo();
    io->open("test.neul.root");
    run.GetRuntimeDb()->setFirstInput(io);

    // Perfect Reco
    run.AddTask(new R3BNeulandMultiplicityCheat("NeulandPrimaryHits", "NeulandMultiplicityCheat"));
    run.AddTask(new R3BNeulandNeutronsCheat("NeulandMultiplicityCheat", "NeulandPrimaryHits", "NeulandNeutronsCheat"));

    // Calorimetric Reco
    run.AddTask(new R3BNeulandMultiplicityCalorimetric("NeulandClusters", "NeulandMultiplicityCalorimetric"));
    run.AddTask(new R3BNeulandNeutronsRValue(
        600, "NeulandMultiplicityCalorimetric", "NeulandClusters", "NeulandNeutronsCalorimetric"));

    // Bayes Multiplicity
    run.AddTask(new R3BNeulandMultiplicityBayes("NeulandClusters", "NeulandMultiplicityBayes"));
    run.AddTask(
        new R3BNeulandNeutronsRValue(600, "NeulandMultiplicityBayes", "NeulandClusters", "NeulandNeutronsBayes"));

    run.Init();
    run.Run(0, 0);

    timer.Stop();
    cout << "Macro finished successfully." << endl;
    cout << "Real time: " << timer.RealTime() << "s, CPU time: " << timer.CpuTime() << "s" << endl;
}
