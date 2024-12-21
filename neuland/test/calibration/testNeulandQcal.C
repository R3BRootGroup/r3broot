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

void testNeulandQcal()
{
    TStopwatch timer;
    timer.Start();

    TString outputFileName = "data_qcal.root"; // name of output file

    // Create analysis run -------------------------------------------------------
    FairRunAna* run = new FairRunAna();
    run->SetOutputFile(outputFileName.Data());
    // ---------------------------------------------------------------------------

    // ----- Runtime DataBase info -----------------------------------------------
    FairRuntimeDb* rtdb = run->GetRuntimeDb();
    // ---------------------------------------------------------------------------

    // ---------------------------------------------------------------------------
    run->AddTask(new R3BNeulandQCalFiller());
    // ---------------------------------------------------------------------------

    // ---------------------------------------------------------------------------
    R3BNeulandTacquilaMapped2QCalPar* pedpar = new R3BNeulandTacquilaMapped2QCalPar();
    pedpar->SetPlanes(1);
    run->AddTask(pedpar);
    // ---------------------------------------------------------------------------

    // Initialize ----------------------------------------------------------------
    FairLogger::GetLogger()->SetLogScreenLevel("info");
    run->Init();
    // ---------------------------------------------------------------------------

    // Run -----------------------------------------------------------------------
    run->Run(0, 10000);
    rtdb->saveOutput();
    // ---------------------------------------------------------------------------

    timer.Stop();
    Double_t rtime = timer.RealTime();
    Double_t ctime = timer.CpuTime();
    cout << endl << endl;
    cout << "Real time " << rtime << " s, CPU time " << ctime << "s" << endl << endl;

    R3BNeulandQCalPar* par = (R3BNeulandQCalPar*)rtdb->getContainer("NeulandQCalPar");
    Bool_t failed = false;

    for (Int_t i = 0; i < 100; i++)
        if (par->GetParAt(1, i / 2 + 1, i % 2 + 1) != i + 100)
        {
            failed = true;
            cout << "[" << 1 << ", " << i / 50 + 1 << ", " << i % 2 + 1 << "] Expected: " << i + 100 << " but value is "
                 << par->GetParAt(1, i / 2 + 1, i % 2 + 1) << endl;
        }

    if (failed)
        cout << " Test failed! " << endl;
    else
        cout << " Test successful! " << endl;
}
