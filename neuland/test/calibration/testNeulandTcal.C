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

void testNeulandTcal()
{
    TStopwatch timer;
    timer.Start();

    Bool_t failed = false;

    Int_t i, start[2] = { 1000, 2000 }, width = 2000, end[2] = { start[0] + 1 + width, start[1] + width };

    // Tacquila ------------------------------------------------------------------
    R3BTCalPar* par = new R3BTCalPar();
    R3BTCalEngine* engine = new R3BTCalEngine(par, 0);

    for (Int_t i = 0; i < 100000; i++)
    {
        engine->Fill(1, 1, 1, i % width + start[0]);
        engine->Fill(1, 1, 3, i % width + start[0]);
        engine->Fill(1, 1, 2, i % width + start[1]);
        engine->Fill(1, 1, 4, i % width + start[1]);
    }

    engine->CalculateParamTacquila();

    for (Int_t pmt = 0; pmt < 2; pmt++)
    {
        R3BTCalModulePar* mpar = par->GetModuleParAt(1, 1, 1 + pmt);
        mpar->print();

        for (i = 0; i < start[pmt] - 1; i++)
            if (mpar->GetTimeTacquila(i) != -10000)
            {
                failed = true;
                cout << "Expected tac-time of channel " << i << " to be -10000 but is " << mpar->GetTimeTacquila(i)
                     << endl;
            }

        Double_t tol = 0.002;
        const Double_t step = 25. / width;
        for (i = 0; i <= width; i++)
        {
            if (mpar->GetTimeTacquila(start[pmt] + i - 1) > ((i * step) + tol) ||
                mpar->GetTimeTacquila(start[pmt] + i - 1) < ((i * step) - tol))
            {
                failed = true;
                cout << "Expected tac-time of channel " << i + start[pmt] - 1 << " to be " << i * step << " ± " << tol
                     << " but is " << mpar->GetTimeTacquila(i + start[pmt] - 1) << endl;
            }
        }

        for (i = end[pmt] + 1; i < 5000; i++)
        {
            if (mpar->GetTimeTacquila(i) != -10000)
            {
                failed = true;
                cout << "Expected tac-time of channel " << i << " to be -10000 but is " << mpar->GetTimeTacquila(i)
                     << endl;
            }
        }
    }
    // ---------------------------------------------------------------------------

    // VFTX ----------------------------------------------------------------------
    delete par;
    delete engine;

    par = new R3BTCalPar();
    engine = new R3BTCalEngine(par, 0);

    for (Int_t i = 0; i < 100000; i++)
    {
        engine->Fill(1, 1, 1, i % width + start[0]);
        engine->Fill(1, 1, 3, i % width + start[0]);
        engine->Fill(1, 1, 2, i % width + start[1]);
        engine->Fill(1, 1, 4, i % width + start[1]);
    }

    engine->CalculateParamVFTX();

    for (Int_t pmt = 0; pmt < 2; pmt++)
    {
        R3BTCalModulePar* mpar = par->GetModuleParAt(1, 1, 1 + pmt);
        mpar->print();

        for (i = 0; i < start[pmt] - 1; i++)
            if (mpar->GetTimeVFTX(i) != -10000)
            {
                failed = true;
                cout << "Expected VFTX-time of channel " << i << " to be -10000 but is " << mpar->GetTimeVFTX(i)
                     << endl;
            }

        Double_t tol = 0.002;
        const Double_t step = 5. / width;
        for (i = 0; i <= width; i++)
        {
            if (mpar->GetTimeVFTX(start[pmt] + i - 1) > ((i * step) + tol) ||
                mpar->GetTimeVFTX(start[pmt] + i - 1) < ((i * step) - tol))
            {
                failed = true;
                cout << "Expected VFTX-time of channel " << i + start[pmt] - 1 << " to be " << i * step << " ± " << tol
                     << " but is " << mpar->GetTimeVFTX(i + start[pmt] - 1) << endl;
            }
        }

        for (i = end[pmt] + 1; i < 5000; i++)
        {
            if (mpar->GetTimeVFTX(i) != -10000)
            {
                failed = true;
                cout << "Expected VFTX-time of channel " << i << " to be -10000 but is " << mpar->GetTimeVFTX(i)
                     << endl;
            }
        }
    }
    // ---------------------------------------------------------------------------

    timer.Stop();
    Double_t rtime = timer.RealTime();
    Double_t ctime = timer.CpuTime();
    cout << endl << endl;
    cout << "Real time " << rtime << " s, CPU time " << ctime << "s" << endl << endl;

    if (failed)
        cout << endl << " Test failed! " << endl;
    else
        cout << endl << " Test successful! " << endl;
}
