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

#include "TClonesArray.h"
#include "TMath.h"

#include "FairLogger.h"
#include "FairRootManager.h"
#include "FairRunAna.h"
#include "FairRuntimeDb.h"
#include <fstream>
#include <iomanip>
#include <iostream>

#include "R3BCalifaCrystalCal2CrystalCalPID.h"
#include "R3BCalifaCrystalCalData.h"

using namespace std;

R3BCalifaCrystalCal2CrystalCalPID::R3BCalifaCrystalCal2CrystalCalPID()
    : FairTask("R3B CALIFA QPID Calibrator")
    , nEvents(0)
    , fCrystalCalDataCA(0)
{
    // counter1=0;
    // counter2=0;

    ifstream fcal("qpidscale.txt");
    int ch = 0;
    int ch_file;
    double s, d;

    cout << "QPID calibration:" << endl << endl;

    while (fcal.good() && ch < 128)
    {
        fcal >> ch_file >> s >> d;
        if (ch != ch_file)
        {
            cerr << "Invalid channel: " << ch_file << " != " << ch << endl;
            continue;
        }
        scale[ch] = s;
        delta[ch++] = d;

        cout << dec << ch << " " << scientific << s << " " << d << endl;
    }
}

R3BCalifaCrystalCal2CrystalCalPID::~R3BCalifaCrystalCal2CrystalCalPID()
{
    LOG(info) << "R3BCalifaCrystalCal2CrystalCalPID: Delete instance";
    delete fCrystalCalDataCA;
}

InitStatus R3BCalifaCrystalCal2CrystalCalPID::Init()
{
    Register();
    return kSUCCESS;
}

void R3BCalifaCrystalCal2CrystalCalPID::SetParContainers()
{
    //  // Get run and runtime database
    //  FairRunAna* run = FairRunAna::Instance();
    //  if (!run) LOG(fatal) << "R3BCalifaCrystalCal2CrystalCalPID::SetParContainers: No analysis run";
    //
    //  FairRuntimeDb* rtdb = run->GetRuntimeDb();
    //  if (!rtdb) LOG(fatal) << "R3BCalifaCrystalCal2CrystalCalPID::SetParContainers: No runtime database";
    //
    //  if ( fVerbose && fCaloCalPar ) {
    //    LOG(info) << "R3BCalifaCrystalCal2CrystalCalPID::SetParContainers() ";
    //    LOG(info) << "Container R3BCalifaCrystalCal2CrystalCalPIDPar loaded ";
    //  }
}

void R3BCalifaCrystalCal2CrystalCalPID::Register()
{
    LOG(debug) << "Registering";
    FairRootManager* fMan = FairRootManager::Instance();
    if (!fMan)
    {
        LOG(fatal) << "Init: No FairRootManager";
        return;
    }
    fCrystalCalDataCA = dynamic_cast<TClonesArray*>(fMan->GetObject("CalifaCrystalCalData"));
    if (NULL == fCrystalCalDataCA)
    {
        LOG(fatal) << "Branch CalifaCrystalCalData not found";
    }
}

InitStatus R3BCalifaCrystalCal2CrystalCalPID::ReInit()
{
    SetParContainers();
    return kSUCCESS;
}

void R3BCalifaCrystalCal2CrystalCalPID::Exec(Option_t* option)
{
    LOG(debug) << "Calibrating CALIFA PID Data";

    if (++nEvents % 10000 == 0)
        LOG(info) << nEvents;

    R3BCalifaCrystalCalData* hit;
    int ch;
    double nf, ns, s, d;

    for (Int_t i = 0; i < fCrystalCalDataCA->GetEntries(); i++)
    {
        hit = dynamic_cast<R3BCalifaCrystalCalData*>(fCrystalCalDataCA->At(i));
        ch = hit->GetCrystalId();
        if (ch < 0 || ch >= 128)
        {
            cerr << "Invalid channel: " << ch << endl;
            continue;
        }
        nf = hit->GetNf();
        ns = hit->GetNs();

        s = 1. + (scale[ch] - 1.) * (nf + ns) / .1;
        d = delta[ch] * (nf + ns) / .1;

        hit->SetNf(nf * (1. + s) / 2. + ns * (1. - s) / 2. - d);
        hit->SetNs(nf * (1. - s) / 2. + ns * (1. + s) / 2. + d);
    }

    return;
}

void R3BCalifaCrystalCal2CrystalCalPID::Finish() {}

void R3BCalifaCrystalCal2CrystalCalPID::Reset()
{
    LOG(debug) << "Clearing CalifaCrystalCal Data Structure";
    if (fCrystalCalDataCA)
        fCrystalCalDataCA->Clear();
}

ClassImp(R3BCalifaCrystalCal2CrystalCalPID)
