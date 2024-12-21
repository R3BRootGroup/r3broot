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

#include "R3BNeulandCheckMapping.h"
#include "R3BPaddleTamexMappedData.h"
#include <FairRootManager.h>
#include <TClonesArray.h>

#include <iostream>
#include <sstream>
#include <stdlib.h>

#define TIME_HVON_MAX 5
#define PLANE_START 12
#define PLANE_END 13

R3BNeulandCheckMapping::R3BNeulandCheckMapping()
    : FairTask("NeulandCheckMapping", 1)
    , timestamp0(0)
    , timestamp1(0)
    , fMapped(NULL)
    , fTrigger(-1)
{
}

R3BNeulandCheckMapping::R3BNeulandCheckMapping(const char* name, Int_t iVerbose)
    : FairTask(name, iVerbose)
    , timestamp0(0)
    , timestamp1(0)
    , fMapped(NULL)
    , fTrigger(-1)
{
}

R3BNeulandCheckMapping::~R3BNeulandCheckMapping() {}

InitStatus R3BNeulandCheckMapping::Init()
{
    std::cout << "Init " << std::endl;

    FairRootManager* rm = FairRootManager::Instance();
    if (!rm)
    {
        return kFATAL;
    }

    header = (R3BEventHeader*)rm->GetObject("EventHeader.");

    fMapped = (TClonesArray*)rm->GetObject("NeulandMappedData");
    if (!fMapped)
    {
        return kFATAL;
    }

    h_countsok = new TH1F("countsok", "countsok", 1200, 0.5, 1200.5);

    std::cout << "Setting all vtargets to -1...\n";
    for (Int_t pln = 0; pln < fNofPlanes; pln++)
    {
        for (Int_t bar = 0; bar < fPaddlesPerPlane; bar++)
        {
            for (Int_t pmt = 0; pmt < 2; pmt++)
            {

                std::ostringstream oss;
                oss << "nl:hv:p" << pln + 1 << "b" << bar + 1 << "t" << pmt + 1;

                auto vmon = oss.str() + ":vmon";
                auto vtarget = oss.str() + ":vtarget.A";

                auto& entry = ca[pln][bar][pmt];
                entry.group = epics.CreateGroup();

                entry.vmon = entry.group->CreateChannel(vmon);
                entry.vtarget = entry.group->CreateChannel(vtarget);

                entry.vtarget->Set(-1);
                entry.group->Commit();

                cntOk[pln][bar][pmt] = 0;
            }
        }
    }
    std::cout << "Waiting for tubes to go down...\n";
    for (Int_t pln = 0; pln < fNofPlanes; pln++)
    {
        for (Int_t bar = 0; bar < fPaddlesPerPlane; bar++)
        {
            for (Int_t pmt = 0; pmt < 2; pmt++)
            {

                auto& entrycheck = ca[pln][bar][pmt];

                for (;;)
                {
                    std::cout << "#";
                    entrycheck.group->Fetch();
                    if (entrycheck.vmon->Get() < 20.)
                        break;
                    sleep(1);
                }
            }
        }
    }
    std::cout << "\nAll off! Enabling the first tube...\n";
    iPlnHvOn = PLANE_START;
    iBarHvOn = 0;
    iPmtHvOn = 0;
    auto& entryon = ca[iPlnHvOn][iBarHvOn][iPmtHvOn];
    entryon.vtarget->Set(1050);
    entryon.group->Commit();
    for (;;)
    {
        entryon.group->Fetch();
        if (fabs(entryon.vmon->Get() - entryon.vtarget->Get()) < 15.)
            break;
        sleep(1);
    }
    timestamp0.Set();
    std::cout << "Done! Ready to rumble.\n";

    finished = false;

    return kSUCCESS;
}

void R3BNeulandCheckMapping::Exec(Option_t* option)
{
    if (finished)
        FinishTask();

    timestamp1.Set();

    // test git ig,  delete after test!

    // std::cout << "Exec " << iPlnHvOn << " B " << iBarHvOn << " T " << iPmtHvOn << std::endl;

    std::cout << "Exec ts0=" << timestamp0 << " ts1=" << timestamp1 << " diff=" << difftime(timestamp1, timestamp0)
              << std::endl;

    if (difftime(timestamp1, timestamp0) > TIME_HVON_MAX)
    {

        std::cout << "5 seconds? " << std::endl;

        auto& entryoff = ca[iPlnHvOn][iBarHvOn][iPmtHvOn];
        entryoff.vtarget->Set(-1);
        entryoff.group->Commit();

        for (;;)
        {
            std::cout << ".";
            entryoff.group->Fetch();
            if (entryoff.vmon->Get() < 20.)
                break;
            sleep(1);
        }

        iBarHvOn++;

        if (iBarHvOn == fPaddlesPerPlane)
        {
            iBarHvOn = 0;
            iPmtHvOn++;
            if (iPmtHvOn == 2)
            {
                iPmtHvOn = 0;
                iPlnHvOn++;
                if (iPlnHvOn > PLANE_END)
                {
                    finished = true;
                    iPlnHvOn = fNofPlanes - 1;
                    iBarHvOn = fPaddlesPerPlane - 1;
                    iPmtHvOn = 1;
                    std::cout << "ALL DONE!" << std::endl;
                }
            }
        }

        if (finished)
            FinishTask();

        std::cout << "Enabling " << iPlnHvOn << ':' << iBarHvOn << ':' << iPmtHvOn << "..." << std::flush;
        auto& entryon = ca[iPlnHvOn][iBarHvOn][iPmtHvOn];
        entryon.vtarget->Set(1050);
        entryon.group->Commit();

        for (;;)
        {
            std::cout << "$";
            entryon.group->Fetch();
            if (fabs(entryon.vmon->Get() - entryon.vtarget->Get()) < 15.)
                break;
            sleep(1);
        }
        std::cout << "Done!\n";

        timestamp0.Set();
    }

    // std::cout << "not 5 seconds? " << std::endl;

    if ((fTrigger >= 0) && (header) && (header->GetTrigger() != fTrigger))
        return;

    Int_t nHits = fMapped->GetEntries();

    // nHits should be 1 * number of planes --- check here?

    // Loop over mapped hits
    for (Int_t i = 0; i < nHits; i++)
    {

        R3BPaddleTamexMappedData* hit = (R3BPaddleTamexMappedData*)fMapped->At(i);

        if (!hit)
            continue; // should not happen

        Int_t iPlane = hit->GetPlaneId();
        Int_t iBar = hit->GetBarId();
        Int_t iSide = -1 == hit->fCoarseTime1LE ? 2 : 1;

        // check if there are signals in proper channel(s)
        //      std::cout << " CHECK!!!!  xP " << iPlane << " B " << iBar << " T " << iSide << std::endl;
        //      std::cout << " HV !!!!  xP " << iPlnHvOn << " B " << iBarHvOn << " T " << iPmtHvOn << std::endl;

        if (iPlane - 1 == iPlnHvOn && iBar - 1 == iBarHvOn && iSide - 1 == iPmtHvOn)
        {
            //	  std::cout << " GOOD!!!!  xP " << iPlane << " B " << iBar << " T " << iSide << std::endl;
            cntOk[iPlnHvOn][iBarHvOn][iPmtHvOn]++;
        }
        else
        {
            // std::cout << " BAD!!!!  xP " << iPlane << " B " << iBar << " T " << iSide << std::endl;
            // std::cout << " HV !!!!  xP " << iPlnHvOn << " B " << iBarHvOn << " T " << iPmtHvOn << std::endl;
            // WriteBadMappedChannel();
        }
    }
}

void R3BNeulandCheckMapping::FinishEvent() {}

void R3BNeulandCheckMapping::FinishTask()
{
    // to do: cleanup resources
    std::cout << " finish task " << std::endl;
    for (Int_t pln = 0; pln < fNofPlanes; pln++)
    {
        for (Int_t bar = 0; bar < fPaddlesPerPlane; bar++)
        {
            for (Int_t pmt = 0; pmt < 2; pmt++)
            {

                h_countsok->SetBinContent(pln * 100 + bar * 2 + pmt + 1, cntOk[pln][bar][pmt]);
            }
        }
    }
    h_countsok->Write();
}
