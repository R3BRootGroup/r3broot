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

// ----------------------------------------------------------------
// -----            R3BLosMapped2CalPar (7ps VFTX)            -----
// -----           Created Feb 4th 2016 by R.Plag             -----
// ----------------------------------------------------------------

/* Some notes:
 *
 * There are different versions of VFTX:
 * 10px delivering 8 leading edges in Ch 1-8 and 8 trailing edges in
 *      Ch 9-16. This one was used for LOS1 but is not used in the analysis
 *  7px delivering 8 leading edges only. Used for LOS2 and this is the
 *      LOS we use for analysis.
 *
 * For s438b we had no synchronisation between 50 MHz tacquila clock
 * and the 200 MHz VFTX clock so we need to always subtract the time
 * of the master trigger from the LOS time.
 * The master trigger is on the last channel of the VFTX and handled
 * as 5th los channel.
 *
 * This file handles 7ps VFTX and TAMEX2, hence we have three times per channel.
 *
 */

#include "R3BLosMapped2CalPar.h"
#include "R3BEventHeader.h"
#include "R3BLogger.h"
#include "R3BLosMappedData.h"
#include "R3BTCalEngine.h"
#include "R3BTCalPar.h"

#include "FairLogger.h"
#include "FairRootManager.h"
#include "FairRtdbRun.h"
#include "FairRunIdGenerator.h"
#include "FairRuntimeDb.h"

#include "TClonesArray.h"
#include "TF1.h"
#include "TH1F.h"
#include "math.h"

#include <iostream>
#include <stdlib.h>

using namespace std;

R3BLosMapped2CalPar::R3BLosMapped2CalPar()
    : R3BLosMapped2CalPar("R3BLosMapped2CalPar", 1)
{
}

R3BLosMapped2CalPar::R3BLosMapped2CalPar(const char* name, Int_t iVerbose)
    : FairTask(name, iVerbose)
    , fUpdateRate(1000000)
    , fMinStats(100000)
    , fTrigger(-1)
    , fNofDetectors(0)
    , fNofChannels(0)
    , fNofTypes(0)
    , fNEvents(0)
    , fCal_Par(NULL)
{
}

R3BLosMapped2CalPar::~R3BLosMapped2CalPar()
{
    R3BLOG(debug1, "Destructor");
    if (fCal_Par)
    {
        delete fCal_Par;
    }
    if (fEngine)
    {
        delete fEngine;
    }
}

InitStatus R3BLosMapped2CalPar::Init()
{
    R3BLOG(info, "");
    for (UInt_t i = 0; i < 16; i++)
    {
        for (UInt_t k = 0; k < 3; k++)
        {
            Icount[i][k] = 0;
            Icounttrig[i][k] = 0;
        }
    }
    FairRootManager* rm = FairRootManager::Instance();
    if (!rm)
    {
        R3BLOG(fatal, "FairRootManager not found");
        return kFATAL;
    }

    header = dynamic_cast<R3BEventHeader*>(rm->GetObject("EventHeader."));
    // may be = NULL!
    if (!header)
    {
        R3BLOG(warn, "EventHeader. not found");
        header = dynamic_cast<R3BEventHeader*>(rm->GetObject("R3BEventHeader"));
    }
    else
        R3BLOG(info, "EventHeader. found");

    fMapped = dynamic_cast<TClonesArray*>(rm->GetObject("LosMapped"));
    if (!fMapped)
    {
        R3BLOG(fatal, "LosMapped not found");
        return kFATAL;
    }

    // get access to Trigger Mapped data
    fMappedTriggerItems = dynamic_cast<TClonesArray*>(rm->GetObject("LosTriggerMapped"));
    R3BLOG_IF(warn, !fMappedTriggerItems, "LosTriggerMapped not found");

    fCal_Par = dynamic_cast<R3BTCalPar*>(FairRuntimeDb::instance()->getContainer("LosTCalPar"));
    fEngine = new R3BTCalEngine(fCal_Par, fMinStats);

    return kSUCCESS;
}

void R3BLosMapped2CalPar::Exec(Option_t* option)
{
    // test for requested trigger (if possible)
    if ((fTrigger >= 0) && (header) && (header->GetTrigger() != fTrigger))
        return;

    UInt_t nHits = fMapped->GetEntries();
    // Loop over mapped hits
    for (UInt_t i = 0; i < nHits; i++)
    {

        R3BLosMappedData* hit = dynamic_cast<R3BLosMappedData*>(fMapped->At(i));
        if (!hit)
        {
            continue; // should not happen
        }

        // channel numbers are supposed to be 1-based (1..n)
        UInt_t iDetector = hit->GetDetector() - 1; // now 0..n-1
        UInt_t iChannel = hit->GetChannel() - 1;   // now 0..n-1
        UInt_t iType = hit->GetType();             // 0,1,2,3

        // cout<<"Mapped2CalPar "<<iDetector<<", "<<iChannel<<", "<<iType<<endl;
        if (iType < 3)
        {
            if (iDetector > (fNofDetectors - 1))
            {
                R3BLOG(error,
                       "More detectors than expected! Det: " << (iDetector + 1) << " allowed are 1.." << fNofDetectors);
                continue;
            }
            if (iChannel > (fNofChannels - 1))
            {
                R3BLOG(
                    error,
                    "More channels than expected! Channel: " << (iChannel + 1) << " allowed are 1.." << fNofChannels);
                continue;
            }

            if (iType > 3)
            {
                R3BLOG(error, "More time-types than expected! Type: " << iType << " allowed are 0..3");
                continue;
            }

            Icount[iChannel][iType]++;

            fEngine->Fill(iDetector + 1, iChannel + 1, iType + 1, hit->GetTimeFine());
        }
    }

    // Calibrate trigger channels -----------------------------------------------
    if (fMappedTriggerItems && fMappedTriggerItems->GetEntriesFast() > 0)
    {
        auto mapped_num = fMappedTriggerItems->GetEntriesFast();
        for (Int_t mapped_i = 0; mapped_i < mapped_num; mapped_i++)
        {
            auto mapped = dynamic_cast<R3BLosMappedData const*>(fMappedTriggerItems->At(mapped_i));

            UInt_t iDetector = mapped->GetDetector() - 1; // now 0..n-1
            UInt_t iChannel = mapped->GetChannel();
            UInt_t iType = mapped->GetType() + 1; // 1,2,3...
            R3BLOG(debug1, "Det: " << iDetector << " channel" << iChannel << " raw " << mapped->GetTimeFine());
            fEngine->Fill(3 + iDetector, iChannel, iType, mapped->GetTimeFine());
            Icounttrig[iChannel - 1][iType - 1]++;
        }
    }

    // Increment events
    fNEvents += 1;
}

void R3BLosMapped2CalPar::FinishTask()
{
    fEngine->CalculateParamVFTX();
    fCal_Par->printParams();
    fCal_Par->setChanged();

    R3BLOG(info, "Calibration of LOS detector");
    for (Int_t i = 0; i < 16; i++)
    {
        for (Int_t k = 0; k < 3; k++)
        {
            if (Icount[i][k] > fMinStats)
            {
                R3BLOG(info, "Channel: " << i + 1 << ", Type: " << k << ", Count: " << Icount[i][k]);
            }
        }
    }

    R3BLOG(info, "Calibration of trigger signals from LOS detector");
    for (Int_t i = 0; i < 16; i++)
    {
        for (Int_t k = 0; k < 3; k++)
        {
            if (Icounttrig[i][k] > fMinStats)
            {
                R3BLOG(info, "Channel: " << i + 1 << ", Type: " << k << ", Count: " << Icounttrig[i][k]);
            }
        }
    }
}

ClassImp(R3BLosMapped2CalPar);
