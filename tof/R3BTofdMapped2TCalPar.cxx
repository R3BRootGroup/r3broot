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
// -----            R3BTofdMapped2TCalPar (7ps VFTX)            -----
// -----           Created Apr 2016 by R.Plag             -----
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
 * This file handles 7ps VFTX only, hence we have one edge per channel.
 *
 */

#include "R3BTofdMapped2TCalPar.h"
#include "R3BEventHeader.h"
#include "R3BPaddleTamexMappedData.h"
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

#include <iostream>
#include <stdlib.h>

using namespace std;

R3BTofdMapped2TCalPar::R3BTofdMapped2TCalPar()
    : FairTask("R3BTofdMapped2TCalPar", 1)
    , fUpdateRate(1000000)
    , fMinStats(100000)
    , fTrigger(-1)
    , fNofPlanes(0)
    , fPaddlesPerPlane(0)
    , fNofModules(0)
    , fNEvents(0)
    , fCal_Par(nullptr)
    , fMapped(nullptr)
    , header(nullptr)
    , fEngine(nullptr)
{
}

R3BTofdMapped2TCalPar::R3BTofdMapped2TCalPar(const char* name, Int_t iVerbose)
    : FairTask(name, iVerbose)
    , fUpdateRate(1000000)
    , fMinStats(100000)
    , fTrigger(-1)
    , fNofPlanes(0)
    , fPaddlesPerPlane(0)
    , fNofModules(0)
    , fNEvents(0)
    , fCal_Par(nullptr)
    , fMapped(nullptr)
    , header(nullptr)
    , fEngine(nullptr)
{
}

R3BTofdMapped2TCalPar::~R3BTofdMapped2TCalPar()
{
    delete fCal_Par;
    delete fEngine;
}

InitStatus R3BTofdMapped2TCalPar::Init()
{
    FairRootManager* rm = FairRootManager::Instance();
    if (!rm)
    {
        return kFATAL;
    }

    header = (R3BEventHeader*)rm->GetObject("R3BEventHeader");
    // may be = NULL!

    fMapped = (TClonesArray*)rm->GetObject("TofdMapped");
    if (!fMapped)
    {
        return kFATAL;
    }

    // container needs to be created in tcal/R3BTCalContFact.cxx AND R3BTCal needs
    // to be set as dependency in CMakelists.txt (in this case in the tof directory)
    fCal_Par = (R3BTCalPar*)FairRuntimeDb::instance()->getContainer("TofdTCalPar");
    if (!fCal_Par)
    {
        LOG(error) << "R3BTofdMapped2TCalPar::Init() Couldn't get handle on TofdTCalPar. ";
        return kFATAL;
    }

    fCal_Par->setChanged();

    if (!fNofModules)
    {
        LOG(error) << "R3BTofdMapped2TCalPar::Init() Number of modules not set. ";
        return kFATAL;
    }

    fEngine = new R3BTCalEngine(fCal_Par, fMinStats);

    return kSUCCESS;
}

void R3BTofdMapped2TCalPar::Exec(Option_t* option)
{
    // test for requested trigger (if possible)
    if ((fTrigger >= 0) && (header) && (header->GetTrigger() != fTrigger))
        return;

    Int_t nHits = fMapped->GetEntries();

    // Loop over mapped hits
    for (Int_t i = 0; i < nHits; i++)
    {
        R3BTofdMappedData* hit = (R3BTofdMappedData*)fMapped->At(i);

        Int_t iDetector = hit->GetDetector(); // 1..n
        Int_t iSide = hit->GetSide();         // 1/2
        Int_t iChannel = hit->GetChannel();   // 1..n
        Int_t iEdge = hit->GetEdge();         // 1/2

        if (iPlane > fNofPlanes)
        {
            LOG(error) << "R3BTofdMapped2TCalPar::Exec() : more planes than expected! Plane: " << iPlane
                       << " allowed are 1.." << fNofPlanes;
            continue;
        }
        if (iBar > fPaddlesPerPlane)
        {
            LOG(error) << "R3BTofdMapped2TCalPar::Exec() : more bars then expected! Det: " << iBar << " allowed are 1.."
                       << fPaddlesPerPlane;
            continue;
        }

        Int_t edge = iSide * 2 + iEdge - 2; // 1..4
        fEngine->Fill(iDetector, iChannel, edge, hit->GetCoarseTime());
    }

    ++fNEvents;
}

void R3BTofdMapped2TCalPar::FinishEvent() {}

void R3BTofdMapped2TCalPar::FinishTask()
{
    fEngine->CalculateParamVFTX();
    fCal_Par->printParams();
}

ClassImp(R3BTofdMapped2TCalPar)
