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
// -----            R3BTofiMapped2CalPar                     -----
// -----           Created March 2021 by M.Heil                 -----
// ----------------------------------------------------------------

#include "R3BTofiMapped2CalPar.h"
#include <FairRootManager.h>

#include "TClonesArray.h"

#include "FairLogger.h"
#include "FairRuntimeDb.h"

#include "R3BTCalEngine.h"
#include "R3BTCalPar.h"
#include "R3BTofiMappedData.h"

R3BTofiMapped2CalPar::R3BTofiMapped2CalPar()
    : FairTask("R3BTofiMapped2CalPar", 1)
    , fUpdateRate(1000000)
    , fMinStats(100000)
    , fNofPlanes(0)
    , fPaddlesPerPlane(0)
    , fNofModules(0)
    , fCalPar(nullptr)
    , fMapped(nullptr)
    , fMappedTrigger(nullptr)
    , fEngine(nullptr)
{
}

R3BTofiMapped2CalPar::R3BTofiMapped2CalPar(const char* name, Int_t iVerbose)
    : FairTask(name, iVerbose)
    , fUpdateRate(1000000)
    , fMinStats(100000)
    , fNofPlanes(0)
    , fPaddlesPerPlane(0)
    , fNofModules(0)
    , fCalPar(nullptr)
    , fMappedTrigger(nullptr)
    , fEngine(nullptr)
{
}

R3BTofiMapped2CalPar::~R3BTofiMapped2CalPar()
{
    delete fCalPar;
    delete fEngine;
}

InitStatus R3BTofiMapped2CalPar::Init()
{
    FairRootManager* rm = FairRootManager::Instance();
    if (!rm)
    {
        return kFATAL;
    }

    fMapped = dynamic_cast<TClonesArray*>(rm->GetObject("TofiMapped"));
    if (!fMapped)
    {
        return kFATAL;
    }
    fMappedTrigger = dynamic_cast<TClonesArray*>(rm->GetObject("TofiTriggerMapped"));
    if (!fMapped)
    {
        return kFATAL;
    }

    // container needs to be created in tcal/R3BTCalContFact.cxx AND R3BTCal needs
    // to be set as dependency in CMakelists.txt (in this case in the tof directory)
    fCalPar = dynamic_cast<R3BTCalPar*>(FairRuntimeDb::instance()->getContainer("TofiTCalPar"));
    if (!fCalPar)
    {
        LOG(error) << "R3BTofiMapped2CalPar::Init() Couldn't get handle on TofiTCalPar. ";
        return kFATAL;
    }

    fCalPar->setChanged();

    if (!fNofModules)
    {
        LOG(error) << "R3BTofiMapped2CalPar::Init() Number of modules not set. ";
        return kFATAL;
    }

    fEngine = new R3BTCalEngine(fCalPar, fMinStats);

    return kSUCCESS;
}

void R3BTofiMapped2CalPar::Exec(Option_t* option)
{
    Int_t nHits = fMapped->GetEntries();

    // Loop over mapped hits
    for (Int_t i = 0; i < nHits; i++)
    {
        auto mapped = dynamic_cast<R3BTofiMappedData const*>(fMapped->At(i));

        if (mapped->GetDetectorId() > fNofPlanes)
        {
            LOG(error) << "R3BTofiMapped2CalPar::Exec() : more planes than expected! Plane: " << mapped->GetDetectorId()
                       << " allowed are 1.." << fNofPlanes;
            continue;
        }
        if (mapped->GetBarId() > fPaddlesPerPlane)
        {
            LOG(error) << "R3BTofiMapped2CalPar::Exec() : more bars then expected! Det: " << mapped->GetBarId()
                       << " allowed are 1.." << fPaddlesPerPlane;
            continue;
        }

        Int_t edge = mapped->GetSideId() * 2 + mapped->GetEdgeId() - 2; // 1..4
        fEngine->Fill(mapped->GetDetectorId(), mapped->GetBarId(), edge, mapped->GetTimeFine());
    }

    nHits = fMappedTrigger->GetEntries();

    // Loop over mapped triggers
    for (Int_t i = 0; i < nHits; i++)
    {
        auto mapped = dynamic_cast<R3BTofiMappedData const*>(fMappedTrigger->At(i));

        if (mapped->GetDetectorId() != fNofPlanes + 1)
        {
            LOG(error) << "R3BTofiMapped2CalPar::Exec() : trigger plane incorrect! Plane: " << mapped->GetDetectorId()
                       << " not " << fNofPlanes + 1;
            continue;
        }

        fEngine->Fill(mapped->GetDetectorId(), mapped->GetBarId(), 1, mapped->GetTimeFine());
    }
}

void R3BTofiMapped2CalPar::FinishEvent() {}

void R3BTofiMapped2CalPar::FinishTask()
{
    fEngine->CalculateParamVFTX();
    fCalPar->printParams();
}

void R3BTofiMapped2CalPar::SetUpdateRate(Int_t rate) { fUpdateRate = rate; }

void R3BTofiMapped2CalPar::SetMinStats(Int_t minStats) { fMinStats = minStats; }

void R3BTofiMapped2CalPar::SetNofModules(Int_t nDets, Int_t nCh)
{
    fNofPlanes = nDets;
    fPaddlesPerPlane = nCh;
    fNofModules = nDets * nCh * 4;
}

ClassImp(R3BTofiMapped2CalPar)
