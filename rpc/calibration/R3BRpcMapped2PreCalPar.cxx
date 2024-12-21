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
#include "TObjArray.h"
#include "TVector3.h"

#include "FairLogger.h"
#include "FairRootManager.h"
#include "FairRunAna.h"
#include "FairRuntimeDb.h"

#include "R3BRpcMapped2PreCalPar.h"

#include "R3BRpcMappedData.h"

#include "R3BTCalEngine.h"
#include "R3BTCalPar.h"

#include <iostream>
#include <stdlib.h>

R3BRpcMapped2PreCalPar::R3BRpcMapped2PreCalPar()
    : R3BRpcMapped2PreCalPar("R3B RPC Calibration Parameters Finder ", 1)
{
}

R3BRpcMapped2PreCalPar::R3BRpcMapped2PreCalPar(const char* name, Int_t iVerbose)
    : FairTask(name, iVerbose)
    , fTCalPar(NULL)
    , fEngine()
    , fMappedDataCA(NULL)
    , fNumChannels(64)
    , fDebugMode(0)
{
}

R3BRpcMapped2PreCalPar::~R3BRpcMapped2PreCalPar()
{
    LOG(info) << "R3BRpcMapped2PreCalPar: Delete instance";

    delete fEngine;
}

InitStatus R3BRpcMapped2PreCalPar::Init()
{
    LOG(info) << "R3BRpcMapped2PreCalPar::Init()";

    FairRootManager* rootManager = FairRootManager::Instance();
    if (!rootManager)
    {
        LOG(error) << "R3BRpcMapped2PreCalPar::Init() FairRootManager not found";
        return kFATAL;
    }

    fMappedDataCA = dynamic_cast<TClonesArray*>(rootManager->GetObject("R3BRpcMappedData"));
    if (!fMappedDataCA)
    {
        LOG(error) << "R3BRpcMapped2PreCalPar::Init() fMappedDataCA not found";
        return kFATAL;
    }

    FairRuntimeDb* rtdb = FairRuntimeDb::instance();
    if (!rtdb)
    {
        LOG(error) << "R3BRpcMapped2PreCalPar::Init() FairRuntimeDb not found";
        return kFATAL;
    }

    fTCalPar = dynamic_cast<R3BTCalPar*>(rtdb->getContainer("RpcTCalPar"));
    if (!fTCalPar)
    {
        LOG(error) << "R3BRpcMapped2PreCalPar::Init() Couldn't get handle on RpcTCalPar container";
        return kFATAL;
    }

    int fMinStats = 10;

    fEngine = new R3BTCalEngine(fTCalPar, fMinStats);

    return kSUCCESS;
}

InitStatus R3BRpcMapped2PreCalPar::ReInit() { return kSUCCESS; }

void R3BRpcMapped2PreCalPar::Exec(Option_t* opt)
{

    // loop over the 3 mapped structures

    // loop over strip data
    Int_t nHits = fMappedDataCA->GetEntries();
    for (Int_t i = 0; i < nHits; i++)
    {
        auto map1 = dynamic_cast<R3BRpcMappedData*>(fMappedDataCA->At(i));

        UInt_t iDetector = map1->GetDetId();
        UInt_t iStrip = map1->GetChannelId();                      // now 1..41
        UInt_t iEdge_Side = map1->GetEdge() * 2 + map1->GetSide(); // 0,3
        fEngine->Fill(iDetector + 1, iStrip, iEdge_Side + 1, map1->GetFineTime());
    }
}

void R3BRpcMapped2PreCalPar::Reset() {}

void R3BRpcMapped2PreCalPar::FinishEvent() {}

void R3BRpcMapped2PreCalPar::FinishTask()
{
    fEngine->CalculateParamVFTX();
    fTCalPar->printParams();
}

ClassImp(R3BRpcMapped2PreCalPar)
