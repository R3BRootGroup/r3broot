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
// -----            R3BPdcMapped2CalPar                     -----
// -----           Created Apr 2016 by R.Plag                 -----
// ----------------------------------------------------------------

#include "R3BPdcMapped2CalPar.h"
#include <FairRootManager.h>

#include "TClonesArray.h"

#include "FairLogger.h"
#include "FairRuntimeDb.h"

#include "R3BPdcMappedData.h"
#include "R3BTCalEngine.h"
#include "R3BTCalPar.h"

R3BPdcMapped2CalPar::R3BPdcMapped2CalPar()
    : FairTask("R3BPdcMapped2CalPar", 1)
    , fUpdateRate(1000000)
    , fMinStats(100000)
    , fCalPar(nullptr)
    , fMapped(nullptr)
    , fEngine(nullptr)
{
}

R3BPdcMapped2CalPar::R3BPdcMapped2CalPar(const char* name, Int_t iVerbose)
    : FairTask(name, iVerbose)
    , fUpdateRate(1000000)
    , fMinStats(100000)
    , fCalPar(nullptr)
    , fEngine(nullptr)
{
}

R3BPdcMapped2CalPar::~R3BPdcMapped2CalPar()
{
    delete fCalPar;
    delete fEngine;
}

InitStatus R3BPdcMapped2CalPar::Init()
{
    FairRootManager* rm = FairRootManager::Instance();
    if (!rm)
    {
        return kFATAL;
    }

    fMapped = dynamic_cast<TClonesArray*>(rm->GetObject("PdcMapped"));
    if (!fMapped)
    {
        return kFATAL;
    }

    // container needs to be created in tcal/R3BTCalContFact.cxx AND R3BTCal needs
    // to be set as dependency in CMakelists.txt (in this case in the tof directory)
    fCalPar = dynamic_cast<R3BTCalPar*>(FairRuntimeDb::instance()->getContainer("PdcTCalPar"));
    if (!fCalPar)
    {
        LOG(error) << "R3BPdcMapped2CalPar::Init() Couldn't get handle on PdcTCalPar. ";
        return kFATAL;
    }

    fCalPar->setChanged();

    fEngine = new R3BTCalEngine(fCalPar, fMinStats);

    return kSUCCESS;
}

void R3BPdcMapped2CalPar::Exec(Option_t* option)
{
    Int_t nHits = fMapped->GetEntries();

    // Loop over mapped hits
    for (Int_t i = 0; i < nHits; i++)
    {
        auto mapped = dynamic_cast<R3BPdcMappedData const*>(fMapped->At(i));

        fEngine->Fill(mapped->GetPlaneId(), mapped->GetWireId(), mapped->GetEdgeId(), mapped->GetTimeFine());
    }
}

void R3BPdcMapped2CalPar::FinishEvent() {}

void R3BPdcMapped2CalPar::FinishTask()
{
    fEngine->CalculateParamClockTDC(R3BTCalEngine::CTDC_16_BWD_150);
    fCalPar->printParams();
}

void R3BPdcMapped2CalPar::SetUpdateRate(Int_t rate) { fUpdateRate = rate; }

void R3BPdcMapped2CalPar::SetMinStats(Int_t minStats) { fMinStats = minStats; }

ClassImp(R3BPdcMapped2CalPar)
