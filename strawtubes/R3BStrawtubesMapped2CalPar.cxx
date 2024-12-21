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

#include "R3BStrawtubesMapped2CalPar.h"
#include "FairLogger.h"
#include "FairRuntimeDb.h"
#include "R3BStrawtubesMappedData.h"
#include "R3BTCalEngine.h"
#include "R3BTCalPar.h"
#include "TClonesArray.h"
#include <FairRootManager.h>

R3BStrawtubesMapped2CalPar::R3BStrawtubesMapped2CalPar(const char* name, Int_t iVerbose)
    : FairTask(name, iVerbose)
    , fMinStats(1e5)
    , fCalPar()
    , fEngine()
    , fMapped()
{
}

R3BStrawtubesMapped2CalPar::~R3BStrawtubesMapped2CalPar()
{
    if (fEngine)
    {
        delete fEngine;
    }
    if (fCalPar)
    {
        delete fCalPar;
    }
}

InitStatus R3BStrawtubesMapped2CalPar::Init()
{
    FairRootManager* rm = FairRootManager::Instance();
    if (!rm)
    {
        return kFATAL;
    }

    fMapped = dynamic_cast<TClonesArray*>(rm->GetObject("StrawtubesMapped"));
    if (!fMapped)
    {
        return kFATAL;
    }

    fCalPar = dynamic_cast<R3BTCalPar*>(FairRuntimeDb::instance()->getContainer("StrawtubesTCalPar"));
    fCalPar->setChanged();

    fEngine = new R3BTCalEngine(fCalPar, fMinStats);

    return kSUCCESS;
}

void R3BStrawtubesMapped2CalPar::Exec(Option_t* option)
{
    Int_t mapped_num = fMapped->GetEntries();
    for (Int_t i = 0; i < mapped_num; i++)
    {
        auto mapped = dynamic_cast<R3BStrawtubesMappedData*>(fMapped->At(i));
        if (!mapped)
        {
            LOG(error) << "R3BStrawtubesMapped2CalPar::Exec(): What is this crap?";
            continue;
        }

        fEngine->Fill(mapped->GetPlane(), mapped->GetStraw(), mapped->GetSide(), mapped->GetTimeFine());
    }
}

void R3BStrawtubesMapped2CalPar::FinishEvent() {}

void R3BStrawtubesMapped2CalPar::FinishTask()
{
    fEngine->CalculateParamVFTX();
    fCalPar->printParams();
}

void R3BStrawtubesMapped2CalPar::SetMinStats(Int_t min_stats) { fMinStats = min_stats; }

ClassImp(R3BStrawtubesMapped2CalPar)
