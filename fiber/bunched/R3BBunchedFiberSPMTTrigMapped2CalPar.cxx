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

#include "FairLogger.h"
#include "FairRuntimeDb.h"
#include <FairRootManager.h>

#include "R3BBunchedFiberSPMTTrigMapped2CalPar.h"
#include "R3BFiberMappedData.h"
#include "R3BLogger.h"
#include "R3BTCalEngine.h"
#include "R3BTCalPar.h"

#include "TClonesArray.h"

R3BBunchedFiberSPMTTrigMapped2CalPar::R3BBunchedFiberSPMTTrigMapped2CalPar(Int_t a_verbose,
                                                                           Int_t a_update_rate,
                                                                           Int_t a_min_stats)
    : FairTask("R3BBunchedFiberSMPTTrigMapped2CalPar", a_verbose)
    , fUpdateRate(a_update_rate)
    , fMinStats(a_min_stats)
{
}

R3BBunchedFiberSPMTTrigMapped2CalPar::~R3BBunchedFiberSPMTTrigMapped2CalPar()
{
    if (fTCalPar)
        delete fTCalPar;
    if (fEngine)
        delete fEngine;
}

InitStatus R3BBunchedFiberSPMTTrigMapped2CalPar::Init()
{
    auto rm = FairRootManager::Instance();
    R3BLOG_IF(fatal, !rm, "FairRootManager not found");

    fMapped = dynamic_cast<TClonesArray*>(rm->GetObject("BunchedFiberSPMTTrigMapped"));
    R3BLOG_IF(fatal, !fMapped, "BunchedFiberSPMTTrigMapped not found");

    // container needs to be created in tcal/R3BTCalContFact.cxx AND R3BTCal needs
    // to be set as dependency in CMakeLists.txt in the detector directory.
    auto name = "BunchedFiberSPMTTrigTCalPar";
    fTCalPar = dynamic_cast<R3BTCalPar*>(FairRuntimeDb::instance()->getContainer(name));
    if (!fTCalPar)
    {
        R3BLOG(error, "Could not get " << name);
        return kFATAL;
    }
    fEngine = new R3BTCalEngine(fTCalPar, fMinStats);

    return kSUCCESS;
}

void R3BBunchedFiberSPMTTrigMapped2CalPar::Exec(Option_t* option)
{
    auto mapped_num = fMapped->GetEntriesFast();
    for (auto i = 0; i < mapped_num; i++)
    {
        auto mapped = dynamic_cast<R3BFiberMappedData*>(fMapped->At(i));
        fEngine->Fill(1, mapped->GetChannel(), 1, mapped->GetFine());
    }
}

void R3BBunchedFiberSPMTTrigMapped2CalPar::FinishTask()
{
    R3BLOG(info, "");
    fEngine->CalculateParamVFTX();
    fTCalPar->printParams();
}

void R3BBunchedFiberSPMTTrigMapped2CalPar::SetUpdateRate(Int_t a_rate) { fUpdateRate = a_rate; }

void R3BBunchedFiberSPMTTrigMapped2CalPar::SetMinStats(Int_t a_min_stats) { fMinStats = a_min_stats; }

ClassImp(R3BBunchedFiberSPMTTrigMapped2CalPar);
