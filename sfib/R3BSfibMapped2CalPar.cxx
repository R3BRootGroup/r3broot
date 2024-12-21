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

#include "R3BSfibMapped2CalPar.h"
#include "FairLogger.h"
#include "FairRuntimeDb.h"
#include "R3BSfibMappedData.h"
#include "R3BTCalEngine.h"
#include "R3BTCalPar.h"
#include "TClonesArray.h"
#include <FairRootManager.h>
#include <cassert>

R3BSfibMapped2CalPar::R3BSfibMapped2CalPar(Int_t a_verbose,
                                           enum R3BTCalEngine::CTDCVariant a_ctdc_variant,
                                           Int_t a_update_rate,
                                           Int_t a_min_stats)
    : FairTask("R3BSfibMapped2CalPar", a_verbose)
    , fCTDCVariant(a_ctdc_variant)
    , fUpdateRate(a_update_rate)
    , fMinStats(a_min_stats)
{
}

R3BSfibMapped2CalPar::~R3BSfibMapped2CalPar()
{
    delete fTCalPar;
    delete fEngine;
}

InitStatus R3BSfibMapped2CalPar::Init()
{
    auto rm = FairRootManager::Instance();
    if (!rm)
    {
        return kFATAL;
    }

    fMapped = dynamic_cast<TClonesArray*>(rm->GetObject("SfibMapped"));
    if (!fMapped)
    {
        return kFATAL;
    }

    // container needs to be created in tcal/R3BTCalContFact.cxx AND R3BTCal needs
    // to be set as dependency in CMakeLists.txt in the detector directory.
    auto name = "SfibTCalPar";
    fTCalPar = dynamic_cast<R3BTCalPar*>(FairRuntimeDb::instance()->getContainer(name));
    if (!fTCalPar)
    {
        LOG(error) << "Could not get SfibTCalPar.";
        return kFATAL;
    }
    fTCalPar->setChanged();
    fEngine = new R3BTCalEngine(fTCalPar, fMinStats);

    return kSUCCESS;
}

void R3BSfibMapped2CalPar::Exec(Option_t* option)
{
    auto mapped_num = fMapped->GetEntriesFast();
    for (auto i = 0; i < mapped_num; i++)
    {
        auto mapped = dynamic_cast<R3BSfibMappedData*>(fMapped->At(i));
        assert(mapped);
        fEngine->Fill(
            1 + mapped->IsTop(), mapped->GetChannel() * 2 - (mapped->IsLeading() ? 1 : 0), 1, mapped->GetFine());
    }
}

void R3BSfibMapped2CalPar::FinishEvent() {}

void R3BSfibMapped2CalPar::FinishTask()
{
    fEngine->CalculateParamClockTDC(fCTDCVariant);
    fTCalPar->printParams();
}

void R3BSfibMapped2CalPar::SetUpdateRate(Int_t a_rate) { fUpdateRate = a_rate; }

void R3BSfibMapped2CalPar::SetMinStats(Int_t a_min_stats) { fMinStats = a_min_stats; }

ClassImp(R3BSfibMapped2CalPar)
