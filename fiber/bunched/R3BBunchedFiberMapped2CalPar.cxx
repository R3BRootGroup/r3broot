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

#include "R3BBunchedFiberMapped2CalPar.h"
#include "FairLogger.h"
#include "FairRuntimeDb.h"
#include "R3BFiberMappedData.h"
#include "R3BLogger.h"
#include "R3BTCalEngine.h"
#include "R3BTCalPar.h"
#include "TClonesArray.h"
#include <FairRootManager.h>
#include <cassert>

R3BBunchedFiberMapped2CalPar::R3BBunchedFiberMapped2CalPar(const char* a_name,
                                                           Int_t a_verbose,
                                                           enum Electronics a_spmt_electronics,
                                                           enum R3BTCalEngine::CTDCVariant a_ctdc_variant,
                                                           Int_t a_update_rate,
                                                           Int_t a_min_stats)
    : FairTask(TString("R3B") + a_name + "Mapped2CalPar", a_verbose)
    , fName(a_name)
    , fSPMTElectronics(a_spmt_electronics)
    , fCTDCVariant(a_ctdc_variant)
    , fUpdateRate(a_update_rate)
    , fMinStats(a_min_stats)
    , fMapped(NULL)
{
}

R3BBunchedFiberMapped2CalPar::~R3BBunchedFiberMapped2CalPar()
{
    if (fMAPMTTCalPar)
        delete fMAPMTTCalPar;
    if (fMAPMTEngine)
        delete fMAPMTEngine;
    if (fMAPMTTrigTCalPar)
        delete fMAPMTTrigTCalPar;
    if (fMAPMTTrigEngine)
        delete fMAPMTTrigEngine;
    if (fSPMTTCalPar)
        delete fSPMTTCalPar;
    if (fSPMTEngine)
        delete fSPMTEngine;
}

InitStatus R3BBunchedFiberMapped2CalPar::Init()
{
    // container needs to be created in tcal/R3BTCalContFact.cxx AND R3BTCal needs
    // to be set as dependency in CMakeLists.txt in the detector directory.
#define GET_TCALPAR(NAME)                                                                            \
    do                                                                                               \
    {                                                                                                \
        auto name = fName + #NAME "TCalPar";                                                         \
        f##NAME##TCalPar = dynamic_cast<R3BTCalPar*>(FairRuntimeDb::instance()->getContainer(name)); \
        if (!f##NAME##TCalPar)                                                                       \
        {                                                                                            \
            R3BLOG(error, "Could not get " << name);                                                 \
            return kFATAL;                                                                           \
        }                                                                                            \
        f##NAME##TCalPar->setChanged();                                                              \
        f##NAME##Engine = new R3BTCalEngine(f##NAME##TCalPar, fMinStats);                            \
    } while (0)

    auto rm = FairRootManager::Instance();
    R3BLOG_IF(fatal, !rm, "FairRootManager not found");

    fMapped = dynamic_cast<TClonesArray*>(rm->GetObject(fName + "Mapped"));
    R3BLOG_IF(fatal, !fMapped, fName + "Mapped not found");

    GET_TCALPAR(MAPMT);
    GET_TCALPAR(MAPMTTrig);
    GET_TCALPAR(SPMT);

    return kSUCCESS;
}

void R3BBunchedFiberMapped2CalPar::Exec(Option_t* option)
{
    auto mapped_num = fMapped->GetEntriesFast();
    for (auto i = 0; i < mapped_num; i++)
    {
        auto mapped = dynamic_cast<R3BFiberMappedData*>(fMapped->At(i));
        assert(mapped);
        auto channel = mapped->GetChannel();
        if (mapped->IsMAPMT())
        {
            fMAPMTEngine->Fill(1, channel * 2 - (mapped->IsLeading() ? 1 : 0), 1, mapped->GetFine());
        }
        else if (mapped->IsSPMT())
        {
            fSPMTEngine->Fill(1, channel * 2 - (mapped->IsLeading() ? 1 : 0), 1, mapped->GetFine());
        }
        else if (mapped->IsMAPMTTrigger())
        {
            fMAPMTTrigEngine->Fill(1, channel, 1, mapped->GetFine());
        }
    }
}

void R3BBunchedFiberMapped2CalPar::FinishTask()
{
    R3BLOG(info, "");
    fMAPMTEngine->CalculateParamClockTDC(fCTDCVariant);
    fMAPMTTrigEngine->CalculateParamClockTDC(fCTDCVariant);
    switch (fSPMTElectronics)
    {
        case CTDC:
            fSPMTEngine->CalculateParamClockTDC(fCTDCVariant);
            break;
        case TAMEX:
            fSPMTEngine->CalculateParamVFTX();
            break;
        default:
            assert(0 && "This should not happen!");
    }
    fMAPMTTCalPar->printParams();
    fMAPMTTrigTCalPar->printParams();
    fSPMTTCalPar->printParams();
}

void R3BBunchedFiberMapped2CalPar::SetUpdateRate(Int_t a_rate) { fUpdateRate = a_rate; }

void R3BBunchedFiberMapped2CalPar::SetMinStats(Int_t a_min_stats) { fMinStats = a_min_stats; }

ClassImp(R3BBunchedFiberMapped2CalPar);
