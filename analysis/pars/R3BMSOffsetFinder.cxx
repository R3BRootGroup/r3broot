/******************************************************************************
 *   Copyright (C) 2019 GSI Helmholtzzentrum für Schwerionenforschung GmbH    *
 *   Copyright (C) 2019-2024 Members of R3B Collaboration                     *
 *                                                                            *
 *             This software is distributed under the terms of the            *
 *                 GNU General Public Licence (GPL) version 3,                *
 *                    copied verbatim in the file "LICENSE".                  *
 *                                                                            *
 * In applying this license GSI does not waive the privileges and immunities  *
 * granted to it by virtue of its status as an Intergovernmental Organization *
 * or submit itself to any jurisdiction.                                      *
 ******************************************************************************/

// ROOT headers
#include "R3BSamplerMappedData.h"
#include "R3BShared.h"
#include "TClonesArray.h"
#include "TF1.h"
#include "TH1F.h"
#include "TMath.h"
#include "TObjArray.h"
#include "TRandom.h"
#include "TVector3.h"
#include <cstdlib>
#include <iostream>
// Fair headers
#include "FairLogger.h"
#include "FairRootManager.h"
#include "FairRunAna.h"
#include "FairRuntimeDb.h"

// R3B headers
#include "R3BEventHeader.h"
#include "R3BMSOffsetFinder.h"
#include "R3BMSOffsetPar.h"

// R3BMSOffsetFinder: Default Constructor --------------------------
R3BMSOffsetFinder::R3BMSOffsetFinder()
    : R3BMSOffsetFinder("R3B MS Offset Finder", 1)
{
}

// R3BMSOffsetFinder: Standard Constructor --------------------------
R3BMSOffsetFinder::R3BMSOffsetFinder(const TString& name, Int_t iVerbose)
    : FairTask(name, iVerbose)
{
}

void R3BMSOffsetFinder::SetParContainers()
{
    // Parameter Container
    if (auto* rtdb = FairRuntimeDb::instance(); rtdb == nullptr)
    {
        LOG(error) << "FairRuntimeDb not opened!";
    }
}

// -----   Public method Init   --------------------------------------------
InitStatus R3BMSOffsetFinder::Init()
{

    LOG(info) << "R3BMSOffsetFinder: Init";

    FairRootManager* rootManager = FairRootManager::Instance();
    if (rootManager == nullptr)
    {
        return kFATAL;
    }

    fSamplerMapped = dynamic_cast<TClonesArray*>(rootManager->GetObject("SamplerMapped"));
    if (fSamplerMapped == nullptr)
    {
        return kFATAL;
    }

    fSamplerMSMapped = dynamic_cast<TClonesArray*>(rootManager->GetObject("SamplerMSMapped"));
    if (fSamplerMSMapped == nullptr)
    {
        return kFATAL;
    }

    FairRuntimeDb* rtdb = FairRuntimeDb::instance();
    if (rtdb == nullptr)
    {
        return kFATAL;
    }

    fMSOffsetPar = dynamic_cast<R3BMSOffsetPar*>(rtdb->getContainer("MSOffsetPar"));
    if (fMSOffsetPar == nullptr)
    {
        LOG(error) << "R3BMSOffsetFinder::Init() Couldn't get handle on MSOffsetPar container";
        return kFATAL;
    }

    fh_Offset_Finder = R3B::root_owned<TH1F>("MSOffsetFinderHisto", "MSOffsetFinderHisto", 2000.0, -1000.0, 1000.0);

    return kSUCCESS;
}

// -----   Public method ReInit   --------------------------------------------
InitStatus R3BMSOffsetFinder::ReInit()
{
    SetParContainers();
    return kSUCCESS;
}

// -----   Public method Exec   --------------------------------------------
void R3BMSOffsetFinder::Exec(Option_t* /*opt*/)
{

    const Int_t sampHits = fSamplerMapped->GetEntriesFast();
    const Int_t sampmsHits = fSamplerMSMapped->GetEntriesFast();
    if ((sampHits == 0) || (sampmsHits != 1))
    {
        return;
    }

    if (auto* SAMPMSMapped = dynamic_cast<R3BSamplerMappedData*>(fSamplerMSMapped->At(0)); SAMPMSMapped != nullptr)
    {
        for (auto* SAMPMapped : TRangeDynCast<R3BSamplerMappedData>(fSamplerMapped))
        {
            fh_Offset_Finder->Fill(SAMPMapped->GetTime() - SAMPMSMapped->GetTime());
        }
    }
    else
    {
        LOG(error) << "SAMPMSMapped is nullptr";
        return;
    }
}

// ---- Public method Reset   --------------------------------------------------
void R3BMSOffsetFinder::Reset() {}

void R3BMSOffsetFinder::FinishEvent() {}

// ---- Public method Finish   --------------------------------------------------
void R3BMSOffsetFinder::FinishTask() { SearchMSOffset(); }

//------------------

void R3BMSOffsetFinder::SearchMSOffset()
{

    LOG(info) << "R3BMSOffsetFinder: Search MS Offset";
    if (fh_Offset_Finder->GetEntries() >= fMinStatistics)
    {
        auto MSpos = fh_Offset_Finder->GetMaximumBin();
        auto offset =
            (fh_Offset_Finder->GetBinCenter(MSpos)) -
            0.5; // The bin is unitary, so the Bin Center will always be X.5 (the center of a unitary bin). To fix that
                 // the -0.5 is added to the equation. The time here is in clock cycles, so thi -0.5 corrects that.
        fMSOffsetPar->SetMSOffset(offset);
    }

    // Fill container:

    fMSOffsetPar->setChanged();
}

ClassImp(R3BMSOffsetFinder) // NOLINT
