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

#include "R3BTprevTnext.h"
#include "FairLogger.h"
#include "FairRootManager.h"
#include "FairRunAna.h"
#include "FairRuntimeDb.h"
#include "R3BEventHeader.h"
#include "R3BMSOffsetPar.h"
#include "R3BSamplerMappedData.h"
#include "R3BShared.h"
#include "TClonesArray.h"
#include "TMath.h"
#include "TObjArray.h"
#include "TRandom.h"
#include <cstdlib>
#include <iostream>
R3BTprevTnext::R3BTprevTnext()
    : R3BTprevTnext("R3B Tprev/Tnext", 1)
{
}

R3BTprevTnext::R3BTprevTnext(const TString& name, Int_t iVerbose)
    : FairTask(name, iVerbose)
{
}

InitStatus R3BTprevTnext::Init()
{

    FairRootManager* rootManager = FairRootManager::Instance();
    if (rootManager == nullptr)
    {
        LOG(error) << "FairRootManager not opened!";
        return kFATAL;
    }

    FairRuntimeDb* rtdb = FairRuntimeDb::instance();
    if (rtdb == nullptr)
    {
        LOG(error) << "FairRuntimeDb not opened!";
        return kFATAL;
    }

    fR3BEventHeader = dynamic_cast<R3BEventHeader*>(rootManager->GetObject("EventHeader."));
    if (fR3BEventHeader == nullptr)
    {
        LOG(error) << "R3BTprevTnext::Init() EventHeader. not found";
        return kFATAL;
    }

    fSamplerMapped = dynamic_cast<TClonesArray*>(rootManager->GetObject("SamplerMapped"));
    if (fSamplerMapped == nullptr)
    {
        LOG(error) << "SamplerMapped not found or wrong type.";
        return kFATAL;
    }

    fSamplerMSMapped = dynamic_cast<TClonesArray*>(rootManager->GetObject("SamplerMSMapped"));
    if (fSamplerMSMapped == nullptr)
    {
        LOG(error) << "SamplerMSMapped not found or wrong type.";
        return kFATAL;
    }

    fMSOffsetPar = dynamic_cast<R3BMSOffsetPar*>(rtdb->getContainer("MSOffsetPar"));
    if (fMSOffsetPar == nullptr)
    {
        LOG(error) << "Could not find MSOffsetPar container!";
        return kFATAL;
    }

    return kSUCCESS;
}

void R3BTprevTnext::Exec(Option_t* /*opt*/)
{
    Double_t fTprev = INVALID_TPTN; // Initialization with invalid values.
    Double_t fTnext = INVALID_TPTN;
    Int_t dMScounter = 0;
    const Int_t sampHits = fSamplerMapped->GetEntriesFast();
    const Int_t sampmsHits = fSamplerMSMapped->GetEntriesFast();
    const Double_t MSOffset = fMSOffsetPar->GetMSOffset();
    Double_t SAMPTime = 0;
    Double_t SAMPMSTime = 0;
    if (sampmsHits == 1 && sampHits > 0)
    {
        R3BSamplerMappedData* SAMPMapped = nullptr;
        R3BSamplerMappedData* SAMPMSMapped = nullptr;
        SAMPMSMapped = dynamic_cast<R3BSamplerMappedData*>(fSamplerMSMapped->At(0));
        SAMPMSTime = SAMPMSMapped->GetTime();
        for (Int_t i = 0; i < sampHits; ++i)
        {
            SAMPMapped = dynamic_cast<R3BSamplerMappedData*>(fSamplerMapped->At(i));
            SAMPTime = SAMPMapped->GetTime();
            auto tpn = SAMPTime - SAMPMSTime - MSOffset;
            if (tpn < -fDelta_clk)
            {
                fTprev = -CLOCK_PERIOD * tpn;
            }
            if (abs(tpn) <= fDelta_clk)
            {
                ++dMScounter;
            }
            if (tpn > fDelta_clk)
            {
                fTnext = CLOCK_PERIOD * tpn;
                break;
            }
        }
    }
    else
    {
        fR3BEventHeader->SetTprev(INVALID_EVENT);
        fR3BEventHeader->SetTnext(INVALID_EVENT);
        return;
    }

    if (dMScounter ==
        0) // This if statement checks if there is no Hit correspondent to the MSHit even after the first if statement
    {
        fTprev = ERROR_NO_MS;
        fTnext = ERROR_NO_MS;
    }
    else if (dMScounter > 1)
    {
        fTprev = ERROR_MULTI_MS;
        fTnext = ERROR_MULTI_MS;
    }
    else if (dMScounter == 1)
    {
        if (fTprev == INVALID_TPTN)
        {
            fTprev = fNo_record_time;
        }
        if (fTnext == INVALID_TPTN)
        {
            fTnext = fNo_record_time;
        }
    }
    fR3BEventHeader->SetTprev(fTprev);
    fR3BEventHeader->SetTnext(fTnext);
}

ClassImp(R3BTprevTnext) // NOLINT
