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
// -----                  R3BNeulandTcalFill                  -----
// -----             Created 27-01-2015 by M.Heil             -----
// ----------------------------------------------------------------

#include "R3BNeulandTcalFill.h"
#include "FairLogger.h"
#include "FairRootManager.h"
#include "FairRtdbRun.h"
#include "FairRunIdGenerator.h"
#include "FairRuntimeDb.h"
#include "R3BEventHeader.h"
#include "R3BPaddleTamexMappedData.h"
#include "R3BTCalEngine.h"
#include "R3BTCalPar.h"
#include "TClonesArray.h"
#include "TF1.h"
#include "TH1F.h"
#include <stdlib.h>

R3BNeulandTcalFill::R3BNeulandTcalFill()
    : FairTask("R3BNeulandTCalFill", 1)
    , fUpdateRate(1000000)
    , fMinStats(100000)
    , fTrigger(-1)
    , fNofPlanes(0)
    , fNofBars(0)
    , fNof17(0)
    , fNEvents(0)
    , fCal_Par(NULL)
{
}

R3BNeulandTcalFill::R3BNeulandTcalFill(const char* name, Int_t iVerbose)
    : FairTask(name, iVerbose)
    , fUpdateRate(1000000)
    , fMinStats(100000)
    , fTrigger(-1)
    , fNofPlanes(0)
    , fNofBars(0)
    , fNof17(0)
    , fNEvents(0)
    , fCal_Par(NULL)
{
}

R3BNeulandTcalFill::~R3BNeulandTcalFill()
{
    if (fCal_Par)
    {
        delete fCal_Par;
    }
    if (fEngine)
    {
        delete fEngine;
    }
}

InitStatus R3BNeulandTcalFill::Init()
{
    FairRootManager* rm = FairRootManager::Instance();
    if (!rm)
    {
        return kFATAL;
    }
    header = dynamic_cast<R3BEventHeader*>(rm->GetObject("EventHeader."));
    /*   if (!header)
       {
           return kFATAL;
       }*/
    fHits = dynamic_cast<TClonesArray*>(rm->GetObject("NeulandMappedData"));
    if (!fHits)
    {
        return kFATAL;
    }

    fCal_Par = dynamic_cast<R3BTCalPar*>(FairRuntimeDb::instance()->getContainer("LandTCalPar"));
    fCal_Par->setChanged();

    fEngine = new R3BTCalEngine(fCal_Par, fMinStats);

    return kSUCCESS;
}

void R3BNeulandTcalFill::Exec(Option_t*)
{
    /*   if (fTrigger >= 0)
       {
           if (header->GetTrigger() != fTrigger)
           {
               return;
           }
       }*/

    Int_t nHits = fHits->GetEntries();
    //    LOG(info) << "number of hits:" << nHits << "   " ;

    /*
        if (nHits > (fNofPMTs / 2))
        {
            return;
        }
    */

    R3BPaddleTamexMappedData* hit;
    Int_t iPlane;
    Int_t iBar;

    // Loop over mapped hits
    for (Int_t i = 0; i < nHits; i++)
    {
        hit = dynamic_cast<R3BPaddleTamexMappedData*>(fHits->At(i));
        if (!hit)
        {
            continue;
        }

        iPlane = hit->GetPlaneId();
        iBar = hit->GetBarId();

        if (iPlane > fNofPlanes)
        {
            LOG(error) << "R3BNeulandTcalFill::Exec() : more planes then expected! Plane: " << iPlane;
            continue;
        }
        if (iBar > fNofBars)
        {
            LOG(error) << "R3BNeulandTcalFill::Exec() : more bars then expected! Plane: " << iBar;
            continue;
        }

        //        if (hit->Is17())
        //        {
        // 17-th channel
        // MH            channel = fNofPMTs + hit->GetGtb() * 20 + hit->GetTacAddr();
        //        }
        //        else
        //        {
        // PMT signal
        // iSide = hit->GetSide();

        //            LOG(info) << "Plane: " << iPlane << " Bar: " << iBar << " Side: " << iSide << " Cal channel: " <<
        //            channel << "   " ;

        //        }

        // Fill TAC histogram
        fEngine->Fill(iPlane, iBar, 0, hit->GetFineTime1LE());
        fEngine->Fill(iPlane, iBar, 0 + 2, hit->GetFineTime1TE());

        fEngine->Fill(iPlane, iBar, 1, hit->GetFineTime2LE());
        fEngine->Fill(iPlane, iBar, 1 + 2, hit->GetFineTime2TE());
    }

    // Increment events
    fNEvents += 1;
}

void R3BNeulandTcalFill::FinishEvent() {}

void R3BNeulandTcalFill::FinishTask() { fEngine->CalculateParamVFTX(); }

ClassImp(R3BNeulandTcalFill)
