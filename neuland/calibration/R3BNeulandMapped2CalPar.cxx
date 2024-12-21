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

#include "R3BNeulandMapped2CalPar.h"
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
#include <iostream>
#include <signal.h>
#include <stdlib.h>

using namespace std;

R3BNeulandMapped2CalPar::R3BNeulandMapped2CalPar()
    : FairTask("R3BNeulandMapped2TCalPar", 1)
    , fMinStats(100000)
    , fTrigger(-1)
    //    , fNofPMTs(0)
    , fNEvents(0)
    , fCal_Par(NULL)
{
}

R3BNeulandMapped2CalPar::R3BNeulandMapped2CalPar(const char* name, Int_t iVerbose)
    : FairTask(name, iVerbose)
    , fMinStats(100000)
    , fTrigger(-1)
    //, fNofPMTs(0)
    , fNEvents(0)
    , fCal_Par(NULL)
{
}

R3BNeulandMapped2CalPar::~R3BNeulandMapped2CalPar()
{
    if (fEngine)
    {
        delete fEngine;
    }
}

InitStatus R3BNeulandMapped2CalPar::Init()
{
    FairRootManager* rm = FairRootManager::Instance();
    if (!rm)
    {
        return kFATAL;
    }
    header = dynamic_cast<R3BEventHeader*>(rm->GetObject("EventHeader."));
    if (!header)
    {
        return kFATAL;
    }
    fHits = dynamic_cast<TClonesArray*>(rm->GetObject("NeulandMappedData"));
    if (!fHits)
    {
        return kFATAL;
    }
    fHitsTrigger = dynamic_cast<TClonesArray*>(rm->GetObject("NeulandTrigMappedData"));
    if (!fHitsTrigger)
    {
        LOG(info) << "Branch NeulandTrigMapped not found";
    }

    // container needs to be created in tcal/R3BTCalContFact.cxx AND R3BTCal needs
    // to be set as dependency in CMakelists.txt (in this case in the land directory)
    fCal_Par = dynamic_cast<R3BTCalPar*>(FairRuntimeDb::instance()->getContainer("LandTCalPar"));
    fCal_Par->setChanged();

    fEngine = new R3BTCalEngine(fCal_Par, fMinStats);

    for (Int_t pln = 0; pln < fNofPlanes; pln++)
    {
        for (Int_t bar = 0; bar < fNofBarsPerPlane; bar++)
        {
            for (Int_t pmt = 0; pmt < 4; pmt++)
            {
                counts[pln][bar][pmt] = 0;
            }
        }
    }

    checkcounts = 0;

    return kSUCCESS;
}

void R3BNeulandMapped2CalPar::Exec(Option_t* option)
{

    if (checkcounts == fNofPMTs)
    {
        // std::cout << "done " << std::endl;
        //    raise(SIGINT);
    }

    if (fTrigger >= 0)
    {
        if (header->GetTrigger() != fTrigger)
        {
            return;
        }
    }

    Int_t nHits = fHits->GetEntries();
    if (nHits > (fNofPMTs / 2))
    {
        // return;
    }

    // Loop over mapped hits
    for (Int_t i = 0; i < nHits; i++)
    {
        auto hit = dynamic_cast<R3BPaddleTamexMappedData*>(fHits->At(i));
        if (!hit)
        {
            continue;
        }

        // Check bar ID
        Int_t iPlane = hit->GetPlaneId();
        Int_t iBar = hit->GetBarId();
        Int_t iSide = -1 == hit->fCoarseTime1LE ? 2 : 1;

        // fill 1    1LE
        //      2    1TE
        //      3    2LE
        //      4    2TE

        Int_t iFine;
        if (1 == iSide)
        {
            iFine = hit->fFineTime1LE;
        }
        else
        {
            iFine = hit->fFineTime2LE;
        }

        fEngine->Fill(iPlane, iBar, (iSide - 1) * 2 + 1, iFine);
        counts[iPlane - 1][iBar - 1][(iSide - 1) * 2]++;
        if (counts[iPlane - 1][iBar - 1][(iSide - 1) * 2] == fMinStats)
        {
            checkcounts++;
            std::cout << iPlane << "a     " << iBar << "   " << iSide << std::endl;
            std::cout << checkcounts << std::endl;
        }

        if (1 == iSide)
        {
            iFine = hit->fFineTime1TE;
        }
        else
        {
            iFine = hit->fFineTime2TE;
        }
        fEngine->Fill(iPlane, iBar, (iSide - 1) * 2 + 2, iFine);
        counts[iPlane - 1][iBar - 1][(iSide - 1) * 2 + 1]++;
        if (counts[iPlane - 1][iBar - 1][(iSide - 1) * 2 + 1] == fMinStats)
        {
            checkcounts++;
            std::cout << iPlane << "b     " << iBar << "   " << iSide << std::endl;
            std::cout << checkcounts << std::endl;
        }
    }

    // Loop over mapped triggers
    if (fHitsTrigger)
    {
        nHits = fHitsTrigger->GetEntriesFast();
        for (Int_t i = 0; i < nHits; i++)
        {
            auto hit = dynamic_cast<R3BPaddleTamexMappedData*>(fHitsTrigger->At(i));
            if (!hit)
            {
                continue;
            }

            // Check bar ID
            auto iBar = hit->GetBarId();
            auto iFine = hit->fFineTime1LE;

            fEngine->Fill(100, iBar, 10, iFine);
        }
    }

    // Increment events
    fNEvents += 1;
}

void R3BNeulandMapped2CalPar::FinishEvent() {}

void R3BNeulandMapped2CalPar::FinishTask()
{
    fEngine->CalculateParamVFTX();
    fCal_Par->printParams();
}

ClassImp(R3BNeulandMapped2CalPar)
