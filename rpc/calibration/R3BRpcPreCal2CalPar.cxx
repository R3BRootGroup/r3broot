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
#include "TF1.h"
#include "TObjArray.h"
#include "TVector3.h"

#include "FairLogger.h"
#include "FairRootManager.h"
#include "FairRunAna.h"
#include "FairRuntimeDb.h"

#include "R3BEventHeader.h"
#include "R3BRpcPreCal2CalPar.h"
#include "R3BRpcTotCalPar.h"

#include "R3BRpcPreCalData.h"
#include "TCanvas.h"

#include <iostream>
#include <stdlib.h>

R3BRpcPreCal2CalPar::R3BRpcPreCal2CalPar()
    : R3BRpcPreCal2CalPar("R3B RPC Calibration Parameters Finder ", 1)
{
}

R3BRpcPreCal2CalPar::R3BRpcPreCal2CalPar(const char* name, Int_t iVerbose)
    : FairTask(name, iVerbose)
    , fDebugMode(0)
    , fNumChannels(64)
    , fR3BEventHeader(NULL)
    , fTotCalPar(NULL)
    , fPreCalDataCA(NULL)
{
    for (Int_t i = 0; i < N_NUM; i++)
    {
        fhTot[i] = NULL;
    }
}

R3BRpcPreCal2CalPar::~R3BRpcPreCal2CalPar() { LOG(info) << "R3BRpcPreCal2CalPar: Delete instance"; }

InitStatus R3BRpcPreCal2CalPar::Init()
{
    LOG(info) << "R3BRpcPreCal2CalPar::Init()";

    FairRootManager* rootManager = FairRootManager::Instance();
    if (!rootManager)
    {
        LOG(error) << "R3BRpcPreCal2CalPar::Init() FairRootManager not found";
        return kFATAL;
    }

    fPreCalDataCA = dynamic_cast<TClonesArray*>(rootManager->GetObject("R3BRpcPreCalData"));
    if (!fPreCalDataCA)
    {
        LOG(error) << "R3BRpcPreCal2CalPar::Init() fPreCalDataCA not found";
        return kFATAL;
    }
    FairRuntimeDb* rtdbPar = FairRuntimeDb::instance();
    if (!rtdbPar)
    {
        LOG(error) << "R3BRpcPreCal2CalPar::Init() FairRuntimeDb not found";
        return kFATAL;
    }
    fR3BEventHeader = (R3BEventHeader*)rootManager->GetObject("EventHeader.");
    if (!fR3BEventHeader)
    {
        LOG(warning) << "EventHeader. not found";
        fR3BEventHeader = (R3BEventHeader*)rootManager->GetObject("R3BEventHeader");
    }
    else
        LOG(error) << "EventHeader. found";

    fTotCalPar = dynamic_cast<R3BRpcTotCalPar*>(rtdbPar->getContainer("RpcTotCalPar"));
    if (!fTotCalPar)
    {
        LOG(error) << "R3BRpcPreCal2CalPar::Init() Couldn't get handle on RpcTotCalPar container";
        return kFATAL;
    }

    // Set container with mapping parameter
    return kSUCCESS;
}

InitStatus R3BRpcPreCal2CalPar::ReInit() { return kSUCCESS; }

void R3BRpcPreCal2CalPar::Exec(Option_t* opt)
{
    // loop over strip data
    Int_t nHits = fPreCalDataCA->GetEntries();
    UInt_t iDetector = 0;
    bool tpat2 = fR3BEventHeader->GetTpat() & 0xf000;
    if (tpat2 > 0)
    {
        for (Int_t i = 0; i < nHits; i++)
        {
            auto map1 = dynamic_cast<R3BRpcPreCalData*>(fPreCalDataCA->At(i));
            iDetector = map1->GetDetId();
            UInt_t inum = (iDetector * 41 + map1->GetChannelId()) * 2 + map1->GetSide() - 2;
            if (iDetector == 0)
            {
                if (NULL == fhTot[inum])
                {
                    char strName[255];
                    sprintf(strName, "%s_totcaldata_%d", fTotCalPar->GetName(), inum);
                    fhTot[inum] = new TH1F(strName, "", 4800, 0, 120);
                }
                fhTot[inum]->Fill(map1->GetTot());
            }
            if (iDetector == 1)
            {
                if (NULL == fhTot[inum])
                {
                    char strName[255];
                    sprintf(strName, "%s_totcaldata_%d", fTotCalPar->GetName(), inum);
                    fhTot[inum] = new TH1F(strName, "", 4800, 0, 120);
                }
            }
        }
    }
}

void R3BRpcPreCal2CalPar::Reset() {}

void R3BRpcPreCal2CalPar::FinishEvent() {}

Double_t fitf(Double_t* x, Double_t* par)
{
    Double_t fitval = par[0] + par[1] * x[0];
    return fitval;
}

void R3BRpcPreCal2CalPar::FinishTask()
{
    for (int t = 0; t < N_NUM; t++)
    {
        if (NULL == fhTot[t])
        {
            continue;
        }
        std::vector<double> v;
        fhTot[t]->Smooth(1, "R");
        for (int s = 2; s <= 4800; s++)
        {
            v.push_back((fhTot[t]->GetBinContent(s) - fhTot[t]->GetBinContent(s - 1)) / (120. / 4800.));
        }
        int maxElementIndex = std::max_element(v.begin(), v.end()) - v.begin();
        int it = maxElementIndex;
        fTotCalPar->SetCalParams(it * (120. / 4800.), t);
    }
    fTotCalPar->setChanged();
    fTotCalPar->printParams();
    fTotCalPar->Write();
}

ClassImp(R3BRpcPreCal2CalPar)
