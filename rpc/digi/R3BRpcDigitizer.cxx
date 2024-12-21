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

#include "R3BRpcDigitizer.h"
#include "FairRootManager.h"
#include "R3BRpc.h"
#include "R3BRpcCalData.h"
#include "R3BRpcPoint.h"
#include "TArrayD.h"
#include "TClonesArray.h"
#include "TMath.h"
#include "TRandom.h"
#include "TVector3.h"
#include <iostream>
#include <stdlib.h>

R3BRpcDigitizer::R3BRpcDigitizer()
    : FairTask("R3B RPC Digitizer")
    , fRpcPointDataCA(NULL)
    , fRpcCalDataCA(NULL)
{
}

R3BRpcDigitizer::~R3BRpcDigitizer()
{
    LOG(info) << "R3BRpcDigitizer: Delete instance";

    if (fRpcPointDataCA)
    {
        fRpcPointDataCA->Delete();
        delete fRpcPointDataCA;
    }
    if (fRpcCalDataCA)
    {
        fRpcCalDataCA->Delete();
        delete fRpcCalDataCA;
    }
}

void R3BRpcDigitizer::SetParContainers()
{
    FairRuntimeDb* rtdb = FairRuntimeDb::instance();
    LOG_IF(error, !rtdb) << "R3BRpcDigitizer::FairRuntimeDb not opened!";

    fSim_Par = (R3BRpcPars4Sim*)rtdb->getContainer("rpcPars4Sim");
    if (!fSim_Par)
    {
        LOG(error) << "R3BRpcDigitizer::Init() Couldn't get handle on "
                      "rpcPars4Sim container";
    }
    else
    {
        LOG(info) << "R3BRpcDigitizer:: rpcPars4Sim container opened";
    }
}

void R3BRpcDigitizer::SetParameter()
{
    //--- Parameter Container ---
    fNumberOfChannels = fSim_Par->GetNumChannels(); // Number of Channels (example!)

    fSim_Par->printParams();

    LOG(info) << "R3BRpcDigitizer:: max channel ID " << fNumberOfChannels;
}

InitStatus R3BRpcDigitizer::Init()
{
    LOG(info) << "R3BRpcDigitizer::Init ";

    FairRootManager* rootManager = FairRootManager::Instance();
    if (!rootManager)
        LOG(fatal) << "Init: No FairRootManager";

    fRpcPointDataCA = (TClonesArray*)rootManager->GetObject("RPCPoint");
    if (!fRpcPointDataCA)
    {
        LOG(fatal) << "Init: No RpcPoint CA";
        return kFATAL;
    }

    fRpcCalDataCA = new TClonesArray("R3BRpcCalData", 10);
    rootManager->Register("RpcCalData", "RPC Cal", fRpcCalDataCA, kTRUE);

    SetParameter();
    return kSUCCESS;
}

void R3BRpcDigitizer::Exec(Option_t* option)
{
    // Reset entries in output arrays, local arrays
    Reset();

    // Reading the Input -- Point data --
    Int_t nHits = fRpcPointDataCA->GetEntries();
    if (!nHits)
        return;

    R3BRpcPoint** pointData = NULL;
    pointData = new R3BRpcPoint*[nHits];
    for (Int_t i = 0; i < nHits; i++)
        pointData[i] = (R3BRpcPoint*)(fRpcPointDataCA->At(i));

    Int_t channelId;
    Double_t time;
    Double_t energy;

    for (Int_t i = 0; i < nHits; i++)
    {
        channelId = pointData[i]->GetChannelId();
        time = pointData[i]->GetTime();
        energy = pointData[i]->GetEnergyLoss();

        Int_t nCals = fRpcCalDataCA->GetEntriesFast();
        Bool_t existHit = 0;
        if (nCals == 0)
            AddCal(channelId, energy, time, 0);
        else
        {
            for (Int_t j = 0; j < nCals; j++)
            {
                if (((R3BRpcCalData*)(fRpcCalDataCA->At(j)))->GetChannelId() == channelId)
                {
                    ((R3BRpcCalData*)(fRpcCalDataCA->At(j)))->AddMoreEnergy(energy);
                    if (((R3BRpcCalData*)(fRpcCalDataCA->At(j)))->GetTime() > time)
                    {
                        ((R3BRpcCalData*)(fRpcCalDataCA->At(j)))->SetTime(time);
                    }
                    existHit = 1; // to avoid the creation of a new Hit
                    break;
                }
            }
            if (!existHit)
                AddCal(channelId, energy, time, 0);
        }
        existHit = 0;
    }

    if (pointData)
        delete[] pointData;

    Int_t nCals = fRpcCalDataCA->GetEntriesFast();
    if (nCals == 0)
        return;
}

void R3BRpcDigitizer::Reset()
{
    // Clear the CA structure
    LOG(debug) << "Clearing RpcCalData Structure";
    if (fRpcCalDataCA)
        fRpcCalDataCA->Clear();

    ResetParameters();
}

R3BRpcCalData* R3BRpcDigitizer::AddCal(Int_t ident, Double_t energy, ULong64_t time, Double_t tot_energy)
{
    TClonesArray& clref = *fRpcCalDataCA;
    Int_t size = clref.GetEntriesFast();
    if (fVerbose > 1)
        LOG(info) << "-I- R3BRpcDigitizer: Adding RpcCalData "
                  << " with unique identifier " << ident << " entering with " << energy * 1e06 << " keV  Time=" << time
                  << " tot_energy=" << tot_energy;

    return new (clref[size]) R3BRpcCalData(ident, energy, time, tot_energy);
}

ClassImp(R3BRpcDigitizer);
