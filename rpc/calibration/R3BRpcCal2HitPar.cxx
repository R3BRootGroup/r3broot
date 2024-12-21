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
#include "TMath.h"
#include "TObjArray.h"
#include "TVector3.h"

#include "FairLogger.h"
#include "FairRootManager.h"
#include "FairRunAna.h"
#include "FairRuntimeDb.h"

#include "R3BRpcCal2HitPar.h"
#include "R3BRpcCalData.h"

#include "R3BRpcHitPar.h"

#include <iostream>
#include <stdlib.h>

R3BRpcCal2HitPar::R3BRpcCal2HitPar()
    : R3BRpcCal2HitPar("R3B Rpc Tot Calibration Parameters Finder ", 1)
{
}

R3BRpcCal2HitPar::R3BRpcCal2HitPar(const char* name, Int_t iVerbose)
    : FairTask(name, iVerbose)
    , fHitPar(NULL)
    , fCalDataCA(NULL)
    , fNumChannels(64)
    , fDebugMode(false)
    , fRpcCalib(true)
{
    for (Int_t i = 0; i < N_STRIP_NB; i++)
    {
        fhPos[i] = NULL;
        fhTime[i] = NULL;
    }
}

R3BRpcCal2HitPar::~R3BRpcCal2HitPar()
{
    LOG(info) << "R3BRpcCal2HitPar: Delete instance";
    if (fCalDataCA)
        delete fCalDataCA;
}

void R3BRpcCal2HitPar::SetParContainers()
{
    // Parameter Container
    // Reading califaMappingPar from FairRuntimeDb
    FairRuntimeDb* rtdb = FairRuntimeDb::instance();
    if (!rtdb)
    {
        LOG(error) << "FairRuntimeDb not opened!";
    }
}

void R3BRpcCal2HitPar::SetParameter() {}

InitStatus R3BRpcCal2HitPar::Init()
{
    LOG(info) << "R3BRpcCal2HitPar::Init()";

    FairRootManager* rootManager = FairRootManager::Instance();
    if (!rootManager)
    {
        LOG(error) << "R3BRpcCal2HitPar::Init() FairRootManager not found";
        return kFATAL;
    }

    fCalDataCA = dynamic_cast<TClonesArray*>(rootManager->GetObject("R3BRpcCalData"));
    if (!fCalDataCA)
    {
        LOG(error) << "R3BRpcCal2HitPar::Init() R3BRpcCalData not found";
        return kFATAL;
    }

    FairRuntimeDb* rtdb = FairRuntimeDb::instance();
    if (!rtdb)
    {
        LOG(error) << "R3BRpcCal2HitPar::Init() FairRuntimeDb not found";
        return kFATAL;
    }

    fHitPar = dynamic_cast<R3BRpcHitPar*>(rtdb->getContainer("RpcHitPar"));
    if (!fHitPar)

    {
        LOG(error) << "R3BRpcCal2HitPar::Init() Couldn't get handle on RPCHitPar container";
        return kFATAL;
    }

    // Set container with mapping parameters
    SetParameter();

    return kSUCCESS;
}

InitStatus R3BRpcCal2HitPar::ReInit()
{
    SetParContainers();
    SetParameter();
    return kSUCCESS;
}

void R3BRpcCal2HitPar::Exec(Option_t* opt)
{

    // loop over strip data
    Int_t nHits = fCalDataCA->GetEntries();
    UInt_t iDetector = 0;
    double charge_left = -1000;
    double charge_right = -1000;
    double time_left = 0;
    double time_right = 0;
    double ichn_right = 0;
    double ichn_left = 0;

    UInt_t inum = 0;

    for (Int_t i = 0; i < nHits; i++)
    {
        auto map1 = dynamic_cast<R3BRpcCalData*>(fCalDataCA->At(i));
        iDetector = map1->GetDetId();
        inum = iDetector * 41 + map1->GetChannelId() - 1;

        if (iDetector == 0)
        {
            if (map1->GetTotR_B() >= charge_right)
            {
                charge_right = map1->GetTotR_B();
                time_right = map1->GetTimeR_B();
                ichn_right = map1->GetChannelId();
            }
            if (map1->GetTotL_T() >= charge_left)
            {

                charge_left = map1->GetTotL_T();
                time_left = map1->GetTimeL_T();
                ichn_left = map1->GetChannelId();
            }

            if (NULL == fhPos[inum])
            {
                char strName[255];
                sprintf(strName, "%s_poscaldata_%d", fHitPar->GetName(), inum);

                if (!fRpcCalib)
                {
                    fhPos[inum] = new TH1F(strName, "", 750, 0., 1500.);
                }
                else
                {
                    fhPos[inum] = new TH1F(strName, "", 400, -2000., 2000.);
                }
            }
            if (NULL == fhTime[inum])
            {
                char strName[255];
                sprintf(strName, "%s_timecaldata_%d", fHitPar->GetName(), inum);
                fhTime[inum] = new TH1F(strName, "", 800, -200, -150);
            }

            if (ichn_left == ichn_right)
            {

                if (!fRpcCalib)
                {
                    fhPos[inum]->Fill((time_left - time_right) / 2. -
                                      10 * (fHitPar->GetCalParams1()->GetAt(inum) - 200));
                }
                else
                {
                    fhPos[inum]->Fill((time_left - time_right) / 2.);
                }
                fhTime[inum]->Fill((time_left + time_right) / 2.);
            }
        }
        if (iDetector == 1)
        {
            if (NULL == fhPos[inum])
            {
                char strName[255];
                sprintf(strName, "%s_poscaldata_%d", fHitPar->GetName(), inum);
                fhPos[inum] = new TH1F(strName, "", 5000, -2500., 2500.);
            }
            fhPos[inum]->Fill((map1->GetTimeR_B() - map1->GetTimeL_T()) * CSCINT / 2.);
        }
    }
    return;
}

void R3BRpcCal2HitPar::Reset() {}

void R3BRpcCal2HitPar::FinishEvent() {}

void R3BRpcCal2HitPar::FinishTask()
{

    if (fRpcCalib)
    {
        CalculateParsStrip();
    }
    else
    {

        CalculateParsPmt();
    }

    fHitPar->setChanged();
    fHitPar->printParams();
    fHitPar->Write();
}

void R3BRpcCal2HitPar::CalculateParsStrip()
{

    for (int t = 0; t < N_STRIP_NB; t++)
    {

        if (NULL == fhPos[t] && t >= 41)
        {
            continue;
        }

        fHitPar->SetCalParams4(0, t);

        int bin_max = fhPos[t]->GetMaximumBin();

        float threshold = (fhPos[t]->GetBinContent(bin_max - 1) + fhPos[t]->GetBinContent(bin_max) +
                           fhPos[t]->GetBinContent(bin_max + 1)) /
                          120.;

        for (int i = 1; i <= bin_max; i++)
        {
            float mean =
                (fhPos[t]->GetBinContent(i - 1) + fhPos[t]->GetBinContent(i) + fhPos[t]->GetBinContent(i + 1)) / 3.;
            if (fhPos[t]->GetBinContent(i) >= threshold && mean >= threshold)
            {
                fHitPar->SetCalParams1(i, t);
                break;
            }
        }
        for (int i = 800; i > bin_max; i--)
        {
            if (fhPos[t]->GetBinContent(i) >= threshold)
            {
                fHitPar->SetCalParams2(i, t);
                break;
            }
        }
    }

    for (int t = 0; t < N_STRIP_NB; t++)
    {

        if (NULL == fhTime[t])
        {
            continue;
        }

        int bin_max = fhTime[t]->GetMaximumBin();

        fhTime[t]->GetXaxis()->SetRange(bin_max - 50, bin_max + 50);

        fhTime[t]->Fit("gaus");

        TF1* g = dynamic_cast<TF1*>(fhTime[t]->GetListOfFunctions()->FindObject("gaus"));
        fHitPar->SetCalParams3(g->GetParameter(1), t);
    }
}

void R3BRpcCal2HitPar::CalculateParsPmt()
{

    double counter = 0;
    double average = 0;
    std::vector<double> v;

    for (int t = 0; t < N_STRIP_NB; t++)
    {

        if (NULL == fhPos[t])
        {
            continue;
        }

        int bin_max = fhPos[t]->GetMaximumBin();
        if (t >= 41)
        {
            float threshold = (fhPos[t]->GetBinContent(bin_max - 1) + fhPos[t]->GetBinContent(bin_max) +
                               fhPos[t]->GetBinContent(bin_max + 1)) /
                              10.;

            for (int i = 1; i <= bin_max; i++)
            {
                if (fhPos[t]->GetBinContent(i) >= threshold)
                {
                    fHitPar->SetCalParams1(i, t);
                    break;
                }
            }
        }
        else
        {

            if (NULL == fhPos[t])
            {
                continue;
            }

            fhPos[t]->GetXaxis()->SetRange(bin_max - 50, bin_max + 50);

            fhPos[t]->Fit("gaus");

            TF1* g = dynamic_cast<TF1*>(fhPos[t]->GetListOfFunctions()->FindObject("gaus"));

            v.push_back(g->GetParameter(1));
            average += g->GetParameter(1);
            counter++;
        }
    }

    for (int i = 0; i < v.size(); i++)
    {
        fHitPar->SetCalParams4(v[i] - average / counter, i);
    }
}

ClassImp(R3BRpcCal2HitPar)
