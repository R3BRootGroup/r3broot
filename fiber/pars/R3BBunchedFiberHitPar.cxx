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

#include "R3BBunchedFiberHitPar.h"

#include "FairLogger.h"
#include "FairParamList.h" // for FairParamList

ClassImp(R3BBunchedFiberHitPar);

R3BBunchedFiberHitPar::R3BBunchedFiberHitPar(const char* name, const char* title, const char* context, Bool_t own)
    : FairParGenericSet(name, title, context, own)
    , fHitParams(new TObjArray(NMODULEMAX))
    , fMapInit(kFALSE)
{
}

R3BBunchedFiberHitPar::~R3BBunchedFiberHitPar()
{
    if (fHitParams)
    {
        delete fHitParams;
        fHitParams = NULL;
    }
}

void R3BBunchedFiberHitPar::putParams(FairParamList* list)
{
    LOG(info) << "R3BBunchedFiberHitPar::putParams() called";
    if (!list)
    {
        return;
    }
    list->addObject(GetName(), fHitParams);
}

Bool_t R3BBunchedFiberHitPar::getParams(FairParamList* list)
{
    if (!list)
    {
        return kFALSE;
    }
    if (!list->fillObject(GetName(), fHitParams))
    {
        return kFALSE;
    }
    return kTRUE;
}

void R3BBunchedFiberHitPar::clear() {}

void R3BBunchedFiberHitPar::printParams()
{

    LOG(info) << " -----------  " << GetName() << " Fiber Hit Parameters -------------  ";

    LOG(info) << " Number of HIT Parameters " << fHitParams->GetEntries();
    for (Int_t i = 0; i < fHitParams->GetEntries(); i++)
    {
        R3BBunchedFiberHitModulePar* t_par = dynamic_cast<R3BBunchedFiberHitModulePar*>(fHitParams->At(i));
        LOG(info) << "----------------------------------------------------------------------";
        if (t_par)
        {
            t_par->printParams();
        }
    }
}

R3BBunchedFiberHitModulePar* R3BBunchedFiberHitPar::GetModuleParAt(Int_t fiber)
{
    if (!fMapInit)
    {
        fIndexMap.clear();
        R3BBunchedFiberHitModulePar* par;
        Int_t tFiber;
        Int_t index;
        for (Int_t i = 0; i < fHitParams->GetEntries(); i++)
        {
            par = dynamic_cast<R3BBunchedFiberHitModulePar*>(fHitParams->At(i));
            if (NULL == par)
            {
                continue;
            }
            tFiber = par->GetFiber();
            if (tFiber < 1 || tFiber > N_FIBER_MAX)
            {
                LOG(error) << "R3BBunchedFiberHitPar::GetModuleParAt : error in fiber indexing. " << tFiber;
                continue;
            }
            index = tFiber - 1;
            if (fIndexMap.find(index) != fIndexMap.end())
            {
                LOG(error) << "R3BBunchedFiberHitPar::GetModuleParAt : parameter found more than once. " << tFiber;
                continue;
            }
            fIndexMap[index] = i;
        }
        fMapInit = kTRUE;
    }

    if (fiber < 1 || fiber > N_FIBER_MAX)
    {
        LOG(error) << "R3BBunchedFiberHitPar::GetModuleParAt : error in fiber indexing. " << fiber;
        return NULL;
    }
    Int_t index = fiber - 1;

    if (fIndexMap.find(index) == fIndexMap.end())
    {
        LOG(warn) << "R3BBunchedFiberHitPar::GetModuleParAt : parameter not found for: " << fiber;
        return NULL;
    }
    Int_t arind = fIndexMap[index];
    R3BBunchedFiberHitModulePar* par = dynamic_cast<R3BBunchedFiberHitModulePar*>(fHitParams->At(arind));
    return par;
}

void R3BBunchedFiberHitPar::AddModulePar(R3BBunchedFiberHitModulePar* tch)
{
    fMapInit = kFALSE;
    fHitParams->Add(tch);
}

void R3BBunchedFiberHitPar::PrintModuleParams(Int_t fiber)
{
    R3BBunchedFiberHitModulePar* par = GetModuleParAt(fiber);
    if (NULL != par)
    {
        par->printParams();
    }
}

void R3BBunchedFiberHitPar::DrawModuleParams(Int_t fiber)
{
    R3BBunchedFiberHitModulePar* par = GetModuleParAt(fiber);
    if (NULL != par)
    {
        par->DrawParams();
    }
}
