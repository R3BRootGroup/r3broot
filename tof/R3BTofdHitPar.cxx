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

#include "R3BTofdHitPar.h"

#include "FairLogger.h"
#include "FairParamList.h" // for FairParamList

ClassImp(R3BTofdHitPar);

R3BTofdHitPar::R3BTofdHitPar(const char* name, const char* title, const char* context, Bool_t own)
    : FairParGenericSet(name, title, context, own)
    , fHitParams(new TObjArray(NMODULEMAX))
    , fMapInit(kFALSE)
{
}

R3BTofdHitPar::~R3BTofdHitPar()
{
    if (fHitParams)
    {
        delete fHitParams;
        fHitParams = NULL;
    }
}

void R3BTofdHitPar::putParams(FairParamList* list)
{
    LOG(info) << "R3BTofdHitPar::putParams() called";
    if (!list)
    {
        return;
    }
    list->addObject(GetName(), fHitParams);
}

Bool_t R3BTofdHitPar::getParams(FairParamList* list)
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

void R3BTofdHitPar::clear() {}

void R3BTofdHitPar::printParams()
{

    LOG(info) << " -----------  " << GetName() << " Time Calib. Parameters -------------  ";

    LOG(info) << " Number of HIT Parameters " << fHitParams->GetEntries();
    for (Int_t i = 0; i < fHitParams->GetEntries(); i++)
    {
        R3BTofdHitModulePar* t_par = dynamic_cast<R3BTofdHitModulePar*>(fHitParams->At(i));
        LOG(info) << "----------------------------------------------------------------------";
        if (t_par)
        {
            t_par->printParams();
        }
    }
}

R3BTofdHitModulePar* R3BTofdHitPar::GetModuleParAt(Int_t plane, Int_t paddle)
{
    if (!fMapInit)
    {
        fIndexMap.clear();
        R3BTofdHitModulePar* par;
        Int_t tplane;
        Int_t tpaddle;
        Int_t index;
        for (Int_t i = 0; i < fHitParams->GetEntries(); i++)
        {
            par = dynamic_cast<R3BTofdHitModulePar*>(fHitParams->At(i));
            if (NULL == par)
            {
                continue;
            }
            tplane = par->GetPlane();
            tpaddle = par->GetPaddle();
            if (tplane < 1 || tplane > N_TOFD_HIT_PLANE_MAX || tpaddle < 1 || tpaddle > N_TOFD_HIT_PADDLE_MAX)
            {
                LOG(error) << "R3BTofdHitPar::GetModuleParAt : error in plane/paddle indexing. " << tplane << " / "
                           << tpaddle;
                continue;
            }
            index = (tplane - 1) * N_TOFD_HIT_PADDLE_MAX + tpaddle - 1;
            if (fIndexMap.find(index) != fIndexMap.end())
            {
                LOG(error) << "R3BTofdHitPar::GetModuleParAt : parameter found more than once. " << tplane << " / "
                           << tpaddle;
                continue;
            }
            fIndexMap[index] = i;
        }
        fMapInit = kTRUE;
    }

    if (plane < 1 || plane > N_TOFD_HIT_PLANE_MAX || paddle < 1 || paddle > N_TOFD_HIT_PADDLE_MAX)
    {
        LOG(error) << "R3BTofdHitPar::GetModuleParAt : error in plane/paddle indexing. " << plane << " / " << paddle;
        return NULL;
    }
    Int_t index = (plane - 1) * N_TOFD_HIT_PADDLE_MAX + paddle - 1;

    if (fIndexMap.find(index) == fIndexMap.end())
    {
        LOG(warn) << "R3BTofdHitPar::GetModuleParAt : parameter not found for: " << plane << " / " << paddle;
        return NULL;
    }
    Int_t arind = fIndexMap[index];
    R3BTofdHitModulePar* par = dynamic_cast<R3BTofdHitModulePar*>(fHitParams->At(arind));
    return par;
}

void R3BTofdHitPar::AddModulePar(R3BTofdHitModulePar* tch)
{
    fMapInit = kFALSE;
    fHitParams->Add(tch);
}

void R3BTofdHitPar::PrintModuleParams(Int_t plane, Int_t paddle)
{
    R3BTofdHitModulePar* par = GetModuleParAt(plane, paddle);
    if (NULL != par)
    {
        par->printParams();
    }
}

void R3BTofdHitPar::DrawModuleParams(Int_t plane, Int_t paddle)
{
    R3BTofdHitModulePar* par = GetModuleParAt(plane, paddle);
    if (NULL != par)
    {
        par->DrawParams();
    }
}
