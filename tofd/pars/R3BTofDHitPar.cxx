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

#include "R3BTofDHitPar.h"
#include "R3BLogger.h"

#include "FairLogger.h"
#include "FairParamList.h"

R3BTofDHitPar::R3BTofDHitPar(const char* name, const char* title, const char* context, Bool_t own)
    : FairParGenericSet(name, title, context, own)
    , fHitParams(new TObjArray(NMODULEMAX))
    , fMapInit(kFALSE)
{
}

R3BTofDHitPar::~R3BTofDHitPar()
{
    clear();
    if (fHitParams)
    {
        delete fHitParams;
        fHitParams = NULL;
    }
}

void R3BTofDHitPar::putParams(FairParamList* list)
{
    R3BLOG(info, "putParams() called");
    if (!list)
    {
        R3BLOG(fatal, "Could not find FairParamList");
        return;
    }
    list->addObject(GetName(), fHitParams);
}

Bool_t R3BTofDHitPar::getParams(FairParamList* list)
{
    if (!list)
    {
        R3BLOG(fatal, "Could not find FairParamList");
        return kFALSE;
    }
    if (!list->fillObject(GetName(), fHitParams))
    {
        R3BLOG(error, "Could not initialize " << GetName());
        return kFALSE;
    }
    return kTRUE;
}

void R3BTofDHitPar::clear()
{
    status = kFALSE;
    resetInputVersions();
}

// ----  Method print ----------------------------------------------------------
void R3BTofDHitPar::print() { printParams(); }

void R3BTofDHitPar::printParams()
{
    R3BLOG(info, GetName() << "Time Calib. Parameters");

    R3BLOG(info, "Number of HIT Parameters " << fHitParams->GetEntries());
    for (Int_t i = 0; i < fHitParams->GetEntries(); i++)
    {
        R3BTofDHitModulePar* t_par = dynamic_cast<R3BTofDHitModulePar*>(fHitParams->At(i));
        LOG(info) << "-------------------------------------------------------";
        if (t_par)
        {
            t_par->printParams();
        }
    }
}

R3BTofDHitModulePar* R3BTofDHitPar::GetModuleParAt(Int_t plane, Int_t paddle)
{
    if (!fMapInit)
    {
        fIndexMap.clear();
        R3BTofDHitModulePar* par;
        Int_t tplane;
        Int_t tpaddle;
        Int_t index;
        for (Int_t i = 0; i < fHitParams->GetEntries(); i++)
        {
            par = dynamic_cast<R3BTofDHitModulePar*>(fHitParams->At(i));
            if (NULL == par)
            {
                continue;
            }
            tplane = par->GetPlane();
            tpaddle = par->GetPaddle();
            if (tplane < 1 || tplane > N_TOFD_HIT_PLANE_MAX || tpaddle < 1 || tpaddle > N_TOFD_HIT_PADDLE_MAX)
            {
                R3BLOG(error, "Error in plane/paddle indexing. " << tplane << " / " << tpaddle);
                continue;
            }
            index = (tplane - 1) * N_TOFD_HIT_PADDLE_MAX + tpaddle - 1;
            if (fIndexMap.find(index) != fIndexMap.end())
            {
                R3BLOG(error, "Parameter found more than once. " << tplane << " / " << tpaddle);
                continue;
            }
            fIndexMap[index] = i;
        }
        fMapInit = kTRUE;
    }

    if (plane < 1 || plane > N_TOFD_HIT_PLANE_MAX || paddle < 1 || paddle > N_TOFD_HIT_PADDLE_MAX)
    {
        R3BLOG(error, "Error in plane/paddle indexing. " << plane << " / " << paddle);
        return NULL;
    }
    Int_t index = (plane - 1) * N_TOFD_HIT_PADDLE_MAX + paddle - 1;

    if (fIndexMap.find(index) == fIndexMap.end())
    {
        R3BLOG(warn, "Parameter not found for: " << plane << " / " << paddle);
        return NULL;
    }
    Int_t arind = fIndexMap[index];
    R3BTofDHitModulePar* par = dynamic_cast<R3BTofDHitModulePar*>(fHitParams->At(arind));
    return par;
}

void R3BTofDHitPar::AddModulePar(R3BTofDHitModulePar* tch)
{
    fMapInit = kFALSE;
    fHitParams->Add(tch);
}

void R3BTofDHitPar::PrintModuleParams(Int_t plane, Int_t paddle)
{
    R3BTofDHitModulePar* par = GetModuleParAt(plane, paddle);
    if (NULL != par)
    {
        par->printParams();
    }
}

void R3BTofDHitPar::DrawModuleParams(Int_t plane, Int_t paddle)
{
    R3BTofDHitModulePar* par = GetModuleParAt(plane, paddle);
    if (NULL != par)
    {
        par->DrawParams();
    }
}

ClassImp(R3BTofDHitPar);
