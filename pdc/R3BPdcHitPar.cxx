/******************************************************************************
 *   Copyright (C) 2019 GSI Helmholtzzentrum für Schwerionenforschung GmbH    *
 *   Copyright (C) 2019 Members of R3B Collaboration                          *
 *                                                                            *
 *             This software is distributed under the terms of the            *
 *                 GNU General Public Licence (GPL) version 3,                *
 *                    copied verbatim in the file "LICENSE".                  *
 *                                                                            *
 * In applying this license GSI does not waive the privileges and immunities  *
 * granted to it by virtue of its status as an Intergovernmental Organization *
 * or submit itself to any jurisdiction.                                      *
 ******************************************************************************/

#include "R3BPdcHitPar.h"
#include "FairLogger.h"
#include "FairParamList.h" // for FairParamList
#include "R3BPdcHitModulePar.h"

ClassImp(R3BPdcHitPar);

R3BPdcHitPar::R3BPdcHitPar(const char* name, const char* title, const char* context, Bool_t own)
    : FairParGenericSet(name, title, context, own)
    , fHitParams(new TObjArray(NMODULEMAX))
    , fMapInit(kFALSE)
{
}
// ----  Destructor ------------------------------------------------------------

R3BPdcHitPar::~R3BPdcHitPar()
{
    if (fHitParams)
    {
        delete fHitParams;
        fHitParams = NULL;
    }
}

void R3BPdcHitPar::putParams(FairParamList* list)
{
    LOG(info) << "R3BPdcHitPar::putParams() called";
    if (!list)
    {
        return;
    }
    list->addObject(GetName(), fHitParams);
}

Bool_t R3BPdcHitPar::getParams(FairParamList* list)
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

void R3BPdcHitPar::clear() {}

void R3BPdcHitPar::printParams()
{

    LOG(info) << " -----------  " << GetName() << " Pdc Hit Parameters -------------  ";

    LOG(info) << " Number of HIT Parameters " << fHitParams->GetEntries();
    for (Int_t i = 0; i < fHitParams->GetEntries(); i++)
    {
        R3BPdcHitModulePar* t_par = (R3BPdcHitModulePar*)fHitParams->At(i);
        LOG(info) << "----------------------------------------------------------------------";
        if (t_par)
        {
            t_par->printParams();
        }
    }
}

R3BPdcHitModulePar* R3BPdcHitPar::GetModuleParAt(Int_t wire)
{
    if (!fMapInit)
    {
        fIndexMap.clear();
        R3BPdcHitModulePar* par;
        Int_t tWire;
        Int_t index;
        for (Int_t i = 0; i < fHitParams->GetEntries(); i++)
        {
            par = (R3BPdcHitModulePar*)fHitParams->At(i);
            if (NULL == par)
            {
                continue;
            }
            tWire = par->GetWire();
            if (tWire < 1 || tWire > N_WIRE_MAX * N_PLANE_MAX_PDC)
            {
                LOG(error) << "R3BPdcHitPar::GetModuleParAt : error in Wire indexing. " << tWire;
                continue;
            }
            index = tWire - 1;
            if (fIndexMap.find(index) != fIndexMap.end())
            {
                LOG(error) << "R3BPdcHitPar::GetModuleParAt : parameter found more than once. " << tWire;
                continue;
            }
            fIndexMap[index] = i;
        }
        fMapInit = kTRUE;
    }

    if (wire < 1 || wire > N_WIRE_MAX * N_PLANE_MAX_PDC)
    {
        LOG(error) << "R3BPdcHitPar::GetModuleParAt : error in wire indexing. " << wire;
        return NULL;
    }
    Int_t index = wire - 1;

    if (fIndexMap.find(index) == fIndexMap.end())
    {
        LOG(warning) << "R3BPdcHitPar::GetModuleParAt : parameter not found for: " << wire;
        return NULL;
    }
    Int_t arind = fIndexMap[index];
    R3BPdcHitModulePar* par = (R3BPdcHitModulePar*)fHitParams->At(arind);
    return par;
}

void R3BPdcHitPar::AddModulePar(R3BPdcHitModulePar* tch)
{
    fMapInit = kFALSE;
    fHitParams->Add(tch);
}

void R3BPdcHitPar::PrintModuleParams(Int_t wire)
{
    R3BPdcHitModulePar* par = GetModuleParAt(wire);
    if (NULL != par)
    {
        par->printParams();
    }
}

void R3BPdcHitPar::DrawModuleParams(Int_t wire)
{
    R3BPdcHitModulePar* par = GetModuleParAt(wire);
    if (NULL != par)
    {
        par->DrawParams();
    }
}
