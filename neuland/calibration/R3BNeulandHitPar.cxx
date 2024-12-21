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

#include "R3BNeulandHitPar.h"
#include "R3BNeulandCommon.h"

#include "FairLogger.h"
#include "FairParamList.h" // for FairParamList

R3BNeulandHitPar::R3BNeulandHitPar(const char* name, const char* title, const char* context, Bool_t own)
    : FairParGenericSet(name, title, context, own)
    , fParams(new TObjArray(Neuland::MaxNumberOfBars))
    , fGlobalTimeOffset(Neuland::NaN)
    , fDistanceToTarget(Neuland::NaN)
    , fEnergyCut(0.0)
{
    fDistancesToFirstPlane.resize(Neuland::MaxNumberOfPlanes, 0);
    for (int p = 1; p < Neuland::MaxNumberOfPlanes; ++p)
    {
        fDistancesToFirstPlane[p] = p * Neuland::BarSize_Z;
    }
}

R3BNeulandHitPar::~R3BNeulandHitPar()
{
    if (fParams)
    {
        delete fParams;
        fParams = NULL;
    }
}

void R3BNeulandHitPar::putParams(FairParamList* list)
{
    LOG(info) << "R3BNeulandHitPar::putParams() called";
    if (!list)
    {
        return;
    }
    list->addObject(GetName(), fParams);
}

Bool_t R3BNeulandHitPar::getParams(FairParamList* list)
{
    LOG(info) << "R3BNeulandHitPar::getParams() called";
    if (!list)
    {
        return kFALSE;
    }
    if (!list->fillObject(GetName(), fParams))
    {
        return kFALSE;
    }
    return kTRUE;
}

void R3BNeulandHitPar::clear() {}

void R3BNeulandHitPar::printParams()
{

    LOG(info) << " -----------  " << GetName() << " NeuLAND Parameters -------------  ";

    LOG(info) << " Number of Parameters " << fParams->GetEntries();
    for (Int_t i = 0; i < fParams->GetEntries(); i++)
    {
        R3BNeulandHitModulePar* t_par = dynamic_cast<R3BNeulandHitModulePar*>(fParams->At(i));
        LOG(info) << "----------------------------------------------------------------------";
        if (t_par)
        {
            t_par->printParams();
        }
    }
}

ClassImp(R3BNeulandHitPar);
