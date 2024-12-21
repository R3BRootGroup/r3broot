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

#include "R3BNeulandQCalPar.h"
#include "FairLogger.h"
#include "FairParamList.h" // for FairParamList

R3BNeulandQCalPar::R3BNeulandQCalPar(const char* name, const char* title, const char* context, Bool_t own)
    : FairParGenericSet(name, title, context, own)
{
}

R3BNeulandQCalPar::~R3BNeulandQCalPar() {}

void R3BNeulandQCalPar::putParams(FairParamList* list)
{
    LOG(info) << "R3BNeulandQCalPar::putParams() called";
    if (!list)
    {
        return;
    }
    list->add("QCalParams", fParams);
}

Bool_t R3BNeulandQCalPar::getParams(FairParamList* list)
{
    LOG(info) << "R3BNeulandQCalPar::getParams() called";
    if (!list)
    {
        return kFALSE;
    }
    if (!list->fill("QCalParams", &fParams))
    {
        return kFALSE;
    }
    return kTRUE;
}

void R3BNeulandQCalPar::SetSize(Int_t size)
{
    fParams.Set(size);
    fParams.Reset();
}

void R3BNeulandQCalPar::printParams()
{

    LOG(info) << " -----------  " << GetName() << " NeuLAND Parameters -------------  ";
    LOG(info) << " Number of Parameters " << fParams.GetSize();
    for (Int_t i = 0; i < fParams.GetSize(); i++)
        LOG(info) << "PMT Nr. " << i << " : pedestal level = " << fParams.At(i);
    LOG(info) << "----------------------------------------------------------------------";
}

ClassImp(R3BNeulandQCalPar);
