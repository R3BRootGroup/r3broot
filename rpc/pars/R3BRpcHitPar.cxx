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

#include "R3BRpcHitPar.h"

#include "FairLogger.h"
#include "FairParamList.h"

#include "TMath.h"
#include "TString.h"

#include <iostream>

// ---- Standard Constructor ---------------------------------------------------
R3BRpcHitPar::R3BRpcHitPar(const char* name, const char* title, const char* context)
    : FairParGenericSet(name, title, context)
    , fNumChannels(N_STRIP_NB)
    , fHitCalParams1(new TArrayF(N_STRIP_NB))
    , fHitCalParams2(new TArrayF(N_STRIP_NB))
    , fHitCalParams3(new TArrayF(N_STRIP_NB))
    , fHitCalParams4(new TArrayF(N_STRIP_NB))
{
}

// ----  Destructor ------------------------------------------------------------
R3BRpcHitPar::~R3BRpcHitPar()
{
    clear();
    delete fHitCalParams1;
    delete fHitCalParams2;
    delete fHitCalParams3;
    delete fHitCalParams4;
}

// ----  Method clear ----------------------------------------------------------
void R3BRpcHitPar::clear()
{
    status = kFALSE;
    resetInputVersions();
}

// ----  Method putParams ------------------------------------------------------
void R3BRpcHitPar::putParams(FairParamList* list)
{
    LOG(info) << "RpcHitPar::putParams() called";
    if (!list)
    {
        return;
    }

    Int_t array_size = fNumChannels;
    LOG(info) << "Array Size: " << array_size;

    fHitCalParams1->Set(array_size);
    fHitCalParams2->Set(array_size);

    list->add("RpcHitPar1", *fHitCalParams1);
    list->add("RpcHitPar2", *fHitCalParams2);
    list->add("RpcHitPar3", *fHitCalParams3);
    list->add("RpcHitPar4", *fHitCalParams4);
    list->add("RPCChannelsNumberPar", fNumChannels);
}

// ----  Method getParams ------------------------------------------------------
Bool_t R3BRpcHitPar::getParams(FairParamList* list)
{
    LOG(info) << "RpcHitPar::getParams() called";
    if (!list)
    {
        return kFALSE;
    }

    if (!list->fill("RPCChannelsNumberPar", &fNumChannels))
    {
        return kFALSE;
    }

    Int_t array_size = fNumChannels;
    LOG(info) << "Nb_channels: " << array_size;
    fHitCalParams1->Set(array_size);
    fHitCalParams2->Set(array_size);
    fHitCalParams3->Set(array_size);
    fHitCalParams4->Set(array_size);

    if (!(list->fill("RpcHitPar1", fHitCalParams1)))
    {
        LOG(info) << "---Could not initialize RpcHitPar1";
        return kFALSE;
    }

    if (!(list->fill("RpcHitPar2", fHitCalParams2)))
    {
        LOG(info) << "---Could not initialize RpcHitPar2";
        return kFALSE;
    }

    if (!(list->fill("RpcHitPar3", fHitCalParams3)))
    {
        LOG(info) << "---Could not initialize RpcHitPar3";
        return kFALSE;
    }

    if (!(list->fill("RpcHitPar4", fHitCalParams4)))
    {
        LOG(info) << "---Could not initialize RpcHitPar4";
        return kFALSE;
    }

    return kTRUE;
}

// ----  Method print ----------------------------------------------------------
void R3BRpcHitPar::print() { printParams(); }

// ----  Method printParams ----------------------------------------------------
void R3BRpcHitPar::printParams()
{
    LOG(info) << "RpcHitPar::RPC Calibration Parameters: ";
    Int_t array_size = fNumChannels;

    for (Int_t i = 0; i < fNumChannels; i++)
    {
        LOG(info) << "Channel number: " << i + 1;
        LOG(info) << "Param1= " << fHitCalParams1->GetAt(i);
        LOG(info) << "Param2= " << fHitCalParams2->GetAt(i);
        LOG(info) << "Param3= " << fHitCalParams3->GetAt(i);
        LOG(info) << "Param4= " << fHitCalParams4->GetAt(i);
    }
}

ClassImp(R3BRpcHitPar)
