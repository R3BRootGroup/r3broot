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

// ------------------------------------------------------------------
// -----            R3BMusicHitPar source file                  -----
// -----       Created 24/11/19  by J.L. Rodriguez-Sanchez      -----
// ------------------------------------------------------------------

#include "R3BMusicHitPar.h"

#include "FairLogger.h"
#include "FairParamList.h"

#include "TArrayF.h"
#include "TMath.h"
#include "TString.h"

#include <iostream>

// ---- Standard Constructor ---------------------------------------------------
R3BMusicHitPar::R3BMusicHitPar(const char* name, const char* title, const char* context)
    : FairParGenericSet(name, title, context)
    , fNumParamsZFit(2)
    , fNumAnodes(8)
{
    fDetZHitParams = new TArrayF(fNumParamsZFit); // 2 Parameters for Z (Linear fits)
    fIn_use = new TArrayI(fNumAnodes);
    fAnode_pos = new TArrayF(fNumAnodes);
    fAngCorParams = new TArrayF(3);
}

// ----  Destructor ------------------------------------------------------------
R3BMusicHitPar::~R3BMusicHitPar()
{
    clear();
    if (fIn_use)
        delete fIn_use;
    if (fAnode_pos)
        delete fAnode_pos;
    if (fDetZHitParams)
        delete fDetZHitParams;
    if (fAngCorParams)
        delete fAngCorParams;
}

// ----  Method clear ----------------------------------------------------------
void R3BMusicHitPar::clear()
{
    status = kFALSE;
    resetInputVersions();
}

// ----  Method putParams ------------------------------------------------------
void R3BMusicHitPar::putParams(FairParamList* list)
{
    LOG(info) << "R3BMusicHitPar::putParams() called";
    if (!list)
    {
        return;
    }

    LOG(info) << "Nb Anodes: " << fNumAnodes;
    list->add("musicAnodeNumberPar", fNumAnodes);

    fAngCorParams->Set(3);
    list->add("musicAngCorHitPar", *fAngCorParams);

    fIn_use->Set(fNumAnodes);
    list->add("musicInUsePar", *fIn_use);

    fAnode_pos->Set(fNumAnodes);
    list->add("musicAnodePosPar", *fAnode_pos);

    LOG(info) << "Nb of parameters for charge-Z: " << fNumParamsZFit;
    list->add("musicZHitParamsFitPar", fNumParamsZFit);
    fDetZHitParams->Set(fNumParamsZFit);
    list->add("musicZHitPar", *fDetZHitParams);
}

// ----  Method getParams ------------------------------------------------------
Bool_t R3BMusicHitPar::getParams(FairParamList* list)
{
    LOG(info) << "R3BMusicHitPar::getParams() called";
    if (!list)
    {
        LOG(error) << "Could not initialize FairParamList";
        return kFALSE;
    }

    if (!list->fill("musicAnodeNumberPar", &fNumAnodes))
    {
        LOG(error) << "Could not initialize musicAnodeNumberPar";
        return kFALSE;
    }
    LOG(info) << "Nb Anodes: " << fNumAnodes;

    fAngCorParams->Set(3);
    if (!(list->fill("musicAngCorHitPar", fAngCorParams)))
    {
        LOG(warn) << "Could not initialize musicAngCorHitPar";
    }

    fIn_use->Set(fNumAnodes);
    if (!(list->fill("musicInUsePar", fIn_use)))
    {
        LOG(error) << "Could not initialize musicInUsePar";
        return kFALSE;
    }

    fAnode_pos->Set(fNumAnodes);
    if (!(list->fill("musicAnodePosPar", fAnode_pos)))
    {
        LOG(error) << "Could not initialize musicAnodePosPar";
        return kFALSE;
    }

    if (!list->fill("musicZHitParamsFitPar", &fNumParamsZFit))
    {
        LOG(error) << "Could not initialize musicZHitParamsFitPar";
        return kFALSE;
    }

    LOG(info) << "Nb of parameters for charge-Z: " << fNumParamsZFit;
    fDetZHitParams->Set(fNumParamsZFit);
    if (!(list->fill("musicZHitPar", fDetZHitParams)))
    {
        LOG(error) << "Could not initialize musicZHitPar";
        return kFALSE;
    }

    return kTRUE;
}

// ----  Method printParams ----------------------------------------------------
void R3BMusicHitPar::printParams()
{
    LOG(info) << "R3BMusicHitPar::Music detector Parameters";
    LOG(info) << "R3BMusicHitPar::Music anodes in use: ";
    for (Int_t j = 0; j < fNumAnodes; j++)
    {
        LOG(info) << "Anode " << j + 1 << " in use " << fIn_use->GetAt(j) << ", Position-Z: " << fAnode_pos->GetAt(j);
    }

    LOG(info) << "FitParam for charge-Z: " << fNumParamsZFit;
    for (Int_t j = 0; j < fNumParamsZFit; j++)
    {
        LOG(info) << "FitParam (" << j + 1 << ") = " << fDetZHitParams->GetAt(j);
    }
}

ClassImp(R3BMusicHitPar);
