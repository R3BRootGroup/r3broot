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

#include "R3BCalifaTotCalPar.h"

#include "FairLogger.h"
#include "FairParamList.h"

#include "TMath.h"
#include "TString.h"

#include <iostream>

// ---- Standard Constructor ---------------------------------------------------
R3BCalifaTotCalPar::R3BCalifaTotCalPar(const char* name, const char* title, const char* context)
    : FairParGenericSet(name, title, context)
    , fNumCrystals(4864)
    , fNumParamsFit(2) // pol1
{
    fCryCalParams = new TArrayF(fNumCrystals * fNumParamsFit);
}

// ----  Destructor ------------------------------------------------------------
R3BCalifaTotCalPar::~R3BCalifaTotCalPar()
{
    clear();
    if (fCryCalParams)
        delete fCryCalParams;
}

void R3BCalifaTotCalPar::clear()
{
    status = kFALSE;
    resetInputVersions();
}

// ----  Method putParams ------------------------------------------------------
void R3BCalifaTotCalPar::putParams(FairParamList* list)
{
    LOG(info) << "R3BCalifaTotCalPar::putParams() called";
    if (!list)
    {
        return;
    }

    Int_t array_size = fNumCrystals * fNumParamsFit;
    LOG(info) << "Array Size: " << array_size;

    fCryCalParams->Set(array_size);

    list->add("califaTotCalPar", *fCryCalParams);
    list->add("califaTotNumberPar", fNumCrystals);
    list->add("califaTotParamsFitPar", fNumParamsFit);
}

// ----  Method getParams ------------------------------------------------------
Bool_t R3BCalifaTotCalPar::getParams(FairParamList* list)
{
    LOG(info) << "R3BCalifaTotCalPar::getParams() called";
    if (!list)
    {
        return kFALSE;
    }

    if (!list->fill("califaTotNumberPar", &fNumCrystals))
    {
        return kFALSE;
    }

    if (!list->fill("califaTotParamsFitPar", &fNumParamsFit))
    {
        return kFALSE;
    }

    Int_t array_size = fNumCrystals * fNumParamsFit;
    LOG(info) << "Array Size: " << array_size;
    fCryCalParams->Set(array_size);

    if (!(list->fill("califaTotCalPar", fCryCalParams)))
    {
        LOG(info) << "---Could not initialize califaTotCalPar";
        return kFALSE;
    }

    return kTRUE;
}

// ----  Method print ----------------------------------------------------------
void R3BCalifaTotCalPar::print() { printParams(); }

// ----  Method printParams ----------------------------------------------------
void R3BCalifaTotCalPar::printParams()
{
    LOG(info) << "R3BCalifaTotCalPar: Califa Crystal Calibration Parameters: ";
    Int_t array_size = fNumCrystals * fNumParamsFit;

    for (Int_t i = 0; i < fNumCrystals; i++)
    {
        LOG(info) << "Califa Crystal number: " << i + 1;
        for (Int_t j = 0; j < fNumParamsFit; j++)
        {
            LOG(info) << "FitParam(" << j << ") = " << fCryCalParams->GetAt(i * fNumParamsFit + j);
        }
    }
}

ClassImp(R3BCalifaTotCalPar)
