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

#include "R3BCalifaCrystalCalPar.h"
#include "R3BLogger.h"

#include <FairLogger.h>
#include <FairParamList.h>

#include <TMath.h>
#include <TString.h>
#include <iostream>

// ---- Standard Constructor ---------------------------------------------------
R3BCalifaCrystalCalPar::R3BCalifaCrystalCalPar(const char* name, const char* title, const char* context)
    : FairParGenericSet(name, title, context)
{
    fCryCalParams = new TArrayF(fNumCrystals * fNumParamsFit);
}

// ----  Destructor ------------------------------------------------------------
R3BCalifaCrystalCalPar::~R3BCalifaCrystalCalPar()
{
    if (fCryCalParams)
        delete fCryCalParams;
}

// ----  Method clear ----------------------------------------------------------
void R3BCalifaCrystalCalPar::clear()
{
    status = kFALSE;
    resetInputVersions();
}

// ----  Method putParams ------------------------------------------------------
void R3BCalifaCrystalCalPar::putParams(FairParamList* list)
{
    LOG(info) << "R3BCalifaCrystalCalPar::putParams() called";
    if (!list)
    {
        return;
    }

    Int_t array_size = fNumCrystals * fNumParamsFit;
    LOG(info) << "Array Size: " << array_size;

    fCryCalParams->Set(array_size);

    list->add("califaCrystalCalPar", *fCryCalParams);
    list->add("califaCrystalNumberPar", fNumCrystals);
    list->add("califaCrystalParamsFitPar", fNumParamsFit);
}

// ----  Method getParams ------------------------------------------------------
Bool_t R3BCalifaCrystalCalPar::getParams(FairParamList* list)
{
    LOG(info) << "R3BCalifaCrystalCalPar::getParams() called";
    if (!list)
    {
        return kFALSE;
    }

    if (!list->fill("califaCrystalNumberPar", &fNumCrystals))
    {
        return kFALSE;
    }

    if (!list->fill("califaCrystalParamsFitPar", &fNumParamsFit))
    {
        return kFALSE;
    }

    Int_t array_size = fNumCrystals * fNumParamsFit;
    LOG(info) << "Nb_crystals * Nb_pars: " << array_size;
    fCryCalParams->Set(array_size);

    if (!(list->fill("califaCrystalCalPar", fCryCalParams)))
    {
        LOG(info) << "---Could not initialize califaCrystalCalPar";
        return kFALSE;
    }

    return kTRUE;
}

// ----  Method print ----------------------------------------------------------
void R3BCalifaCrystalCalPar::print() { printParams(); }

// ----  Method printCalCrystalInfo --------------------------------------------
void R3BCalifaCrystalCalPar::printCalCrystalInfo(const int cryID)
{
    assert(std::clamp(cryID, 1, fNumCrystals) == cryID && "crystal_id must be given in 1-base");
    R3BLOG_IF(fatal, cryID > fNumCrystals, "crystal_id does not exist, crystal_id<=" << fNumCrystals);

    auto index = cryID - 1;
    LOG(info) << "Califa Crystal number: " << cryID;
    for (Int_t j = 0; j < fNumParamsFit; j++)
    {
        LOG(info) << "FitParam(" << j + 1 << ") = " << fCryCalParams->GetAt(index * fNumParamsFit + j);
    }
}

// ----  Method printParams ----------------------------------------------------
void R3BCalifaCrystalCalPar::printParams()
{
    LOG(info) << "R3BCalifaCrystalCalPar::Califa Crystal Calibration Parameters: ";
    Int_t array_size = fNumCrystals * fNumParamsFit;

    for (Int_t i = 0; i < fNumCrystals; i++)
    {
        LOG(info) << "Califa Crystal number: " << i + 1;
        for (Int_t j = 0; j < fNumParamsFit; j++)
        {
            LOG(info) << "FitParam(" << j + 1 << ") = " << fCryCalParams->GetAt(i * fNumParamsFit + j);
        }
    }
}

ClassImp(R3BCalifaCrystalCalPar)
