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
// -----         R3BAmsStripCalPar source file                  -----
// -----         Created 29/05/18  by J.L. Rodriguez-Sanchez    -----
// ------------------------------------------------------------------

#include "R3BAmsStripCalPar.h"

#include "FairLogger.h"
#include "FairParamList.h"

#include "TArrayF.h"
#include "TMath.h"
#include "TString.h"

#include <iostream>

// ---- Standard Constructor ---------------------------------------------------
R3BAmsStripCalPar::R3BAmsStripCalPar(const char* name, const char* title, const char* context)
    : FairParGenericSet(name, title, context)
{
    fStripCalParams = new TArrayF(3072); // 1024strips*3CalibrationParameters
    fNumDets = 1;
    fNumStrips = 1024;
    fNumStripsS = 640;
    fNumStripsK = 384;
    fNumParamsFit = 3; // gaus
}

// ----  Destructor ------------------------------------------------------------
R3BAmsStripCalPar::~R3BAmsStripCalPar()
{
    clear();
    delete fStripCalParams;
}

// ----  Method clear ----------------------------------------------------------
void R3BAmsStripCalPar::clear()
{
    status = kFALSE;
    resetInputVersions();
}

// ----  Method putParams ------------------------------------------------------
void R3BAmsStripCalPar::putParams(FairParamList* list)
{
    LOG(info) << "R3BAmsStripCalPar::putParams() called";
    if (!list)
    {
        return;
    }

    Int_t array_size = fNumDets * fNumStrips * fNumParamsFit;
    LOG(info) << "Array Size: " << array_size;

    fStripCalParams->Set(array_size);

    list->add("amsStripCalPar", *fStripCalParams);
    list->add("amsDetNumberPar", fNumDets);
    list->add("amsStripNumberPar", fNumStrips);
    list->add("amsStripSNumberPar", fNumStripsS);
    list->add("amsStripKNumberPar", fNumStripsK);
    list->add("amsStripParamsFitPar", fNumParamsFit);
}

// ----  Method getParams ------------------------------------------------------
Bool_t R3BAmsStripCalPar::getParams(FairParamList* list)
{
    LOG(info) << "R3BAmsStripCalPar::getParams() called";
    if (!list)
    {
        return kFALSE;
    }

    if (!list->fill("amsDetNumberPar", &fNumDets))
    {
        return kFALSE;
    }

    if (!list->fill("amsStripNumberPar", &fNumStrips))
    {
        return kFALSE;
    }

    if (!list->fill("amsStripSNumberPar", &fNumStripsS))
    {
        return kFALSE;
    }

    if (!list->fill("amsStripKNumberPar", &fNumStripsK))
    {
        return kFALSE;
    }

    if (!list->fill("amsStripParamsFitPar", &fNumParamsFit))
    {
        return kFALSE;
    }

    Int_t array_size = fNumDets * fNumStrips * fNumParamsFit;
    LOG(info) << "Array Size: " << array_size;
    fStripCalParams->Set(array_size);

    if (!(list->fill("amsStripCalPar", fStripCalParams)))
    {
        LOG(info) << "---Could not initialize amsStripCalPar";
        return kFALSE;
    }

    return kTRUE;
}

// ----  Method print ----------------------------------------------------------
void R3BAmsStripCalPar::print() { printParams(); }

// ----  Method printParams ----------------------------------------------------
void R3BAmsStripCalPar::printParams()
{
    LOG(info) << "R3BAmsStripCalPar: ams strip Parameters: ";
    Int_t array_size = fNumDets * fNumStrips * fNumParamsFit;

    for (Int_t d = 0; d < fNumDets; d++)
    {
        LOG(info) << "AMS detector number: " << d;
        for (Int_t i = 0; i < fNumStrips; i++)
        {
            LOG(info) << "AMS Strip number: " << i;
            for (Int_t j = 0; j < fNumParamsFit; j++)
            {
                LOG(info) << "FitParam(" << j
                          << ") = " << fStripCalParams->GetAt(d * fNumParamsFit * fNumStrips + i * fNumParamsFit + j);
            }
        }
    }
}

ClassImp(R3BAmsStripCalPar)
