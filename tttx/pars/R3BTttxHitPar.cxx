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
// -----                  R3BTttxHitPar source file             -----
// -----         Created 18/01/24 by J.L. Rodriguez-Sanchez     -----
// ------------------------------------------------------------------

#include <FairLogger.h>
#include <FairParamList.h>

#include "R3BLogger.h"
#include "R3BTttxHitPar.h"

#include <TArrayF.h>
#include <TMath.h>
#include <TString.h>

// ---- Standard Constructor ---------------------------------------------------
R3BTttxHitPar::R3BTttxHitPar(const char* name, const char* title, const char* context)
    : FairParGenericSet(name, title, context)
    , fNumDets(2)
{
    detName = "TttxHit";
    fZfitpar = new TArrayF(fNumDets * fNumParsZfit);
}

// ----  Destructor ------------------------------------------------------------
R3BTttxHitPar::~R3BTttxHitPar()
{
    this->clear(); // NOLINT
}

// ----  Method clear ----------------------------------------------------------
void R3BTttxHitPar::clear()
{
    status = kFALSE;
    resetInputVersions();
}

// ----  Method putParams ------------------------------------------------------
void R3BTttxHitPar::putParams(FairParamList* list)
{
    R3BLOG(info, "called");
    if (!list)
    {
        R3BLOG(error, "Could not initialize FairParamList");
        return;
    }

    list->add("tttxDetNumberPar", fNumDets);
    list->add("tttxNumberParsZFit", fNumParsZfit);

    auto array_size = fNumDets * fNumParsZfit;
    LOG(info) << "Size for array fZfitpar: " << array_size;
    fZfitpar->Set(array_size);
    list->add("tttxZfitPar", *fZfitpar);
}

// ----  Method getParams ------------------------------------------------------
Bool_t R3BTttxHitPar::getParams(FairParamList* list)
{
    R3BLOG(info, "called");
    if (!list)
    {
        R3BLOG(error, "Could not initialize FairParamList");
        return kFALSE;
    }

    if (!list->fill("tttxDetNumberPar", &fNumDets))
    {
        R3BLOG(fatal, "Could not initialize tttxDetNumberPar");
        return kFALSE;
    }

    if (!list->fill("tttxNumberParsZFit", &fNumParsZfit))
    {
        R3BLOG(fatal, "Could not initialize tttxNumberParsZFit");
        return kFALSE;
    }

    auto array_size = fNumDets * fNumParsZfit;
    R3BLOG(info, "Total number of Z-fit parameters: " << array_size);

    fZfitpar->Set(array_size);
    if (!(list->fill("tttxZfitPar", fZfitpar)))
    {
        R3BLOG(fatal, "Could not initialize tttxZfitPar");
        return kFALSE;
    }
    return kTRUE;
}

// ----  Method print ----------------------------------------------------------
void R3BTttxHitPar::print()
{
    R3BLOG(info, "Number of detectors: " << fNumDets);
    R3BLOG(info, "Number of Z-fit parameters: " << fNumParsZfit);
    R3BLOG(info, "TTTX Hit Parameters");
    for (int d = 0; d < fNumDets; d++)
    {
        R3BLOG(info, "Detector: " << d + 1);
        for (Int_t j = 0; j < fNumParsZfit; j++)
        {
            LOG(info) << "FitParam(" << j + 1 << ") = " << fZfitpar->GetAt(d * fNumParsZfit + j);
        }
    }
}

ClassImp(R3BTttxHitPar)
