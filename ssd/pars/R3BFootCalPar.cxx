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
// -----              R3BFootCalPar source file                 -----
// -----         Created 07/09/21 by J.L. Rodriguez-Sanchez     -----
// ------------------------------------------------------------------

#include "R3BFootCalPar.h"
#include "R3BLogger.h"

#include <FairLogger.h>
#include <FairParamList.h>

#include <TArrayF.h>
#include <TMath.h>
#include <TString.h>

// ---- Standard Constructor ---------------------------------------------------
R3BFootCalPar::R3BFootCalPar(const char* name, const char* title, const char* context)
    : FairParGenericSet(name, title, context)
{
    detName = "FootCal";
    fStripCalParams = new TArrayF(fNumDets * fNumStrips * fNumParsFit);
}

// ----  Destructor ------------------------------------------------------------
R3BFootCalPar::~R3BFootCalPar()
{
    clear();
    delete fStripCalParams;
}

// ----  Method clear ----------------------------------------------------------
void R3BFootCalPar::clear()
{
    status = kFALSE;
    resetInputVersions();
}

// ----  Method putParams ------------------------------------------------------
void R3BFootCalPar::putParams(FairParamList* list)
{
    R3BLOG(info, "called");
    if (!list)
    {
        R3BLOG(error, "Could not initialize FairParamList");
        return;
    }

    Int_t array_size = fNumDets * fNumStrips;
    LOG(info) << "Array Size: " << array_size;

    fStripCalParams->Set(array_size);

    list->add("footDetNumberPar", fNumDets);
    list->add("footStripNumberPar", fNumStrips);
    list->add("footNumberParsFit", fNumParsFit);
    list->add("footStripCalPar", *fStripCalParams);
}

// ----  Method getParams ------------------------------------------------------
Bool_t R3BFootCalPar::getParams(FairParamList* list)
{
    R3BLOG(info, "called");
    if (!list)
    {
        R3BLOG(error, "Could not initialize FairParamList");
        return kFALSE;
    }

    if (!list->fill("footDetNumberPar", &fNumDets))
    {
        LOG(fatal) << "R3BFootCalPar::Could not initialize footDetNumberPar";
        return kFALSE;
    }

    if (!list->fill("footStripNumberPar", &fNumStrips))
    {
        LOG(fatal) << "R3BFootCalPar::Could not initialize footStripNumberPar";
        return kFALSE;
    }

    if (!list->fill("footNumberParsFit", &fNumParsFit))
    {
        LOG(fatal) << "R3BFootCalPar::Could not initialize footNumberParsFit";
        return kFALSE;
    }

    Int_t array_size = fNumDets * fNumStrips;
    LOG(info) << "Total number of strips: " << array_size;
    fStripCalParams->Set(array_size * fNumParsFit);

    if (!(list->fill("footStripCalPar", fStripCalParams)))
    {
        LOG(fatal) << "R3BFootCalPar::Could not initialize footStripCalPar";
        return kFALSE;
    }

    return kTRUE;
}

// ----  Method print ----------------------------------------------------------
void R3BFootCalPar::print()
{
    R3BLOG(info, "Foot strip Parameters");

    for (Int_t d = 0; d < fNumDets; d++)
    {
        LOG(info) << "Foot detector number: " << d + 1;
        for (Int_t i = 0; i < fNumStrips; i++)
        {
            LOG(info) << "Foot strip number: " << i + 1;
            for (Int_t j = 0; j < fNumParsFit; j++)
            {
                LOG(info) << "FitParam(" << j + 1
                          << ") = " << fStripCalParams->GetAt(d * fNumParsFit * fNumStrips + i * fNumParsFit + j);
            }
        }
    }
}

ClassImp(R3BFootCalPar)
