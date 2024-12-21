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
// -----              R3BTttxStripCalPar source file            -----
// -----         Created 18/01/24 by J.L. Rodriguez-Sanchez     -----
// ------------------------------------------------------------------

#include <FairLogger.h>
#include <FairParamList.h>

#include "R3BLogger.h"
#include "R3BTttxStripCalPar.h"

#include <TMath.h>
#include <TString.h>

// ---- Standard Constructor ---------------------------------------------------
R3BTttxStripCalPar::R3BTttxStripCalPar(const char* name, const char* title, const char* context)
    : FairParGenericSet(name, title, context)
{
    detName = "TttxCal";
    fStripCalParams = new TArrayF(fNumDets * fNumStrips * fNumParsFit);
    fIn_use = new TArrayI(fNumDets * fNumStrips);
}

// ----  Destructor ------------------------------------------------------------
R3BTttxStripCalPar::~R3BTttxStripCalPar()
{
    this->clear(); // NOLINT
    delete fStripCalParams;
}

// ----  Method clear ----------------------------------------------------------
void R3BTttxStripCalPar::clear()
{
    status = kFALSE;
    resetInputVersions();
}

// ----  Method putParams ------------------------------------------------------
void R3BTttxStripCalPar::putParams(FairParamList* list)
{
    R3BLOG(info, "called");
    if (!list)
    {
        R3BLOG(error, "Could not initialize FairParamList");
        return;
    }

    auto array_size = fNumDets * fNumStrips * fNumParsFit;
    LOG(info) << "Size for array fStripCalParams: " << array_size;
    fStripCalParams->Set(array_size);

    auto array_size_use = fNumDets * fNumStrips;
    LOG(info) << "Size for array fIn_use: " << array_size_use;
    fIn_use->Set(array_size_use);

    list->add("tttxDetNumberPar", fNumDets);
    list->add("tttxStripNumberPar", fNumStrips);
    list->add("tttxNumberParsFit", fNumParsFit);
    list->add("tttxInUsePar", *fIn_use);
    list->add("tttxStripCalPar", *fStripCalParams);
}

// ----  Method getParams ------------------------------------------------------
Bool_t R3BTttxStripCalPar::getParams(FairParamList* list)
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

    if (!list->fill("tttxStripNumberPar", &fNumStrips))
    {
        R3BLOG(fatal, "Could not initialize tttxStripNumberPar");
        return kFALSE;
    }

    if (!list->fill("tttxNumberParsFit", &fNumParsFit))
    {
        R3BLOG(fatal, "Could not initialize tttxNumberParsFit");
        return kFALSE;
    }

    auto array_size = fNumDets * fNumStrips;
    R3BLOG(info, "Total number of strips: " << array_size);

    fIn_use->Set(array_size);
    if (!(list->fill("tttxInUsePar", fIn_use)))
    {
        R3BLOG(fatal, "Could not initialize tttxInUsePar");
        return kFALSE;
    }

    fStripCalParams->Set(array_size * fNumParsFit);
    if (!(list->fill("tttxStripCalPar", fStripCalParams)))
    {
        R3BLOG(fatal, "Could not initialize tttxStripCalPar");
        return kFALSE;
    }

    return kTRUE;
}

// ----  Method print ----------------------------------------------------------
void R3BTttxStripCalPar::print()
{
    R3BLOG(info, "TTTX strip Parameters");

    for (int d = 0; d < fNumDets; d++)
    {
        R3BLOG(info, "Tttx detector: " << d + 1);
        for (Int_t i = 0; i < fNumStrips; i++)
        {
            R3BLOG(info, "Tttx strip number: " << i + 1);
            for (Int_t j = 0; j < fNumParsFit; j++)
            {
                LOG(info) << "FitParam(" << j + 1
                          << ") = " << fStripCalParams->GetAt(d * fNumParsFit * fNumStrips + i * fNumParsFit + j);
            }
        }
    }
}

ClassImp(R3BTttxStripCalPar)
