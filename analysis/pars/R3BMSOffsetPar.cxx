/******************************************************************************
 *   Copyright (C) 2019 GSI Helmholtzzentrum für Schwerionenforschung GmbH    *
 *   Copyright (C) 2019-2024 Members of R3B Collaboration                     *
 *                                                                            *
 *             This software is distributed under the terms of the            *
 *                 GNU General Public Licence (GPL) version 3,                *
 *                    copied verbatim in the file "LICENSE".                  *
 *                                                                            *
 * In applying this license GSI does not waive the privileges and immunities  *
 * granted to it by virtue of its status as an Intergovernmental Organization *
 * or submit itself to any jurisdiction.                                      *
 ******************************************************************************/

#include "R3BMSOffsetPar.h"
#include "R3BLogger.h"

#include "FairDetParIo.h"
#include "FairLogger.h"
#include "FairParamList.h"

#include "TArrayD.h"
#include "TMath.h"
#include "TString.h"

// ---- Standard Constructor ---------------------------------------------------
// R3BMSOffsetPar::R3BMSOffsetPar(const char* name, const char* title, const char* context)
//    : FairParGenericSet(name, title, context)
R3BMSOffsetPar::R3BMSOffsetPar(std::string_view name, std::string_view title, std::string_view context)
    : FairParGenericSet(name.data(), title.data(), context.data())
    , fMSOffset(0)
{
}

// ----  Method clear ----------------------------------------------------------
void R3BMSOffsetPar::clear()
{
    status = kFALSE;
    resetInputVersions();
}

// ----  Method putParams ------------------------------------------------------
void R3BMSOffsetPar::putParams(FairParamList* list)
{
    R3BLOG(info, "called");
    if (list == nullptr)
    {
        R3BLOG(fatal, "Could not find FairParamList");
        return;
    }

    list->add("MSOffsetPar", fMSOffset);
}

// ----  Method getParams ------------------------------------------------------
Bool_t R3BMSOffsetPar::getParams(FairParamList* list)
{
    R3BLOG(info, "called");
    if (list == nullptr)
    {
        R3BLOG(fatal, "Could not find FairParamList");
        return kFALSE;
    }

    if (!list->fill("MSOffsetPar", &fMSOffset))
    {
        LOG(error) << "Could not initialize MSOffsetPar";
        return kFALSE;
    }

    return kTRUE;
}

// ----  Method print ----------------------------------------------------------
void R3BMSOffsetPar::print() { printParams(); }

// ----  Method printParams ----------------------------------------------------
void R3BMSOffsetPar::printParams() { R3BLOG(info, "Master Start Offset: " << fMSOffset); }

ClassImp(R3BMSOffsetPar); // NOLINT
