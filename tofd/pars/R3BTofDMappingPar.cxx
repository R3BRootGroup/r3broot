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

// -------------------------------------------------------------
// -----        R3BTofDMappingPar source file              -----
// -----    Created 18/03/22 by J.L. Rodriguez-Sanchez     -----
// -------------------------------------------------------------

#include "R3BTofDMappingPar.h"

#include "FairParamList.h"
#include "R3BLogger.h"

#include "TMath.h"
#include "TString.h"

// ---- Standard Constructor ---------------------------------------------------
R3BTofDMappingPar::R3BTofDMappingPar(const TString& name, const TString& title, const TString& context)
    : FairParGenericSet(name, title, context)
    , fNumPlanes(4)
    , fNumPaddles(44)
    , fNumPmts(2)
{
    for (Int_t p = 0; p < fNumPmts; p++)
    {
        fTrigmap[p].resize(fNumPlanes);
    }
    for (Int_t plane = 0; plane < fNumPlanes; plane++)
        for (Int_t p = 0; p < fNumPmts; p++)
        {
            fTrigmap[p][plane] = new TArrayI(fNumPaddles);
            for (Int_t paddle = 0; paddle < fNumPaddles; paddle++)
                fTrigmap[p][plane]->AddAt(0, paddle);
        }
}

// ----  Destructor ------------------------------------------------------------
R3BTofDMappingPar::~R3BTofDMappingPar()
{
    clear();
    for (Int_t plane = 0; plane < fNumPlanes; plane++)
        for (Int_t p = 0; p < fNumPmts; p++)
        {
            if (fTrigmap[p][plane])
                delete fTrigmap[p][plane];
        }
}

// ----  Method clear ----------------------------------------------------------
void R3BTofDMappingPar::clear()
{
    status = kFALSE;
    resetInputVersions();
}

// ----  Method putParams ------------------------------------------------------
void R3BTofDMappingPar::putParams(FairParamList* list)
{
    R3BLOG(info, "called");
    if (!list)
    {
        R3BLOG(fatal, "FairParamList not found");
        return;
    }
    list->add("tofdPlanesPar", fNumPlanes);
    list->add("tofdPaddlesPar", fNumPaddles);

    R3BLOG(info, "Nb of planes: " << fNumPlanes);
    R3BLOG(info, "Nb of paddles: " << fNumPaddles);

    for (Int_t p = 0; p < fNumPmts; p++)
    {
        fTrigmap[p].resize(fNumPlanes);
    }
    char name[300];
    for (Int_t plane = 0; plane < fNumPlanes; plane++)
        for (Int_t p = 0; p < fNumPmts; p++)
        {
            fTrigmap[p][plane]->Set(fNumPaddles);
            sprintf(name, "tofdplane%dPmt%dPar", plane + 1, p + 1);
            list->add(name, *fTrigmap[p][plane]);
        }
}

// ----  Method getParams ------------------------------------------------------
Bool_t R3BTofDMappingPar::getParams(FairParamList* list)
{
    R3BLOG(info, "called");
    if (!list)
    {
        R3BLOG(fatal, "FairParamList not found");
        return kFALSE;
    }
    if (!list->fill("tofdPlanesPar", &fNumPlanes))
    {
        R3BLOG(info, "Could not initialize tofdPlanesPar");
        return kFALSE;
    }
    if (!list->fill("tofdPaddlesPar", &fNumPaddles))
    {
        R3BLOG(info, "Could not initialize tofdPaddlesPar");
        return kFALSE;
    }

    for (Int_t p = 0; p < fNumPmts; p++)
    {
        fTrigmap[p].resize(fNumPlanes);
    }
    char name[300];
    for (Int_t plane = 0; plane < fNumPlanes; plane++)
        for (Int_t p = 0; p < fNumPmts; p++)
        {
            fTrigmap[p][plane]->Set(fNumPaddles);
            sprintf(name, "tofdplane%dPmt%dPar", plane + 1, p + 1);

            if (!(list->fill(name, fTrigmap[p][plane])))
            {
                R3BLOG(error, "Could not initialize " << name);
                return kFALSE;
            }
        }

    return kTRUE;
}

// ----  Method print ----------------------------------------------------------
void R3BTofDMappingPar::print() { printParams(); }

// ----  Method printParams ----------------------------------------------------
void R3BTofDMappingPar::printParams()
{
    R3BLOG(info, "Mapping params for TofD: Num of planes: " << fNumPlanes << " and paddles: " << fNumPaddles);

    for (Int_t plane = 0; plane < fNumPlanes; plane++)
        for (Int_t p = 0; p < fNumPmts; p++)
            for (Int_t paddle = 0; paddle < fNumPaddles; paddle++)
            {
                R3BLOG(info,
                       "Plane: " << plane + 1 << ", pmt: " << p + 1 << ", paddle: " << paddle + 1
                                 << ", value: " << fTrigmap[p][plane]->GetAt(paddle));
            }
}

ClassImp(R3BTofDMappingPar);
