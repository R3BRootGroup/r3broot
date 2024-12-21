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

// --------------------------------------------------------------------
// -----		 R3BMwpc0CalPar 			    -----
// ----- 	Created 7/10/2019 by G.García Jiménez		    -----
// --------------------------------------------------------------------

#include "R3BMwpc0CalPar.h"
#include "FairLogger.h"
#include "FairParamList.h"

#include "TMath.h"
#include "TString.h"

#include <iostream>

// ---- Standard Constructor ---------------------------------------------------
R3BMwpc0CalPar::R3BMwpc0CalPar(const char* name, const char* title, const char* context)
    : FairParGenericSet(name, title, context)
    , fNumPadsX(64)
    , fNumPadsY(64)
    , fNumParamsFit(1)
{
    fPadCalParams = new TArrayF((fNumPadsX + fNumPadsY) *
                                fNumParamsFit); //(64 Pads in X + 64 Pads in Y) x 1 Calibration Parameter (Pedestal)
}

// ----  Destructor ------------------------------------------------------------
R3BMwpc0CalPar::~R3BMwpc0CalPar()
{
    clear();
    if (fPadCalParams)
        delete fPadCalParams;
}

// ----  Method clear ----------------------------------------------------------
void R3BMwpc0CalPar::clear()
{
    status = kFALSE;
    resetInputVersions();
}

// ----  Method putParams ------------------------------------------------------
void R3BMwpc0CalPar::putParams(FairParamList* list)
{

    LOG(info) << "R3BMwpc0CalPar::putParams() called";
    if (!list)
    {
        return;
    }

    Int_t array_size = (fNumPadsX + fNumPadsY) * fNumParamsFit;
    LOG(info) << "Array Size: " << array_size;

    fPadCalParams->Set(array_size);

    list->add("mwpc0CalPar", *fPadCalParams);
    list->add("mwpc0PadXNumberPar", fNumPadsX);
    list->add("mwpc0PadYNumberPar", fNumPadsY);
    list->add("mwpc0ParamsFitPar", fNumParamsFit);
}

// ----  Method getParams ------------------------------------------------------
Bool_t R3BMwpc0CalPar::getParams(FairParamList* list)
{
    LOG(info) << "R3BMwpc0CalPar::getParams() called";
    if (!list)
    {
        return kFALSE;
    }

    if (!list->fill("mwpc0PadXNumberPar", &fNumPadsX))
    {
        return kFALSE;
    }
    if (!list->fill("mwpc0PadYNumberPar", &fNumPadsY))
    {
        return kFALSE;
    }
    if (!list->fill("mwpc0ParamsFitPar", &fNumParamsFit))
    {
        return kFALSE;
    }

    Int_t array_size = (fNumPadsX + fNumPadsY) * fNumParamsFit;
    LOG(info) << "Array Size: " << array_size;
    fPadCalParams->Set(array_size);

    if (!(list->fill("mwpc0CalPar", fPadCalParams)))
    {
        LOG(info) << "---Could not initialize mwpc0CalPar";
        return kFALSE;
    }

    return kTRUE;
}

// ----  Method printParams ----------------------------------------------------
void R3BMwpc0CalPar::printParams()
{
    LOG(info) << "R3BMwpc0CalPar: SOFIA MWPC0 Calibration Parameters: ";
    Int_t array_size = (fNumPadsX + fNumPadsY) * fNumParamsFit;

    for (Int_t i = 0; i < (fNumPadsX + fNumPadsY); i++)
    {
        if (i < fNumPadsX)
        {
            LOG(info) << "MWPC0 Plane X Pad Number: " << i + 1;
        }
        else
        {
            LOG(info) << "MWPC0 Plane Y Pad Number: " << i + 1 - fNumPadsX;
        }
        for (Int_t j = 0; j < fNumParamsFit; j++)
        {
            LOG(info) << "FitParam(" << j << ") = " << fPadCalParams->GetAt(i * fNumParamsFit + j);
        }
    }
}

ClassImp(R3BMwpc0CalPar)
