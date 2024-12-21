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

// ---------------------------------------------------------------------
// -----		      R3BMwpc1CalPar        	             -----
// -----          Created 15/10/19 by G. García Jiménez            -----
// -----  by modifying J.L. Rodriguez-Sanchez classes for Mwpc2    -----
// ---------------------------------------------------------------------

#include "R3BMwpc1CalPar.h"
#include "FairLogger.h"
#include "FairParamList.h"

#include "TMath.h"
#include "TString.h"

#include <iostream>

// ---- Standard Constructor ---------------------------------------------------
R3BMwpc1CalPar::R3BMwpc1CalPar(const char* name, const char* title, const char* context)
    : FairParGenericSet(name, title, context)
    , fNumPadsX(128)
    , fNumPadsY(40)
    , fNumParamsFit(1)
{
    fPadCalParams = new TArrayF((fNumPadsX + fNumPadsY) * fNumParamsFit);
}

// ----  Destructor ------------------------------------------------------------
R3BMwpc1CalPar::~R3BMwpc1CalPar()
{
    clear();
    delete fPadCalParams;
}

// ----  Method clear ----------------------------------------------------------
void R3BMwpc1CalPar::clear()
{
    status = kFALSE;
    resetInputVersions();
}

// ----  Method putParams ------------------------------------------------------
void R3BMwpc1CalPar::putParams(FairParamList* list)
{

    LOG(info) << "R3BMwpc1CalPar::putParams() called";
    if (!list)
    {
        return;
    }

    Int_t array_size = (fNumPadsX + fNumPadsY) * fNumParamsFit;
    LOG(info) << "Array Size: " << array_size;

    fPadCalParams->Set(array_size);

    list->add("mwpc1CalPar", *fPadCalParams);
    list->add("mwpc1PadXNumberPar", fNumPadsX);
    list->add("mwpc1PadYNumberPar", fNumPadsY);
    list->add("mwpc1ParamsFitPar", fNumParamsFit);
}

// ----  Method getParams ------------------------------------------------------
Bool_t R3BMwpc1CalPar::getParams(FairParamList* list)
{
    LOG(info) << "R3BMwpc1CalPar::getParams() called";
    if (!list)
    {
        return kFALSE;
    }

    if (!list->fill("mwpc1PadXNumberPar", &fNumPadsX))
    {
        return kFALSE;
    }
    if (!list->fill("mwpc1PadYNumberPar", &fNumPadsY))
    {
        return kFALSE;
    }
    if (!list->fill("mwpc1ParamsFitPar", &fNumParamsFit))
    {
        return kFALSE;
    }

    Int_t array_size = (fNumPadsX + fNumPadsY) * fNumParamsFit;
    LOG(info) << "Array Size: " << array_size;
    fPadCalParams->Set(array_size);

    if (!(list->fill("mwpc1CalPar", fPadCalParams)))
    {
        LOG(info) << "---Could not initialize mwpc1CalPar";
        return kFALSE;
    }

    return kTRUE;
}

// ----  Method printParams ----------------------------------------------------
void R3BMwpc1CalPar::printParams()
{
    LOG(info) << "R3BMwpc1CalPar: SOFIA MWPC1 Calibration Parameters: ";
    Int_t array_size = (fNumPadsX + fNumPadsY) * fNumParamsFit;

    for (Int_t i = 0; i < (fNumPadsX + fNumPadsY); i++)
    {
        if (i < fNumPadsX / 2)
        {
            LOG(info) << "MWPC1 Plane X down, Pad Number: " << i + 1;
        }
        else if (i < fNumPadsX)
        {
            LOG(info) << "MWPC1 Plane X up, Pad Number: " << i + 1 - fNumPadsX / 2;
        }
        else
        {
            LOG(info) << "MWPC1 Plane Y, Pad Number: " << i + 1 - fNumPadsX;
        }
        for (Int_t j = 0; j < fNumParamsFit; j++)
        {
            LOG(info) << "FitParam(" << j << ") = " << fPadCalParams->GetAt(i * fNumParamsFit + j);
        }
    }
}

ClassImp(R3BMwpc1CalPar);
