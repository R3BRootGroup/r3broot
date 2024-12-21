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

// --------------------------------------------------------
// -----		 R3BMwpc3CalPar                -----
// -----      Created 11/10/2019 by G.García Jiménez  -----
// --------------------------------------------------------

#include "R3BMwpc3CalPar.h"
#include "FairLogger.h"
#include "FairParamList.h"
#include "TMath.h"
#include "TString.h"

/* ---- Standard Constructor ---- */
R3BMwpc3CalPar::R3BMwpc3CalPar(const char* name, const char* title, const char* context)
    : FairParGenericSet(name, title, context)
    , fNumPadsX(288)
    , fNumPadsY(120)
    , fNumParamsFit(2)
{
    fPadCalParams = new TArrayF((fNumPadsX + fNumPadsY) * fNumParamsFit);
}

/* ----  Destructor ---- */
R3BMwpc3CalPar::~R3BMwpc3CalPar()
{
    clear();
    delete fPadCalParams;
}

/* ----  Method clear ---- */
void R3BMwpc3CalPar::clear()
{
    status = kFALSE;
    resetInputVersions();
}

/* ----  Method putParams ---- */
void R3BMwpc3CalPar::putParams(FairParamList* list)
{

    LOG(info) << "R3BMwpc3CalPar::putParams() called";
    if (!list)
    {
        return;
    }

    Int_t array_size = (fNumPadsX + fNumPadsY) * fNumParamsFit;
    LOG(info) << "Array Size: " << array_size;

    fPadCalParams->Set(array_size);

    list->add("mwpc3CalPar", *fPadCalParams);
    list->add("mwpc3PadXNumberPar", fNumPadsX);
    list->add("mwpc3PadYNumberPar", fNumPadsY);
    list->add("mwpc3ParamsFitPar", fNumParamsFit);
}

/* ----  Method getParams ---- */
Bool_t R3BMwpc3CalPar::getParams(FairParamList* list)
{
    LOG(info) << "R3BMwpc3CalPar::getParams() called";
    if (!list)
    {
        return kFALSE;
    }

    if (!list->fill("mwpc3PadXNumberPar", &fNumPadsX))
    {
        return kFALSE;
    }
    if (!list->fill("mwpc3PadYNumberPar", &fNumPadsY))
    {
        return kFALSE;
    }
    if (!list->fill("mwpc3ParamsFitPar", &fNumParamsFit))
    {
        return kFALSE;
    }

    Int_t array_size = (fNumPadsX + fNumPadsY) * fNumParamsFit;
    LOG(info) << "Array Size: " << array_size;
    fPadCalParams->Set(array_size);

    if (!(list->fill("mwpc3CalPar", fPadCalParams)))
    {
        LOG(info) << "---Could not initialize mwpc3CalPar";
        return kFALSE;
    }

    return kTRUE;
}

/* ----  Method printParams ---- */
void R3BMwpc3CalPar::printParams()
{
    LOG(info) << "R3BMwpc3CalPar: SOFIA MWPC3 Calibration Parameters: ";
    Int_t array_size = (fNumPadsX + fNumPadsY) * fNumParamsFit;
    for (Int_t i = 0; i < (fNumPadsX + fNumPadsY); i++)
    {
        if (i < fNumPadsX)
        {
            LOG(info) << "MWPC3 Plane X Pad Number: " << i + 1;
        }
        else
        {
            LOG(info) << "MWPC3 Plane Y Pad Number: " << i + 1 - fNumPadsX;
        }
        for (Int_t j = 0; j < fNumParamsFit; j++)
        {
            LOG(info) << "FitParam(" << j << ") = " << fPadCalParams->GetAt(i * fNumParamsFit + j);
        }
    }
}

ClassImp(R3BMwpc3CalPar);
