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

#include "R3BMusliHitPar.h"
#include "R3BLogger.h"

#include "FairDetParIo.h"
#include "FairLogger.h"
#include "FairParamList.h"

#include "TArrayD.h"
#include "TMath.h"
#include "TString.h"

// ---- Standard Constructor ---------------------------------------------------
R3BMusliHitPar::R3BMusliHitPar(const char* name, const char* title, const char* context)
    : FairParGenericSet(name, title, context)
    , fNumGroups(15)
    , fNumTypes(4)
    , fMaxMult(20)
    , fNumParamsEcorrBetaFit(2)
    , fNumParamsZFit(3)
{
    fEaveVsBetaHitParams = new TArrayD(fNumTypes);
    fEcorrBetaHitParams = new TArrayD(fNumTypes * fNumParamsEcorrBetaFit);
    fZHitParams = new TArrayD(fNumTypes * fNumParamsZFit);
}

// ----  Destructor ------------------------------------------------------------
R3BMusliHitPar::~R3BMusliHitPar() { clear(); }

// ----  Method clear ----------------------------------------------------------
void R3BMusliHitPar::clear()
{
    status = kFALSE;
    resetInputVersions();
}

// ----  Method putParams ------------------------------------------------------
void R3BMusliHitPar::putParams(FairParamList* list)
{
    R3BLOG(info, "called");
    if (!list)
    {
        R3BLOG(fatal, "Could not find FairParamList");
        return;
    }

    list->add("musliNumGroupsAnodesPar", fNumGroups);
    list->add("musliNumTypesPar", fNumTypes);
    list->add("musliMaxMultPar", fMaxMult);
    list->add("musliNumParamsEcorrBetaFitPar", fNumParamsEcorrBetaFit);
    list->add("musliNumParamsZFitPar", fNumParamsZFit);

    fEaveVsBetaHitParams->Set(fNumTypes);
    list->add("musliEaveVsBetaHitPar", fEaveVsBetaHitParams);

    fEaveVsBetaHitParams->Set(fNumTypes * fNumParamsEcorrBetaFit);
    list->add("musliEcorrBetaHitPar", fEcorrBetaHitParams);

    fZHitParams->Set(fNumTypes * fNumParamsZFit);
    list->add("musliZHitPar", fZHitParams);
}

// ----  Method getParams ------------------------------------------------------
Bool_t R3BMusliHitPar::getParams(FairParamList* list)
{
    R3BLOG(info, "called");
    if (!list)
    {
        R3BLOG(fatal, "Could not find FairParamList");
        return kFALSE;
    }

    if (!list->fill("musliNumGroupsAnodesPar", &fNumGroups))
    {
        LOG(error) << "Could not initialize musliNumTypesPar";
        return kFALSE;
    }

    if (!list->fill("musliNumTypesPar", &fNumTypes))
    {
        LOG(error) << "Could not initialize musliNumTypesPar";
        return kFALSE;
    }

    if (!list->fill("musliMaxMultPar", &fMaxMult))
    {
        LOG(error) << "Could not initialize musliMaxMultPar";
        return kFALSE;
    }

    if (!list->fill("musliNumParamsEcorrBetaFitPar", &fNumParamsEcorrBetaFit))
    {
        LOG(error) << "Could not initialize musliNumParamsEcorrBetaFit";
        return kFALSE;
    }

    if (!list->fill("musliNumParamsZFitPar", &fNumParamsZFit))
    {
        LOG(error) << "Could not initialize musliNumParamsZFit";
        return kFALSE;
    }

    fEaveVsBetaHitParams->Set(fNumTypes);
    if (!list->fill("musliEaveVsBetaHitPar", fEaveVsBetaHitParams))
    {
        LOG(error) << "Could not initialize musliEaveVsBetaHitPar";
        return kFALSE;
    }

    fEcorrBetaHitParams->Set(fNumTypes * fNumParamsEcorrBetaFit);
    if (!list->fill("musliEcorrBetaHitPar", fEcorrBetaHitParams))
    {
        LOG(error) << "Could not initialize musliEcorrBetaHitPar";
        return kFALSE;
    }

    fZHitParams->Set(fNumTypes * fNumParamsZFit);
    if (!list->fill("musliZHitPar", fZHitParams))
    {
        LOG(error) << "Could not initialize musliZHitPar";
        return kFALSE;
    }

    return kTRUE;
}

// ----  Method print ----------------------------------------------------------
void R3BMusliHitPar::print() { printParams(); }

// ----  Method printParams ----------------------------------------------------
void R3BMusliHitPar::printParams()
{

    R3BLOG(info, "Nb of musli groups anodes signals: " << fNumGroups);
    R3BLOG(info, "Nb of musli data types: " << fNumTypes);
    R3BLOG(info, "Nb of musli data types: " << fNumTypes);
    R3BLOG(info, "Nb of parameters in the EcorrBeta Fit: " << fNumParamsEcorrBetaFit);
    R3BLOG(info, "Nb of parameters in the Z vs Sqrt(e) Fit: " << fNumParamsZFit);

    for (Int_t type = 0; type < fNumTypes; type++)
    {
        LOG(info) << "Types of Musli Hit Data : " << type + 1;
        LOG(info) << "fEaveVsBetaHitParams : " << fEaveVsBetaHitParams->GetAt(type);
    }
    for (Int_t type = 0; type < fNumTypes; type++)
    {
        LOG(info) << "Types of Musli Hit Data : " << type + 1;
        for (Int_t i = 0; i < fNumParamsEcorrBetaFit; i++)
        {
            LOG(info) << "fEcorrBetaHitParams : FitParam(" << i
                      << ") = " << fEcorrBetaHitParams->GetAt(type * fNumParamsEcorrBetaFit + i);
        }
    }
    for (Int_t type = 0; type < fNumTypes; type++)
    {
        LOG(info) << "Types of Musli Hit Data : " << type + 1;
        for (Int_t i = 0; i < fNumParamsZFit; i++)
        {
            LOG(info) << "fZHitParams : FitParam(" << i << ") = " << fZHitParams->GetAt(type * fNumParamsZFit + i);
        }
    }
}

ClassImp(R3BMusliHitPar);
