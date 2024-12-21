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

#include "R3BLosHitPar.h"
#include "FairLogger.h"
#include "FairParamList.h"
#include "TArrayF.h"
#include "TMath.h"
#include "TString.h"

#include <iostream>

// ---- Standard Constructor ---------------------------------------------------
R3BLosHitPar::R3BLosHitPar(const TString& name, const TString& title, const TString& context)
    : FairParGenericSet(name, title, context)
    , fNumParamsTamexLE(2)
    , fNumParamsTamexTE(2)
{
    fLEMatchParams = new TArrayF(8 * fNumParamsTamexLE);
    fTEMatchParams = new TArrayF(8 * fNumParamsTamexTE);
}

// ----  Destructor ------------------------------------------------------------
R3BLosHitPar::~R3BLosHitPar()
{
    clear();
    if (fLEMatchParams)
        delete fLEMatchParams;
    if (fTEMatchParams)
        delete fTEMatchParams;
}

// ----  Method clear ----------------------------------------------------------
void R3BLosHitPar::clear()
{
    status = kFALSE;
    resetInputVersions();
}

// ----  Method putParams ------------------------------------------------------
void R3BLosHitPar::putParams(FairParamList* list)
{
    LOG(info) << "R3BLosHitPar::putParams() called";
    if (!list)
    {
        return;
    }
    list->add("p0", fp0);
    list->add("p1", fp1);
    list->add("xoffset_MCFD", fxoffset_MCFD);
    list->add("yoffset_MCFD", fyoffset_MCFD);
    list->add("xveff_MCFD", fxveff_MCFD);
    list->add("yveff_MCFD", fyveff_MCFD);
    list->add("NumParamsTamexLE", fNumParamsTamexLE);
    list->add("NumParamsTamexTE", fNumParamsTamexTE);

    Int_t array_LE = 8 * fNumParamsTamexLE;
    LOG(info) << "Array size LOS Tamex LE matching: " << array_LE;
    fLEMatchParams->Set(array_LE);
    list->add("LosLEMatchPar", *fLEMatchParams);

    Int_t array_TE = 8 * fNumParamsTamexTE;
    LOG(info) << "Array size LOS Tamex TE matching: " << array_TE;
    fTEMatchParams->Set(array_TE);
    list->add("LosTEMatchPar", *fTEMatchParams);
}

// ----  Method getParams ------------------------------------------------------
Bool_t R3BLosHitPar::getParams(FairParamList* list)
{
    LOG(info) << "R3BLosHitPar::getParams() called";
    if (!list)
    {
        return kFALSE;
    }

    if (!(list->fill("p0", &fp0)))
    {
        LOG(info) << "---Could not initialize LosHit p0 Par";
        return kFALSE;
    }

    if (!(list->fill("p1", &fp1)))
    {
        LOG(info) << "---Could not initialize LosHit p1 Par";
        return kFALSE;
    }

    if (!(list->fill("xoffset_MCFD", &fxoffset_MCFD)))
    {
        LOG(info) << "---Could not initialize LosHit xoffset_MCFD Par";
        return kFALSE;
    }

    if (!(list->fill("yoffset_MCFD", &fyoffset_MCFD)))
    {
        LOG(info) << "---Could not initialize LosHit yoffset_MCFD Par";
        return kFALSE;
    }

    if (!(list->fill("xveff_MCFD", &fxveff_MCFD)))
    {
        LOG(info) << "---Could not initialize LosHit xveff_MCFD Par";
        return kFALSE;
    }

    if (!(list->fill("yveff_MCFD", &fyveff_MCFD)))
    {
        LOG(info) << "---Could not initialize LosHit yveff_MCFD Par";
        return kFALSE;
    }

    if (!list->fill("NumParamsTamexLE", &fNumParamsTamexLE))
    {
        LOG(info) << "Could not initialize NumParamsTamexLE";
        return kFALSE;
    }

    if (!list->fill("NumParamsTamexTE", &fNumParamsTamexTE))
    {
        LOG(info) << "Could not initialize NumParamsTamexTE";
        return kFALSE;
    }

    Int_t array_LE = 8 * fNumParamsTamexLE;
    fLEMatchParams->Set(array_LE);
    if (!(list->fill("LosLEMatchPar", fLEMatchParams)))
    {
        LOG(info) << "Could not initialize LosLEMatchPar";
        return kFALSE;
    }

    Int_t array_TE = 8 * fNumParamsTamexTE;
    fLEMatchParams->Set(array_TE);
    if (!(list->fill("LosTEMatchPar", fTEMatchParams)))
    {
        LOG(info) << "Could not initialize LosTEMatchPar";
        return kFALSE;
    }

    return kTRUE;
}

// ----  Method print ----------------------------------------------------------
void R3BLosHitPar::print() { printParams(); }

// ----  Method printParams ----------------------------------------------------
void R3BLosHitPar::printParams()
{
    LOG(info) << "R3BLosHitPar: Los calibration parameters (Z-charge):";

    LOG(info) << "p0"
              << " "
              << "p1"
              << " ";

    LOG(info) << fp0 << "\t" << fp1;

    LOG(info) << "R3BLosHitPar: LOS Tamex LE match par:";
    LOG(info) << "Nb of Tamex LE match Par: ";
    LOG(info) << fNumParamsTamexLE;
    for (Int_t i = 0; i < 8; i++)
    {
        for (Int_t j = 0; j < fNumParamsTamexLE; j++)
        {
            LOG(info) << "FitParam(" << j << ") = " << fLEMatchParams->GetAt(i * fNumParamsTamexLE + j);
        }
    }

    LOG(info) << "R3BLosHitPar: LOS Tamex TE match par:";
    LOG(info) << "Nb of Tamex TE match Par: ";
    LOG(info) << fNumParamsTamexTE;
    for (Int_t i = 0; i < 8; i++)
    {
        for (Int_t j = 0; j < fNumParamsTamexTE; j++)
        {
            LOG(info) << "FitParam(" << j << ") = " << fTEMatchParams->GetAt(i * fNumParamsTamexTE + j);
        }
    }
}

ClassImp(R3BLosHitPar)
