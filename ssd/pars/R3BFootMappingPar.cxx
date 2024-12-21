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

// ----------------------------------------------------------------
// -----        R3BFootMappingPar source file                 -----
// -----    Created 05/11/21  by J.L. Rodriguez-Sanchez       -----
// -----    for the info from the lookup table	        -----
// ----------------------------------------------------------------

#include "R3BFootMappingPar.h"

#include "FairLogger.h"
#include "FairParamList.h"

#include "TMath.h"
#include "TString.h"

#include <iostream>

// ---- Standard Constructor ---------------------------------------------------
R3BFootMappingPar::R3BFootMappingPar(const TString& name, const TString& title, const TString& context)
    : FairParGenericSet(name, title, context)
    , fNumDet(20)
    , fPolPar(5)
    , fGeo(202111)
{
    fDistance2target = new TArrayF(fNumDet);
    fAngleTheta = new TArrayF(fNumDet);
    fAnglePhi = new TArrayF(fNumDet);
    fOffsetX = new TArrayF(fNumDet);
    fOffsetY = new TArrayF(fNumDet);
    fEnevsPosCorr = new TArrayF(fPolPar * fNumDet);
}

// ----  Destructor ------------------------------------------------------------
R3BFootMappingPar::~R3BFootMappingPar()
{
    clear();
    if (fDistance2target)
        delete fDistance2target;
    if (fAngleTheta)
        delete fAngleTheta;
    if (fAnglePhi)
        delete fAnglePhi;
    if (fOffsetX)
        delete fOffsetX;
    if (fOffsetY)
        delete fOffsetY;
    if (fEnevsPosCorr)
        delete fEnevsPosCorr;
}

// ----  Method clear ----------------------------------------------------------
void R3BFootMappingPar::clear()
{
    status = kFALSE;
    resetInputVersions();
}

// ----  Method putParams ------------------------------------------------------
void R3BFootMappingPar::putParams(FairParamList* list)
{
    LOG(info) << "R3BFootMappingPar::putParams() called";
    if (!list)
    {
        return;
    }
    list->add("footGeoPar", fGeo);
    list->add("footDetPar", fNumDet);
    fDistance2target->Set(fNumDet);
    list->add("footDistance2targetPar", *fDistance2target);
    fAngleTheta->Set(fNumDet);
    list->add("footAngleThetaPar", *fAngleTheta);
    fAnglePhi->Set(fNumDet);
    list->add("footAnglePhiPar", *fAnglePhi);
    fOffsetX->Set(fNumDet);
    list->add("footOffsetXPar", *fOffsetX);
    fOffsetY->Set(fNumDet);
    list->add("footOffsetYPar", *fOffsetY);
    list->add("footPolPar", fPolPar);
    fEnevsPosCorr->Set(fNumDet * fPolPar);
    list->add("footEnevsPosCorrPar", *fEnevsPosCorr);
}

// ----  Method getParams ------------------------------------------------------
Bool_t R3BFootMappingPar::getParams(FairParamList* list)
{
    LOG(info) << "R3BFootMappingPar::getParams() called";
    if (!list)
    {
        return kFALSE;
    }
    if (!list->fill("footGeoPar", &fGeo))
    {
        LOG(info) << "Could not initialize footGeoPar";
        return kFALSE;
    }
    if (!list->fill("footDetPar", &fNumDet))
    {
        LOG(info) << "Could not initialize footDetPar";
        return kFALSE;
    }
    fDistance2target->Set(fNumDet);
    if (!(list->fill("footDistance2targetPar", fDistance2target)))
    {
        LOG(info) << "Could not initialize footDistance2targetPar";
        return kFALSE;
    }
    fAngleTheta->Set(fNumDet);
    if (!(list->fill("footAngleThetaPar", fAngleTheta)))
    {
        LOG(info) << "Could not initialize footAngleThetaPar";
        return kFALSE;
    }
    fAnglePhi->Set(fNumDet);
    if (!(list->fill("footAnglePhiPar", fAnglePhi)))
    {
        LOG(info) << "Could not initialize footAnglePhiPar";
        return kFALSE;
    }
    fOffsetX->Set(fNumDet);
    if (!(list->fill("footOffsetXPar", fOffsetX)))
    {
        LOG(info) << "Could not initialize footOffsetXPar";
        return kFALSE;
    }
    fOffsetY->Set(fNumDet);
    if (!(list->fill("footOffsetYPar", fOffsetY)))
    {
        LOG(info) << "Could not initialize footOffsetYPar";
        return kFALSE;
    }
    if (!list->fill("footPolPar", &fPolPar))
    {
        LOG(info) << "Could not initialize footPolPar";
        return kFALSE;
    }
    fEnevsPosCorr->Set(fNumDet * fPolPar);
    if (!(list->fill("footEnevsPosCorrPar", fEnevsPosCorr)))
    {
        LOG(info) << "Could not initialize footEnevsPosCorrPar";
        return kFALSE;
    }

    return kTRUE;
}

// ----  Method print ----------------------------------------------------------
void R3BFootMappingPar::print() { printParams(); }

// ----  Method printParams ----------------------------------------------------
void R3BFootMappingPar::printParams()
{
    LOG(info) << "R3BFootMappingPar::FOOT mapping parameters for geometry " << fGeo << ": ";

    LOG(info) << "Detector"
              << " "
              << "Distance_to_target"
              << " "
              << "Angle_theta"
              << " "
              << "Angle_Phi"
              << " "
              << "Offset_X"
              << " "
              << "Offset_Y";

    for (Int_t i = 0; i < fNumDet; i++)
    {
        LOG(info) << i + 1 << "\t"
                  << "\t" << fDistance2target->GetAt(i) << "\t\t" << fAngleTheta->GetAt(i) << "\t"
                  << fAnglePhi->GetAt(i) << "\t" << fOffsetX->GetAt(i) << "\t" << fOffsetY->GetAt(i);
    }

    LOG(info) << "Energy vs Position Correction Parameters : ";

    for (Int_t i = 0; i < fNumDet; i++)
    {
        LOG(info) << "FOOT " << i + 1 << " : ";

        for (Int_t j = 0; j < fPolPar; j++)
        {
            LOG(info) << "Pol" << j << "\t"
                      << "\t" << fEnevsPosCorr->GetAt(i * fPolPar + j);
        }
    }
}

ClassImp(R3BFootMappingPar);
