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
// -----        R3BAmsMappingPar source file                  -----
// -----    Created 20/01/20  by J.L. Rodriguez-Sanchez       -----
// -----    for the info from the lookup table	        -----
// ----------------------------------------------------------------

#include "R3BAmsMappingPar.h"

#include "FairLogger.h"
#include "FairParamList.h"

#include "TMath.h"
#include "TString.h"

#include <iostream>

// ---- Standard Constructor ---------------------------------------------------
R3BAmsMappingPar::R3BAmsMappingPar(const TString& name, const TString& title, const TString& context)
    : FairParGenericSet(name, title, context)
    , fNumDet(6)
    , fGeo(2020)
{
    fSam = new TArrayI(fNumDet);
    fGtb = new TArrayI(fNumDet);
    fSiderem = new TArrayI(fNumDet);
    fDistance2target = new TArrayF(fNumDet);
    fAngleTheta = new TArrayF(fNumDet);
    fOffsetY = new TArrayF(fNumDet);
}

// ----  Destructor ------------------------------------------------------------
R3BAmsMappingPar::~R3BAmsMappingPar()
{
    clear();
    if (fSam)
        delete fSam;
    if (fGtb)
        delete fGtb;
    if (fSiderem)
        delete fSiderem;
    if (fDistance2target)
        delete fDistance2target;
    if (fAngleTheta)
        delete fAngleTheta;
    if (fOffsetY)
        delete fOffsetY;
}

// ----  Method clear ----------------------------------------------------------
void R3BAmsMappingPar::clear()
{
    status = kFALSE;
    resetInputVersions();
}

// ----  Method putParams ------------------------------------------------------
void R3BAmsMappingPar::putParams(FairParamList* list)
{
    LOG(info) << "R3BAmsMappingPar::putParams() called";
    if (!list)
    {
        return;
    }
    list->add("amsGeoPar", fGeo);
    list->add("amsDetPar", fNumDet);
    fSam->Set(fNumDet);
    list->add("amsSamPar", *fSam);
    fGtb->Set(fNumDet);
    list->add("amsGtbPar", *fGtb);
    fSiderem->Set(fNumDet);
    list->add("amsSideremPar", *fSiderem);
    fDistance2target->Set(fNumDet);
    list->add("amsDistance2targetPar", *fDistance2target);
    fAngleTheta->Set(fNumDet);
    list->add("amsAngleThetaPar", *fAngleTheta);
    fOffsetY->Set(fNumDet);
    list->add("amsOffsetYPar", *fOffsetY);
}

// ----  Method getParams ------------------------------------------------------
Bool_t R3BAmsMappingPar::getParams(FairParamList* list)
{
    LOG(info) << "R3BAmsMappingPar::getParams() called";
    if (!list)
    {
        return kFALSE;
    }
    if (!list->fill("amsGeoPar", &fGeo))
    {
        return kFALSE;
    }
    if (!list->fill("amsDetPar", &fNumDet))
    {
        return kFALSE;
    }
    fSam->Set(fNumDet);
    if (!(list->fill("amsSamPar", fSam)))
    {
        LOG(info) << "---Could not initialize amsSamPar";
        return kFALSE;
    }
    fGtb->Set(fNumDet);
    if (!(list->fill("amsGtbPar", fGtb)))
    {
        LOG(info) << "---Could not initialize amsGtbPar";
        return kFALSE;
    }
    fSiderem->Set(fNumDet);
    if (!(list->fill("amsSideremPar", fSiderem)))
    {
        LOG(info) << "---Could not initialize amsSideremPar";
        return kFALSE;
    }
    fDistance2target->Set(fNumDet);
    if (!(list->fill("amsDistance2targetPar", fDistance2target)))
    {
        LOG(info) << "---Could not initialize amsDistance2targetPar";
        return kFALSE;
    }
    fAngleTheta->Set(fNumDet);
    if (!(list->fill("amsAngleThetaPar", fAngleTheta)))
    {
        LOG(info) << "---Could not initialize amsAngleThetaPar";
        return kFALSE;
    }
    fOffsetY->Set(fNumDet);
    if (!(list->fill("amsOffsetYPar", fOffsetY)))
    {
        LOG(info) << "---Could not initialize amsOffsetYPar";
        return kFALSE;
    }

    return kTRUE;
}

// ----  Method print ----------------------------------------------------------
void R3BAmsMappingPar::print() { printParams(); }

// ----  Method printParams ----------------------------------------------------
void R3BAmsMappingPar::printParams()
{
    LOG(info) << "R3BAmsMappingPar: AMS Mapping Parameters for Geometry " << fGeo << ": ";

    LOG(info) << "Detector"
              << " "
              << "Sam"
              << " "
              << "    Gtb"
              << " "
              << "Siderem"
              << " "
              << "Distance_to_target"
              << " "
              << "Angle_theta"
              << " "
              << "Offset_Y";

    for (Int_t i = 0; i < fNumDet; i++)
    {
        LOG(info) << i + 1 << "\t" << fSam->GetAt(i) << "\t" << fGtb->GetAt(i) << "\t" << fSiderem->GetAt(i) << "\t"
                  << fDistance2target->GetAt(i) << "\t\t" << fAngleTheta->GetAt(i) << "\t" << fOffsetY->GetAt(i);
    }
}

ClassImp(R3BAmsMappingPar)
