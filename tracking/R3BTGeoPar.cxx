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

#include "R3BTGeoPar.h"

#include "FairLogger.h"
#include "FairParamList.h" // for FairParamList
#include "FairRtdbRun.h"

ClassImp(R3BTGeoPar);

R3BTGeoPar::R3BTGeoPar(const char* name, const char* title, const char* context, Bool_t own)
    : FairParGenericSet(name, title, context, own)
    , fPosX(0.)
    , fPosY(0.)
    , fPosZ(0.)
    , fRotX(0.)
    , fRotY(0.)
    , fRotZ(0.)
    , fDimX(0.)
    , fDimY(0.)
    , fDimZ(0.)
    , fSigmaX(0.)
    , fSigmaY(0.)
    , fZ(0.)
    , fA(0.)
    , fDensity(0.)
    , fI(0.)
{
}

R3BTGeoPar::~R3BTGeoPar() { clear(); }

void R3BTGeoPar::putParams(FairParamList* list)
{
    LOG(info) << "R3BTGeoPar::putParams() called";
    if (!list)
    {
        return;
    }
    list->add("PositionX", fPosX);
    list->add("PositionY", fPosY);
    list->add("PositionZ", fPosZ);
    list->add("RotationX", fRotX);
    list->add("RotationY", fRotY);
    list->add("RotationZ", fRotZ);
    list->add("DimensionX", fDimX);
    list->add("DimensionY", fDimY);
    list->add("DimensionZ", fDimZ);
    list->add("SigmaX", fSigmaX);
    list->add("SigmaY", fSigmaY);
    list->add("Z", fZ);
    list->add("A", fA);
    list->add("Density", fDensity);
    list->add("IonisationEnergy", fI);
}

Bool_t R3BTGeoPar::getParams(FairParamList* list)
{
    if (!list)
    {
        return kFALSE;
    }
    if (!list->fill("PositionX", &fPosX))
    {
        return kFALSE;
    }
    if (!list->fill("PositionY", &fPosY))
    {
        return kFALSE;
    }
    if (!list->fill("PositionZ", &fPosZ))
    {
        return kFALSE;
    }
    if (!list->fill("RotationX", &fRotX))
    {
        return kFALSE;
    }
    if (!list->fill("RotationY", &fRotY))
    {
        return kFALSE;
    }
    if (!list->fill("RotationZ", &fRotZ))
    {
        return kFALSE;
    }
    if (!list->fill("DimensionX", &fDimX))
    {
        return kFALSE;
    }
    if (!list->fill("DimensionY", &fDimY))
    {
        return kFALSE;
    }
    if (!list->fill("DimensionZ", &fDimZ))
    {
        return kFALSE;
    }
    if (!list->fill("SigmaX", &fSigmaX))
    {
        return kFALSE;
    }
    if (!list->fill("SigmaY", &fSigmaY))
    {
        return kFALSE;
    }
    if (!list->fill("Z", &fZ))
    {
        return kFALSE;
    }
    if (!list->fill("A", &fA))
    {
        return kFALSE;
    }
    if (!list->fill("Density", &fDensity))
    {
        return kFALSE;
    }
    if (!list->fill("IonisationEnergy", &fI))
    {
        return kFALSE;
    }
    return kTRUE;
}

void R3BTGeoPar::clear()
{
    status = kFALSE;
    resetInputVersions();
}

void R3BTGeoPar::SetMaterial(const Double_t& Z, const Double_t& A, const Double_t& density, const Double_t& I)
{
    fZ = Z;
    fA = A;
    fDensity = density;
    fI = I;
}

void R3BTGeoPar::printParams()
{
    LOG(info) << " -----------  " << GetName() << " Geometry Parameters -------------  ";
    LOG(info) << " Position in cave: " << fPosX << " " << fPosY << " " << fPosZ;
    LOG(info) << " Rotation in deg: " << fRotX << " " << fRotY << " " << fRotZ;
    LOG(info) << " Dimensions: " << fDimX << " " << fDimY << " " << fDimZ;
    LOG(info) << " SigmaXY: " << fSigmaX << " " << fSigmaY;
    LOG(info) << " Z=" << fZ << " A=" << fA << " Density=" << fDensity << " Ionisation=" << fI;
    LOG(info) << " ------------------------------------------------------------------  ";
}
