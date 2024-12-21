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
// -----------------------------------------------------------------
// -----           R3BPspxMappedPar header file                -----
// -----           Created 16/06/22  by I.Syndikus             -----
// -----           Modified Dec 2019 by M. Holl                -----
// -----------------------------------------------------------------

#include "R3BPspxHitPar.h"

#include "FairDetParIo.h"
#include "FairLogger.h"
#include "FairParIo.h"
#include "FairParamList.h"
#include "TMath.h"
#include "TString.h"

R3BPspxHitPar::R3BPspxHitPar(const char* name, const char* title, const char* context)
    : FairParGenericSet(name, title, context)
    , fNumDetectors(-1)
    , fHitPosPar()
    , fHitEPar()
{
    detName = "Pspx";
}

R3BPspxHitPar::~R3BPspxHitPar() { clear(); }

void R3BPspxHitPar::clear()
{
    status = kFALSE;
    resetInputVersions();
}

void R3BPspxHitPar::printParams()
{
    LOG(info) << "R3BPspxHitPar::printParams";
    LOG(info) << "fNumDetectors: " << fNumDetectors;
    Int_t size = fHitPosPar.GetSize();
    LOG(info) << "fHitPosPar size: " << size;
    for (Int_t i = 0; i < size; i++)
    {
        size = fHitEPar.GetSize();
        LOG(info) << i << " :" << fHitPosPar.GetAt(i);
    }
    LOG(info) << "fHitEPar size: " << size;
    for (Int_t i = 0; i < size; i++)
    {
        LOG(info) << i << " :" << fHitEPar.GetAt(i);
    }
}

void R3BPspxHitPar::putParams(FairParamList* l)
{
    LOG(info) << "I am in R3BPspxHitPar::putParams ";
    if (!l)
        return;
    l->add("R3BPspxHitDetectors", fNumDetectors);

    // count entries for position calibration,
    // can change when more than a linear polynomial function is used => more parameters
    Int_t array_size = (fNumDetectors * 6 + 1); // 6 entries per detector + 1 header word
    LOG(info) << "R3BPspxHitPosPar Array Size: " << array_size;
    fHitPosPar.Set(array_size);
    l->add("fHitPosPar", fHitPosPar);

    // count entries for energy calibration,
    array_size = (fNumDetectors * 5 + 1); // 5 entries per detector + 1 header word
    LOG(info) << "R3BPspxHitPosPar Array Size: " << array_size;
    fHitEPar.Set(array_size);
    l->add("fHitEPar", fHitEPar);
}

Bool_t R3BPspxHitPar::getParams(FairParamList* l)
{
    LOG(info) << "I am in R3BPspxCalPar::getParams ";

    if (!l)
        return kFALSE;
    if (!l->fill("R3BPspxHitDetectors", &fNumDetectors))
        return kFALSE;

    // count entries for position calibration,
    // can change when more than a linear polynomial function is used => more parameters
    Int_t array_size = (fNumDetectors * 7 + 1); // 6 entries per detector + 1 header word
    LOG(info) << "R3BPspxHitPosPar Array Size: " << array_size;
    fHitPosPar.Set(array_size);
    if (!(l->fill("R3BPspxHitPosPar", &fHitPosPar)))
    {
        LOG(warn) << "Could not initialize R3BPspxHitPosPar";
        return kFALSE;
    }

    // count entries for energy calibration,
    array_size = (fNumDetectors * 5 + 1); // 4 entries per detector + 1 header word
    LOG(info) << "R3BPspxHitEPar Array Size: " << array_size;
    fHitEPar.Set(array_size);
    if (!(l->fill("R3BPspxHitEPar", &fHitEPar)))
    {
        LOG(warn) << "Could not initialize R3BPspxHitEPar";
        return kFALSE;
    }

    return kTRUE;
}

ClassImp(R3BPspxHitPar)
