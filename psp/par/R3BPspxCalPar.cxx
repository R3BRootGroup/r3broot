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
// -----           Created 16/05/12  by I.Syndikus             -----
// -----           Modified Dec 2019 by M. Holl                -----
// -----------------------------------------------------------------

#include "R3BPspxCalPar.h"

#include "FairDetParIo.h"
#include "FairLogger.h"
#include "FairParIo.h"
#include "FairParamList.h"
#include "TMath.h"
#include "TString.h"

R3BPspxCalPar::R3BPspxCalPar(const char* name, const char* title, const char* context)
    : FairParGenericSet(name, title, context)
    , fNumDetectors(-1)
    , fNumStrips()
    , fCalPar()
{
    detName = "Pspx";
}

R3BPspxCalPar::~R3BPspxCalPar() { clear(); }

void R3BPspxCalPar::clear()
{
    status = kFALSE;
    resetInputVersions();
}

void R3BPspxCalPar::printParams()
{

    LOG(info) << "R3BPspxCalPar::printParams";
    LOG(info) << "fNumDetectors: " << fNumDetectors;
    Int_t size = fNumStrips.GetSize();
    LOG(info) << "fNumStrips size: " << size;
    LOG(info) << "Detectorno.: No. of Strips";
    for (Int_t i = 0; i < size; i++)
    {
        LOG(info) << i << " :" << fNumStrips.GetAt(i);
    }

    size = fCalPar.GetSize();
    LOG(info) << "fCalPar size: " << size;
    for (Int_t i = 0; i < size; i++)
    {
        LOG(info) << i << " :" << fCalPar.GetAt(i);
    }
}

void R3BPspxCalPar::putParams(FairParamList* l)
{
    LOG(info) << "I am in R3BPspxCalPar::putParams ";
    if (!l)
        return;
    l->add("R3BPspxCalDetectors", fNumDetectors);
    l->add("R3BPspxCalStrips", fNumStrips);

    Int_t count_strips = 0;
    for (Int_t i = 0; i < fNumDetectors; i++)
    {
        count_strips += fNumStrips[i];
    }

    // count all entries: lines with strip info (2 entries) + lines with detector info (3 entries)
    Int_t array_size = (count_strips * 2 + fNumDetectors * 3);
    LOG(info) << "R3BPspxCalPar Array Size: " << array_size;
    fCalPar.Set(array_size);
    l->add("R3BPspxCalPar", fCalPar);
}

Bool_t R3BPspxCalPar::getParams(FairParamList* l)
{
    LOG(info) << "I am in R3BPspxCalPar::getParams ";

    if (!l)
        return kFALSE;
    if (!l->fill("R3BPspxCalDetectors", &fNumDetectors))
        return kFALSE;

    fNumStrips.Set(fNumDetectors);

    if (!l->fill("R3BPspxCalStrips", &fNumStrips))
        return kFALSE;

    Int_t count_strips = 0;
    for (Int_t i = 0; i < fNumDetectors; i++)
    {
        count_strips += fNumStrips[i];
    }
    LOG(info) << "Total number of strips: " << count_strips;

    // count all entries: lines with strip info (2 entries) + lines with detector info (3 entries)
    Int_t array_size = (count_strips * 2 + fNumDetectors * 3);
    LOG(info) << "R3BPspxCalPar Array Size: " << array_size;
    fCalPar.Set(array_size);
    if (!(l->fill("R3BPspxCalPar", &fCalPar)))
    {
        LOG(warn) << "Could not initialize R3BPspxCalPar";
        return kFALSE;
    }

    return kTRUE;
}

ClassImp(R3BPspxCalPar)
