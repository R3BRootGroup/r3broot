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
// ------------------------------------------------------------
// -----           R3BPspxPrecalPar header file           -----
// -----           Created 17/03/13  by I.Syndikus        -----
// -----           Modified Dec 2019 by M. Holl           -----
// ------------------------------------------------------------

#include "R3BPspxPrecalPar.h"

#include "FairDetParIo.h"
#include "FairLogger.h"
#include "FairParIo.h"
#include "FairParamList.h"
#include "TMath.h"
#include "TString.h"

R3BPspxPrecalPar::R3BPspxPrecalPar(const char* name, const char* title, const char* context)
    : FairParGenericSet(name, title, context)
    , fNumDetectors(-1)
    , fNumStrips()
    , fPrecalPar()
{
    detName = "Pspx";
}

R3BPspxPrecalPar::~R3BPspxPrecalPar() { clear(); }

void R3BPspxPrecalPar::clear()
{
    status = kFALSE;
    resetInputVersions();
}

void R3BPspxPrecalPar::printParams()
{
    LOG(info) << "R3BPspxPrecalPar::printParams";
    LOG(info) << "PSPX: fNumDetectors: " << fNumDetectors;
    Int_t size = fNumStrips.GetSize();
    LOG(info) << "fNumStrips size: " << size;
    LOG(info) << "Detectorsno.: No. of Strips";
    for (Int_t i = 0; i < size; i++)
    {
        LOG(info) << i << " :" << fNumStrips.GetAt(i);
    }

    size = fPrecalPar.GetSize();
    LOG(info) << "fPrecalPar size: " << size;
    for (Int_t i = 0; i < size; i++)
    {
        LOG(info) << i << " :" << fPrecalPar.GetAt(i);
    }
}

void R3BPspxPrecalPar::putParams(FairParamList* l)
{
    LOG(info) << "I am in R3BPspxPrecalPar::putParams ";
    if (!l)
        return;

    l->add("R3BPspxPrecalDetectors", fNumDetectors);
    l->add("R3BPspxPrecalStrips", fNumStrips);

    Int_t count_strips = 0;
    for (Int_t i = 0; i < fNumDetectors; i++)
    {
        count_strips += fNumStrips[i];
    }

    // count all entries: lines with strip info (4 entries) + lines with detector info (3 entries)
    Int_t array_size = (count_strips * 4 + fNumDetectors * 3);
    LOG(info) << "R3BPspxPrecalPar Array Size: " << array_size;
    fPrecalPar.Set(array_size);
    l->add("R3BPspxPrecalPar", fPrecalPar);
}

Bool_t R3BPspxPrecalPar::getParams(FairParamList* l)
{
    // print();
    LOG(info) << "I am in R3BPspxPrecalPar::getParams ";

    if (!l)
    {
        return kFALSE;
    }

    fNumStrips.Set(fNumDetectors);

    if (!l->fill("R3BPspxPrecalDetectors", &fNumDetectors))
    {
        return kFALSE;
    }
    if (!l->fill("R3BPspxPrecalStrips", &fNumStrips))
    {
        return kFALSE;
    }

    Int_t count_strips = 0;
    for (Int_t i = 0; i < fNumDetectors; i++)
    {
        count_strips += fNumStrips[i];
    }

    LOG(info) << "Total number of strips: " << count_strips;
    // count all entries: lines with strip info (4 entries) + lines with detector info (3 entries)
    Int_t array_size = (count_strips * 4 + fNumDetectors * 3);
    LOG(info) << "R3BPspxPrecalPar Array Size: " << array_size;
    fPrecalPar.Set(array_size);
    if (!(l->fill("R3BPspxPrecalPar", &fPrecalPar)))
    {
        LOG(warn) << "Could not initialize R3BPspxPrecalPar";
        return kFALSE;
    }

    return kTRUE;
}

ClassImp(R3BPspxPrecalPar)
