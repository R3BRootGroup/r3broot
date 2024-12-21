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

#include "R3BPspxMappedPar.h"

#include "FairDetParIo.h"
#include "FairLogger.h"
#include "FairParIo.h"
#include "FairParamList.h"
#include "TMath.h"
#include "TString.h"

R3BPspxMappedPar::R3BPspxMappedPar(const char* name, const char* title, const char* context)
    : FairParGenericSet(name, title, context)
    , pspxmappedparstrip()
    , pspxmappedpardetector(-1)
{
    detName = "Pspx";
}

R3BPspxMappedPar::~R3BPspxMappedPar() { clear(); }

void R3BPspxMappedPar::clear()
{
    status = kFALSE;
    resetInputVersions();
}

void R3BPspxMappedPar::printparams()
{

    LOG(info) << "Print";
    LOG(info) << "pspxmappedpardetector: " << pspxmappedpardetector;
    Int_t size = pspxmappedparstrip.GetSize();
    LOG(info) << "size: " << size;
    for (Int_t i = 0; i < size; i++)
    {
        LOG(info) << i << " :" << pspxmappedparstrip.GetAt(i);
    }
}

void R3BPspxMappedPar::putParams(FairParamList* l)
{

    LOG(info) << "I am in R3BPspxMappedPar::putParams ";
    if (!l)
    {
        return;
    }
    l->add("R3BPspxMappedDetectors", pspxmappedpardetector);
    l->add("R3BPspxMappedStripsPerDetector", pspxmappedparstrip);
}

Bool_t R3BPspxMappedPar::getParams(FairParamList* l)
{

    LOG(info) << "I am in R3BPspxMappedPar::getParams ";

    if (!l)
    {
        return kFALSE;
    }
    if (!l->fill("R3BPspxMappedDetectors", &pspxmappedpardetector))
    {
        return kFALSE;
    }
    pspxmappedparstrip.Set(pspxmappedpardetector);

    if (!l->fill("R3BPspxMappedStripsPerDetector", &pspxmappedparstrip))
    {
        return kFALSE;
    }

    return kTRUE;
}

ClassImp(R3BPspxMappedPar)
