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

#include "R3BTarget2pDigiPar.h"
#include "FairLogger.h"
#include "TMath.h"

R3BTarget2pDigiPar::R3BTarget2pDigiPar(const char* name, const char* title, const char* context)
    : FairParGenericSet(name, title, context)
{
    clear();
}

void R3BTarget2pDigiPar::putParams(FairParamList* list)
{
    LOG(info) << "R3BTarget2pDigiPar::putParams() called";
    if (!list)
        return;
    //  list->add("max_paddle", (Int_t)nMaxPaddle);
    //  list->add("max_plane", (Int_t)nMaxPlane);
}

Bool_t R3BTarget2pDigiPar::getParams(FairParamList* list)
{
    LOG(info) << "R3BTarget2pDigiPar::getParams() called";
    if (!list)
        return kFALSE;
    LOG(info) << "R3BTarget2pDigiPar::getParams() 1 ";

    //  if (!list->fill("max_paddle", &nMaxPaddle)) return kFALSE;
    //  if (!list->fill("max_plane", &nMaxPlane)) return kFALSE;
    return kTRUE;
}

void R3BTarget2pDigiPar::Print(Option_t* option) const
{
    LOG(info) << "Target2p Digi Parameters:";
    //    LOG(info)<<"   Max Paddle   = "<<nMaxPaddle;
    //    LOG(info)<<"   Max Plane   = "<<nMaxPlane;
}

ClassImp(R3BTarget2pDigiPar);
