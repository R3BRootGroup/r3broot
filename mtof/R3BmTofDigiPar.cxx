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

#include "R3BmTofDigiPar.h"
#include "TMath.h"
#include <iostream>

R3BmTofDigiPar::R3BmTofDigiPar(const char* name, const char* title, const char* context)
    : FairParGenericSet(name, title, context)
{
    clear();
}

void R3BmTofDigiPar::putParams(FairParamList* list)
{
    std::cout << "-I- R3BmTofDigiPar::putParams() called" << std::endl;
    if (!list)
        return;
    //  list->add("max_paddle", (Int_t)nMaxPaddle);
    //  list->add("max_plane", (Int_t)nMaxPlane);
}

Bool_t R3BmTofDigiPar::getParams(FairParamList* list)
{
    std::cout << "-I- R3BmTofDigiPar::getParams() called" << std::endl;
    if (!list)
        return kFALSE;
    std::cout << "-I- R3BmTofDigiPar::getParams() 1 ";

    //  if (!list->fill("max_paddle", &nMaxPaddle)) return kFALSE;
    //  if (!list->fill("max_plane", &nMaxPlane)) return kFALSE;
    return kTRUE;
}

void R3BmTofDigiPar::Print(Option_t* option) const
{
    std::cout << "-I- mTof Digi Parameters:" << std::endl;
    //    std::cout<<"   Max Paddle   = "<<nMaxPaddle<<std::endl;
    //    std::cout<<"   Max Plane   = "<<nMaxPlane<<std::endl;
}

ClassImp(R3BmTofDigiPar);
