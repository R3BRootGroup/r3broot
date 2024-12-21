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

#include "R3BPtofHitPar.h"
#include <iostream>

R3BPtofHitPar::R3BPtofHitPar(const char* name, const char* title, const char* context)
    : FairParGenericSet(name, title, context)
    , fTOffset1(TArrayD(24))
    , fTOffset2(TArrayD(24))
    , fZScale(TArrayD(24))
    , fVEff(TArrayD(24))
{
}

void R3BPtofHitPar::putParams(FairParamList* list)
{
    std::cout << "-I- R3BPtofHitPar::putParams() called" << std::endl;
    if (!list)
        return;

    list->add("TOffset1Params", fTOffset1);
    list->add("TOffset2Params", fTOffset2);
    list->add("ZScaleParams", fZScale);
    list->add("VEffParams", fVEff);
}

Bool_t R3BPtofHitPar::getParams(FairParamList* list)
{
    std::cout << "-I- R3BPtofHitPar::getParams() called" << std::endl;
    if (!list)
        return kFALSE;
    std::cout << "-I- R3BPtofHitPar::getParams() 1 " << std::endl;

    if (!list->fill("TOffset1Params", &fTOffset1))
        return kFALSE;
    if (!list->fill("TOffset2Params", &fTOffset2))
        return kFALSE;
    if (!list->fill("ZScaleParams", &fZScale))
        return kFALSE;
    if (!list->fill("VEffParams", &fVEff))
        return kFALSE;

    return kTRUE;
}

void R3BPtofHitPar::Print(Option_t* option) const
{
    std::cout << "-I- Ptof Hit Parameters:" << std::endl;

    for (Int_t i = 0; i < 12; i++)
    {
        std::cout << "Paddle Nr. " << i << " : TOffset1 = " << fTOffset1.At(i) << "   TOffset2 = " << fTOffset2.At(i)
                  << "   ZScale = " << fZScale.At(i) << "   VEff = " << fVEff.At(i) << std::endl;
    }
}

ClassImp(R3BPtofHitPar)
