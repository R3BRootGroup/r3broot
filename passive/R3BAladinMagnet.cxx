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

#include "R3BAladinMagnet.h"
#include "FairLogger.h"

const Double_t __ALADIN_POS_DX = 0.;
const Double_t __ALADIN_POS_DY = 0.;
const Double_t __ALADIN_POS_DZ = 350. - 95.;
const Double_t __ALADIN_ROT = +7.;
const TString __ALADIN_NAME = "Aladin Magnet";

R3BAladinMagnet::R3BAladinMagnet()
    // NOTE: Delegate Constructor, such that geometry operations are set up.
    : R3BAladinMagnet("")
{
}

R3BAladinMagnet::R3BAladinMagnet(const TString& geoFile)
    : R3BModule(__ALADIN_NAME,
                __ALADIN_NAME,
                kFALSE,
                geoFile,
                { __ALADIN_POS_DX, __ALADIN_POS_DY, __ALADIN_POS_DZ },
                { "", -90., __ALADIN_ROT, 90. })
{
}

void R3BAladinMagnet::SetPosition(const TGeoTranslation&)
{
    LOG(fatal) << "Moving " << GetName() << " (which is a " << ClassName() << ") is not allowed!";
}

void R3BAladinMagnet::SetRotation(const TGeoRotation&)
{
    LOG(fatal) << "Rotating " << GetName() << " (which is a " << ClassName() << ") is not allowed!";
}

ClassImp(R3BAladinMagnet)
