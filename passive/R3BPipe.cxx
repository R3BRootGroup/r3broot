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

#include "R3BPipe.h"
#include "TGeoManager.h"

R3BPipe::R3BPipe()
    : R3BModule()
{
}

R3BPipe::R3BPipe(const TString& name)
    : R3BModule(name, "R3BPipe", kFALSE)
{
}

R3BPipe::R3BPipe(const TString& geoFile, const TGeoTranslation& trans, const TGeoRotation& rot)
    : R3BPipe(geoFile, { trans, rot })
{
}

R3BPipe::R3BPipe(const TString& geoFile, const TGeoCombiTrans& combi)
    : R3BModule("R3BPipe", "Pipe", kFALSE, geoFile, combi)
{
}

void R3BPipe::ConstructGeometry() { R3BModule::ConstructGeometry(); }

ClassImp(R3BPipe);
