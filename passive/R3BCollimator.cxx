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

#include "R3BCollimator.h"

#include "TGeoManager.h"
#include "TGeoTube.h"

R3BCollimator::R3BCollimator()
    : R3BModule()
{
}

R3BCollimator::R3BCollimator(const TString& name)
    : R3BModule(name, "Collimator", kFALSE)
{
}

R3BCollimator::R3BCollimator(const TString& name,
                             const TString& geoFile,
                             const TGeoTranslation& trans,
                             const TGeoRotation& rot)
    : R3BCollimator(name, geoFile, { trans, rot })
{
}

R3BCollimator::R3BCollimator(const TString& name, const TString& geoFile, const TGeoCombiTrans& combi)
    : R3BModule(name, "Collimator", kFALSE, geoFile, combi)
{
}

void R3BCollimator::ConstructGeometry() { R3BModule::ConstructGeometry(); }

void R3BCollimator::SetParContainers() {}

ClassImp(R3BCollimator)
