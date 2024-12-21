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

#include "R3BVacVesselCool.h"
#include "FairLogger.h"

R3BVacVesselCool::R3BVacVesselCool()
    : R3BModule()
{
}

R3BVacVesselCool::R3BVacVesselCool(const TString& name)
    : R3BModule(name, "VacVesselCool", kFALSE)
{
}

R3BVacVesselCool::R3BVacVesselCool(const TString& name,
                                   const TString& geoFile,
                                   const TGeoTranslation& trans,
                                   const TGeoRotation& rot)
    : R3BVacVesselCool(name, geoFile, { trans, rot })
{
}

R3BVacVesselCool::R3BVacVesselCool(const TString& name, const TString& geoFile, const TGeoCombiTrans& combi)
    : R3BModule(name, "VacVesselCool", kFALSE, geoFile, combi)
{
}

void R3BVacVesselCool::ConstructGeometry()
{

    if (TString(GetName()) != TString("LiH"))
    {
        LOG(fatal) << "VacVesselCool is only supposed to work with LiH Target";
        // At least if the previous code was interpreted correctly
    }
    R3BModule::ConstructGeometry();
}

ClassImp(R3BVacVesselCool)
