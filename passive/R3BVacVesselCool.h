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

#ifndef R3BVACVESSELCOOL_H
#define R3BVACVESSELCOOL_H

#include "R3BModule.h"

class R3BVacVesselCool : public R3BModule
{
  public:
    R3BVacVesselCool();

    R3BVacVesselCool(const TString& name);

    R3BVacVesselCool(const TString& name,
                     const TString& geoFile,
                     const TGeoTranslation& trans,
                     const TGeoRotation& rot = TGeoRotation());

    R3BVacVesselCool(const TString& name, const TString& geoFile, const TGeoCombiTrans& combi = TGeoCombiTrans());

    void ConstructGeometry(); // override;

    ClassDef(R3BVacVesselCool, 3)
    // ClassDefOverride(R3BVacVesselCool, 3)
};

#endif // R3BVACVESSELCOOL_H
