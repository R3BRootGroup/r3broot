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

#ifndef R3BCOLLIMATOR_H
#define R3BCOLLIMATOR_H

#include "R3BModule.h"

class R3BTGeoPar;

class R3BCollimator : public R3BModule
{
  public:
    R3BCollimator();

    R3BCollimator(const TString& name);

    R3BCollimator(const TString& name,
                  const TString& geoFile,
                  const TGeoTranslation& trans,
                  const TGeoRotation& rot = TGeoRotation());

    R3BCollimator(const TString& name, const TString& geoFile, const TGeoCombiTrans& combi = TGeoCombiTrans());

    virtual void ConstructGeometry();

    ClassDef(R3BCollimator, 1) void SetParContainers();

  protected:
    R3BTGeoPar* fTGeoPar; //!
};

#endif // R3BCOLLIMATOR_H
