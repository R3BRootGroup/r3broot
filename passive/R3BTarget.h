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

#ifndef R3BTARGET_H
#define R3BTARGET_H 1

#include "R3BModule.h"

class R3BTGeoPar;

class R3BTarget : public R3BModule
{
  public:
    R3BTarget();

    R3BTarget(const TString& name);

    R3BTarget(const TString& name,
              const TString& geoFile,
              const TGeoTranslation& trans,
              const TGeoRotation& rot = TGeoRotation());

    R3BTarget(const TString& name, const TString& geoFile, const TGeoCombiTrans& combi = TGeoCombiTrans());

    virtual void ConstructGeometry();

    ClassDef(R3BTarget, 3) void SetParContainers();

  protected:
    R3BTGeoPar* fTGeoPar; //!
};

#endif // R3BTARGET_H
