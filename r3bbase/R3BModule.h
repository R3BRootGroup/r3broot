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

#ifndef R3BMODULE_H
#define R3BMODULE_H

#include "FairModule.h"
#include "TGeoMatrix.h"
#include "TString.h"

class R3BModule : public FairModule
{
  public:
    R3BModule();

    R3BModule(const TString& name, const TString& title, const Bool_t active = kFALSE);

    R3BModule(const TString& name,
              const TString& title,
              const Bool_t active,
              const TString& geoFile,
              const TGeoTranslation& trans,
              const TGeoRotation& rot = TGeoRotation());

    R3BModule(const TString& name,
              const TString& title,
              const Bool_t active,
              const TString& geoFile,
              const TGeoCombiTrans& combi = TGeoCombiTrans());

    /** Construct geometry from root files */
    virtual void ConstructGeometry() override;
    virtual void ConstructRootGeometry(TGeoMatrix* p = nullptr) override;

    /** Method to setup the position of the detector-origin. */
    virtual void SetPosition(const TGeoTranslation& trans) { fCombiTrans.SetTranslation(trans); }

    /** Method to setup the rotation of the detector */
    virtual void SetRotation(const TGeoRotation& rot) { fCombiTrans.SetRotation(rot); }

  protected:
    TGeoCombiTrans fCombiTrans;

    // ClassDef(R3BModule, 3)
    ClassDefOverride(R3BModule, 3)
};

#endif // R3BMODULE_H
