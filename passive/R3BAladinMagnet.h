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

#ifndef R3BALADINMAGNET_H
#define R3BALADINMAGNET_H

#include "R3BModule.h"

class R3BAladinMagnet : public R3BModule
{
  public:
    // NOTE:
    // - There is no scenario where name and title should be set to something different than "Aladin Magnet"
    // - Aladin should not be moved or rotated from its given position
    // ... thus don't offer the option to.
    // TODO: Consider if the default constructor should be given at all, as the only valid operation is setting the
    // geometry file name
    R3BAladinMagnet();
    R3BAladinMagnet(const TString& geoFile);

    void SetPosition(const TGeoTranslation&); // override;
    void SetRotation(const TGeoRotation&);    // override;

  private:
    ClassDef(R3BAladinMagnet, 3)
    // ClassDefOverride(R3BAladinMagnet, 3)
};

#endif // R3BALADINMAGNET_H
