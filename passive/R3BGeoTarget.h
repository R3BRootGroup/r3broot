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

#ifndef R3BGEOTARGET_H
#define R3BGEOTARGET_H

#include "FairGeoSet.h"

class R3BGeoTarget : public FairGeoSet
{
  protected:
    char modName[2]; // name of module
    char eleName[2]; // substring for elements in module
  public:
    R3BGeoTarget();
    ~R3BGeoTarget() {}
    const char* getModuleName(Int_t) { return modName; }
    const char* getEleName(Int_t) { return eleName; }
    ClassDef(R3BGeoTarget, 0) // Class for geometry of Target
};

#endif /* !R3BGEOTARGET_H */
