/******************************************************************************
 *   Copyright (C) 2019 GSI Helmholtzzentrum für Schwerionenforschung GmbH    *
 *   Copyright (C) 2019 Members of R3B Collaboration                          *
 *                                                                            *
 *             This software is distributed under the terms of the            *
 *                 GNU General Public Licence (GPL) version 3,                *
 *                    copied verbatim in the file "LICENSE".                  *
 *                                                                            *
 * In applying this license GSI does not waive the privileges and immunities  *
 * granted to it by virtue of its status as an Intergovernmental Organization *
 * or submit itself to any jurisdiction.                                      *
 ******************************************************************************/

#ifndef R3BGEOLAND_H
#define R3BGEOLAND_H

#include "FairGeoSet.h"

class R3BGeoLand : public FairGeoSet
{
  protected:
    char modName[30]; // name of module
    char eleName[20]; // substring for elements in module
  public:
    R3BGeoLand();
    ~R3BGeoLand() {}
    const char* getModuleName(Int_t);
    const char* getEleName(Int_t);
    inline Int_t getModNumInMod(const TString&);
    ClassDef(R3BGeoLand, 0) // Class for STS
};

inline Int_t R3BGeoLand::getModNumInMod(const TString& name)
{
    // returns the module index from module name
    return (Int_t)(name[3] - '0') - 1;
}

#endif /* !R3BGEOLAND_H */
