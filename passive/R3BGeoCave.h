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

#ifndef R3BGEOCAVE_H
#define R3BGEOCAVE_H 1

#include "FairGeoSet.h"

#include <Rtypes.h>
#include <TString.h>
#include <fstream>
#include <iosfwd>

class FairGeoMedia;

class R3BGeoCave : public FairGeoSet
{
  protected:
    TString name;

  public:
    R3BGeoCave();
    ~R3BGeoCave() {}
    const char* getModuleName(Int_t) { return name.Data(); }
    Bool_t read(std::fstream&, FairGeoMedia*);
    void addRefNodes();
    void write(std::fstream&);
    void print();
    ClassDef(R3BGeoCave, 0) // Class for the geometry of CAVE
};

#endif /* !R3BGEOCAVE_H */
