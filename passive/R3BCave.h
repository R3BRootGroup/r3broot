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

#ifndef R3BCAVE_H
#define R3BCAVE_H 1

#include "FairDetector.h"
#include "FairModule.h"

#include "Rtypes.h"

class R3BCave : public FairModule
{
  public:
    R3BCave(const char* name, const char* Title = "Exp Cave");
    R3BCave();
    virtual ~R3BCave();
    virtual void ConstructGeometry();

    virtual FairModule* CloneModule() const;

  private:
    Double_t world[3];

    R3BCave(const R3BCave&);
    R3BCave& operator=(const R3BCave&);

    ClassDef(R3BCave, 2);
};

#endif // R3BCave_H
