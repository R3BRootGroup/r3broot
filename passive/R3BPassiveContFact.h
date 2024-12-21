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

#ifndef R3BPASSIVECONTFACT_H
#define R3BPASSIVECONTFACT_H

#include "FairContFact.h"
#include "FairLogger.h"
#include "FairRuntimeDb.h"

class FairContainer;

class R3BPassiveContFact : public FairContFact
{
  private:
    void setAllContainers();

  public:
    R3BPassiveContFact();
    ~R3BPassiveContFact() {}
    FairParSet* createContainer(FairContainer*);
    ClassDef(R3BPassiveContFact, 0) // Factory for all Passive parameter containers
};

#endif /* !R3BPASSIVECONTFACT_H */
