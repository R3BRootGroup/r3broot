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

#ifndef R3BGFICONTFACT_H
#define R3BGFICONTFACT_H

#include "FairContFact.h"

class FairContainer;

class R3BGfiContFact : public FairContFact
{
  private:
    void setAllContainers();

  public:
    R3BGfiContFact();
    ~R3BGfiContFact() {}
    FairParSet* createContainer(FairContainer*);
    void activateParIo(FairParIo* io);
    ClassDef(R3BGfiContFact, 0) // Factory for all GFI parameter containers
};

#endif /* !R3BGFICONTFACT_H */
