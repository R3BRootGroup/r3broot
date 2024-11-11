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

#ifndef R3BPDCCONTFACT_H
#define R3BPDCCONTFACT_H

#include "FairContFact.h"

class FairContainer;

class R3BPdcContFact : public FairContFact
{
  private:
    void setAllContainers();

  public:
    R3BPdcContFact();
    ~R3BPdcContFact() {}
    FairParSet* createContainer(FairContainer*);
    void activateParIo(FairParIo* io);
    ClassDef(R3BPdcContFact, 1) // Factory for all Pdc parameter containers
};

#endif /* !R3BPDCCONTFACT_H */
