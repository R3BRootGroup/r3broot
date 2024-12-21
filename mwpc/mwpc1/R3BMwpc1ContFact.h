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

#ifndef R3BMwpc1ContFact_H
#define R3BMwpc1ContFact_H 1

#include "FairContFact.h"

class FairContainer;

class R3BMwpc1ContFact : public FairContFact
{
  private:
    void setAllContainers();

  public:
    R3BMwpc1ContFact();
    ~R3BMwpc1ContFact() {}
    FairParSet* createContainer(FairContainer*);
    ClassDef(R3BMwpc1ContFact, 0) // Factory for all MWPC1 parameter containers
};

#endif /* !R3BMwpc1ContFact_H */
