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

// ------------------------------------------------------------------
// -----             R3BMusicContFact source file               -----
// -----         Created 24/11/19  by J.L. Rodriguez-Sanchez    -----
// ------------------------------------------------------------------

#ifndef R3BMusicContFact_H
#define R3BMusicContFact_H

#include "FairContFact.h"

class FairContainer;

class R3BMusicContFact : public FairContFact
{
  private:
    void setAllContainers();

  public:
    R3BMusicContFact();
    ~R3BMusicContFact() {}
    FairParSet* createContainer(FairContainer*);
    ClassDef(R3BMusicContFact, 0) // Factory for all R3BMUSIC parameter containers
};

#endif /* R3BMusicContFact_H */
