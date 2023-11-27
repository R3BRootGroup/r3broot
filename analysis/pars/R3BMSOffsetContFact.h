/******************************************************************************
 *   Copyright (C) 2019 GSI Helmholtzzentrum für Schwerionenforschung GmbH    *
 *   Copyright (C) 2019-2024 Members of R3B Collaboration                     *
 *                                                                            *
 *             This software is distributed under the terms of the            *
 *                 GNU General Public Licence (GPL) version 3,                *
 *                    copied verbatim in the file "LICENSE".                  *
 *                                                                            *
 * In applying this license GSI does not waive the privileges and immunities  *
 * granted to it by virtue of its status as an Intergovernmental Organization *
 * or submit itself to any jurisdiction.                                      *
 ******************************************************************************/
#pragma once

#include "FairContFact.h"

class FairContainer;

class R3BMSOffsetContFact : public FairContFact
{
  private:
    void setAllContainers();

  public:
    R3BMSOffsetContFact();
    R3BMSOffsetContFact(const R3BMSOffsetContFact&) = delete;
    R3BMSOffsetContFact(R3BMSOffsetContFact&&) = delete;
    R3BMSOffsetContFact& operator=(const R3BMSOffsetContFact&) = delete;
    R3BMSOffsetContFact& operator=(R3BMSOffsetContFact&&) = delete;
    ~R3BMSOffsetContFact() override = default;
    FairParSet* createContainer(FairContainer* /*unused*/) override;
    ClassDefOverride(R3BMSOffsetContFact, 0) // Factory for all MSOffset parameter containers
};
