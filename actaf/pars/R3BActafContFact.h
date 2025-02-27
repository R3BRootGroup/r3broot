/******************************************************************************
 *   Copyright (C) 2025 GSI Helmholtzzentrum für Schwerionenforschung GmbH    *
 *   Copyright (C) 2025 Members of R3B Collaboration                          *
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

#include <FairContFact.h>

class FairContainer;

class R3BActafContFact : public FairContFact
{
  public:
    R3BActafContFact();
    virtual ~R3BActafContFact() = default;

    FairParSet* createContainer(FairContainer*) override;

  private:
    void setAllContainers();

  public:
    ClassDefOverride(R3BActafContFact, 0)
};
