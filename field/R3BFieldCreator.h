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

// -------------------------------------------------------------------------
// -----                    R3BFieldCreator header file                -----
// -----                Created 15/01/07 by M. Al-Turany               -----
// -------------------------------------------------------------------------

#pragma once

#include "R3BFieldPar.h"
#include <FairField.h>
#include <FairFieldFactory.h>

class R3BFieldCreator : public FairFieldFactory
{
  public:
    R3BFieldCreator();
    virtual ~R3BFieldCreator();
    virtual FairField* createFairField();
    virtual void SetParm();
    ClassDef(R3BFieldCreator, 1);

  protected:
    R3BFieldPar* fFieldPar;
};
