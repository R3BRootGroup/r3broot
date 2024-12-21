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
// -----                    R3BFieldContFact header file               -----
// -------------------------------------------------------------------------

/** R3BFieldContFact.h
 ** @author
 ** @since 20.02.2006
 ** @version 1.0
 **
 ** Parameter container factory for field parameter container
 **/

#pragma once

#include <FairContFact.h>

class FairContainer;
class FairParSet;

class R3BFieldContFact : public FairContFact
{
  public:
    /** Constructor **/
    R3BFieldContFact();

    /** Destructor **/
    ~R3BFieldContFact();

    /** Create containers
     ** Creates the requested parameter sets (R3BFieldPar)
     **/
    FairParSet* createContainer(FairContainer* container);

  private:
    /** Set all containers
     ** Creates container objects with all accepted contexts and adds them
     ** to the list of containers for the field library.
     **/
    void SetAllContainers();

  public:
    ClassDef(R3BFieldContFact, 1);
};
