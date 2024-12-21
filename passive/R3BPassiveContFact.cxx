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

//*-- AUTHOR : Denis Bertini
//*-- Created : 21/06/2005

/////////////////////////////////////////////////////////////
//
//  R3BPassiveContFact
//
//  Factory for the parameter containers in libPassive
//
/////////////////////////////////////////////////////////////

#include "R3BPassiveContFact.h"
#include "R3BGeoPassivePar.h"
#include "R3BTGeoPar.h"

static R3BPassiveContFact gR3BPassiveContFact;

R3BPassiveContFact::R3BPassiveContFact()
{
    // Constructor (called when the library is loaded)
    fName = "R3BPassiveContFact";
    fTitle = "Factory for parameter containers in libPassive";
    setAllContainers();
    FairRuntimeDb::instance()->addContFactory(this);
}

void R3BPassiveContFact::setAllContainers()
{
    /** Creates the Container objects with all accepted contexts and adds them to
     *  the list of containers for the STS library.*/

    FairContainer* p = new FairContainer("R3BGeoPassivePar", "Passive Geometry Parameters", "TestDefaultContext");
    p->addContext("TestNonDefaultContext");

    containers->Add(p);

    FairContainer* p2 = new FairContainer("TargetGeoPar", "Target geometry parameters", "GeometryParameterContext");
    p2->addContext("GeometryParameterContext");
    containers->Add(p2);
}

FairParSet* R3BPassiveContFact::createContainer(FairContainer* c)
{
    /** Calls the constructor of the corresponding parameter container.
     * For an actual context, which is not an empty string and not the default context
     * of this container, the name is concatinated with the context. */

    const char* name = c->GetName();
    FairParSet* p = NULL;
    if (strcmp(name, "R3BGeoPassivePar") == 0)
    {
        p = new R3BGeoPassivePar(c->getConcatName().Data(), c->GetTitle(), c->getContext());
    }
    if (strcmp(name, "TargetGeoPar") == 0)
    {
        p = new R3BTGeoPar(c->getConcatName().Data(), c->GetTitle(), c->getContext());
    }
    return p;
}

ClassImp(R3BPassiveContFact);
