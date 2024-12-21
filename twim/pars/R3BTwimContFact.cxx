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

// --------------------------------------------------------------------
// -----           R3BTwimContFact source file                    -----
// -----     Created 24/11/19 by J.L. Rodriguez-Sanchez           -----
// --------------------------------------------------------------------
//
//  R3BTwimContFact
//
//  Factory for the parameter containers in libR3BTwim
//

#include "R3BTwimContFact.h"

#include <FairLogger.h>
#include <FairRuntimeDb.h>

#include "R3BLogger.h"
#include "R3BTGeoPar.h"
#include "R3BTwimCalPar.h"
#include "R3BTwimHitPar.h"

#include <TClass.h>

static R3BTwimContFact gR3BTwimContFact;

R3BTwimContFact::R3BTwimContFact()
{
    // Constructor (called when the library is loaded)
    fName = "R3BTwimContFact";
    fTitle = "Factory for parameter containers in libR3BTwim";
    setAllContainers();
    FairRuntimeDb::instance()->addContFactory(this);
}

void R3BTwimContFact::setAllContainers()
{
    // Creates the Container objects with all accepted contexts and adds them to
    // the list of containers for the STS library.

    FairContainer* p1 = new FairContainer("twimCalPar", "TWIM Cal Parameters", "TwimCalParContext");
    p1->addContext("TwimCalParContext");
    containers->Add(p1);

    FairContainer* p2 = new FairContainer("twimHitPar", "TWIM Hit Parameters", "TwimHitParContext");
    p2->addContext("TwimHitParContext");
    containers->Add(p2);

    FairContainer* p3 = new FairContainer("TwimGeoPar", "Twim-Music geometry parameters", "GeometryParameterContext");
    p3->addContext("GeometryParameterContext");
    containers->Add(p3);
}

FairParSet* R3BTwimContFact::createContainer(FairContainer* c)
{
    // Trals the constructor of the corresponding parameter container.
    // For an actual context, which is not an empty string and not the default context
    // of this container, the name is concatinated with the context.

    const char* name = c->GetName();
    R3BLOG(info, "Create container name: " << name);
    FairParSet* p = 0;
    if (strcmp(name, "twimCalPar") == 0)
    {
        p = new R3BTwimCalPar(c->getConcatName().Data(), c->GetTitle(), c->getContext());
    }
    else if (strcmp(name, "twimHitPar") == 0)
    {
        p = new R3BTwimHitPar(c->getConcatName().Data(), c->GetTitle(), c->getContext());
    }
    else if (strcmp(name, "TwimGeoPar") == 0)
    {
        p = new R3BTGeoPar(c->getConcatName().Data(), c->GetTitle(), c->getContext());
    }
    return p;
}

ClassImp(R3BTwimContFact)
