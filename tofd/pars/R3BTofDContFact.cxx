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

// ------------------------------------------------------------
// -----                 R3BTofDContFact                  -----
// -----    Created 18/03/22 by J.L. Rodriguez-Sanchez    -----
// ------------------------------------------------------------

#include "FairRuntimeDb.h"

#include "R3BLogger.h"
#include "R3BTGeoPar.h"
#include "R3BTofDContFact.h"
#include "R3BTofDHitPar.h"
#include "R3BTofDMappingPar.h"
#include "TClass.h"

#include <string.h>

static R3BTofDContFact gR3BTofDContFact;

R3BTofDContFact::R3BTofDContFact()
{
    // Constructor (called when the library is loaded)
    fName = "R3BTofDContFact";
    fTitle = "Factory for parameter containers in libR3BTofD";
    setAllContainers();
    FairRuntimeDb::instance()->addContFactory(this);
}

void R3BTofDContFact::setAllContainers()
{
    /** Creates the Container objects with all accepted contexts and adds them to
     *  the list of containers for the STS library.*/

    FairContainer* p1 = new FairContainer("tofdGeoPar", "Tofd geometry parameters", "tofdGeoContext");
    p1->addContext("tofdGeoContext");
    containers->Add(p1);

    FairContainer* p2 = new FairContainer("tofdMappingPar", "Tofd Mapping parameters", "tofdMappingContext");
    p2->addContext("tofdMappingContext");
    containers->Add(p2);

    FairContainer* p3 = new FairContainer("tofdHitPar", "Tofd Hit parameters", "tofdHitContext");
    p3->addContext("tofdHitContext");
    containers->Add(p3);
}

FairParSet* R3BTofDContFact::createContainer(FairContainer* c)
{
    /** Tofd the constructor of the corresponding parameter container.
     * For an actual context, which is not an empty string and not the default context
     * of this container, the name is concatinated with the context. */

    const char* name = c->GetName();
    R3BLOG(info, name);
    FairParSet* p = 0;
    if (strcmp(name, "tofdGeoPar") == 0)
    {
        p = new R3BTGeoPar(c->getConcatName().Data(), c->GetTitle(), c->getContext());
    }
    else if (strcmp(name, "tofdMappingPar") == 0)
    {
        p = new R3BTofDMappingPar(c->getConcatName().Data(), c->GetTitle(), c->getContext());
    }
    else if (strcmp(name, "tofdHitPar") == 0)
    {
        p = new R3BTofDHitPar(c->getConcatName().Data(), c->GetTitle(), c->getContext());
    }

    return p;
}

ClassImp(R3BTofDContFact);
