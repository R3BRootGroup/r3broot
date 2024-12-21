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

// ---------------------------------------------------------------
// -----            R3BSsdContFact source file               -----
// -----    Created 29/05/2018 by J.L. Rodriguez-Sanchez     -----
// ----- Factory for the parameter containers in libR3BSsd   -----
//----------------------------------------------------------------

#include "FairLogger.h"
#include "FairRuntimeDb.h"

#include "R3BGeoTraPar.h"
#include "R3BSsdContFact.h"

#include "R3BAmsMappingPar.h"
#include "R3BAmsStripCalPar.h"

#include "R3BFootCalPar.h"
#include "R3BFootMappingPar.h"

#include "R3BTGeoPar.h"

#include "TClass.h"

static R3BSsdContFact gR3BSsdContFact;

R3BSsdContFact::R3BSsdContFact()
{
    // Constructor (called when the library is loaded)
    fName = "R3BSsdContFact";
    fTitle = "Factory for parameter containers in libR3BSsd";
    setAllContainers();
    FairRuntimeDb::instance()->addContFactory(this);
}

void R3BSsdContFact::setAllContainers()
{
    // Creates the Container objects with all accepted contexts and adds them to
    // the list of containers for the SSD library.

    FairContainer* p1 = new FairContainer("amsMappingPar", "AMS Mapping Parameters", "AmsMappingParContext");
    p1->addContext("AmsMappingParContext");
    containers->Add(p1);

    FairContainer* p2 = new FairContainer("amsStripCalPar", "AMS StripCal Parameters", "AmsCalParContext");
    p2->addContext("AmsCalParContext");
    containers->Add(p2);

    FairContainer* p3 = new FairContainer("amsGeoPar", "AMS Geometry Parameters", "GeometryParameterContext");
    p3->addContext("GeometryParameterContext");
    containers->Add(p3);

    FairContainer* p4 = new FairContainer("footMappingPar", "FOOT Mapping Parameters", "FootMappingParContext");
    p4->addContext("FootMappingParContext");
    containers->Add(p4);

    FairContainer* p5 = new FairContainer("footCalPar", "FOOT StripCal Parameters", "FootCalParContext");
    p5->addContext("FootCalParContext");
    containers->Add(p5);

    FairContainer* p6 = new FairContainer("footGeoPar", "FOOT Geometry Parameters", "GeometryParameterContext");
    p6->addContext("GeometryParameterContext");
    containers->Add(p6);
}

FairParSet* R3BSsdContFact::createContainer(FairContainer* c)
{
    // For an actual context, which is not an empty string and not the default context
    // of this container, the name is concatinated with the context.

    const char* name = c->GetName();
    LOG(info) << "R3BSsdContFact: Create container name: " << name;
    FairParSet* p = 0;
    if (strcmp(name, "amsMappingPar") == 0)
    {
        p = new R3BAmsMappingPar(c->getConcatName().Data(), c->GetTitle(), c->getContext());
    }
    else if (strcmp(name, "amsStripCalPar") == 0)
    {
        p = new R3BAmsStripCalPar(c->getConcatName().Data(), c->GetTitle(), c->getContext());
    }
    else if (strcmp(name, "footMappingPar") == 0)
    {
        p = new R3BFootMappingPar(c->getConcatName().Data(), c->GetTitle(), c->getContext());
    }
    else if (strcmp(name, "footCalPar") == 0)
    {
        p = new R3BFootCalPar(c->getConcatName().Data(), c->GetTitle(), c->getContext());
    }
    else if (strcmp(name, "amsGeoPar") == 0 || strcmp(name, "footGeoPar") == 0)
    {
        p = new R3BTGeoPar(c->getConcatName().Data(), c->GetTitle(), c->getContext());
    }
    return p;
}

ClassImp(R3BSsdContFact);
