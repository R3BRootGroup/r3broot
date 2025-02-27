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

// ----------------------------------------------------------------
// -----            R3BActafContFact source file              -----
// -----    Created 27/02/2025 by J.L. Rodriguez-Sanchez      -----
// ----- Factory for the parameter containers in libR3BActaf  -----
//-----------------------------------------------------------------

#include <FairRuntimeDb.h>

#include "R3BActafContFact.h"
#include "R3BLogger.h"
#include "R3BTGeoPar.h"

#include <TClass.h>

static R3BActafContFact gR3BActafContFact;

R3BActafContFact::R3BActafContFact()
{
    // Constructor (called when the library is loaded)
    fName = "R3BActafContFact";
    fTitle = "Factory for parameter containers in libR3BActaf";
    setAllContainers();
    FairRuntimeDb::instance()->addContFactory(this);
}

void R3BActafContFact::setAllContainers()
{
    // Creates the container objects with all accepted contexts and adds them to
    // the list of containers for the Actaf library.

    auto* p1 = new FairContainer("ActafGeoPar", "ACTAF Geometry Parameters", "GeometryParContext");
    p1->addContext("GeometryParContext");
    containers->Add(p1);
}

FairParSet* R3BActafContFact::createContainer(FairContainer* c)
{
    // For an actual context, which is not an empty string and not the default context
    // of this container, the name is concatinated with the context.
    const std::string name(c->GetName());
    R3BLOG(info, "Create container name: " << name.c_str());

    FairParSet* p = nullptr;
    if (name == "ActafGeoPar")
    {
        p = new R3BTGeoPar(c->getConcatName().Data(), c->GetTitle(), c->getContext());
    }
    return p;
}

ClassImp(R3BActafContFact)
