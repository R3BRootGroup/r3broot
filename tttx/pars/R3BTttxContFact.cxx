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
// -----            R3BTttxContFact source file              -----
// -----    Created 18/01/2024 by J.L. Rodriguez-Sanchez     -----
//----------------------------------------------------------------

#include <FairLogger.h>
#include <FairRuntimeDb.h>

#include <TClass.h>
#include <string>

#include "R3BTGeoPar.h"
#include "R3BTttxContFact.h"
// #include "R3BTttxMappingPar.h"
#include "R3BTttxHitPar.h"
#include "R3BTttxStripCalPar.h"

static R3BTttxContFact gR3BTttxContFact;

R3BTttxContFact::R3BTttxContFact()
{
    // Constructor (called when the library is loaded)
    fName = "R3BTttxContFact";
    fTitle = "Factory for parameter containers in libR3BTttx";
    setAllContainers();
    FairRuntimeDb::instance()->addContFactory(this);
}

void R3BTttxContFact::setAllContainers()
{
    // Creates the Container objects with all accepted contexts and adds them to
    // the list of containers for the TTTX library.

    FairContainer* p1 = new FairContainer("tttxGeoPar", "TTTX Geometry Parameters", "GeometryParameterContext");
    p1->addContext("GeometryParameterContext");
    containers->Add(p1);

    // FairContainer* p2 = new FairContainer("tttxMappingPar", "TTTX Mapping Parameters", "TttxMappingParContext");
    // p2->addContext("TttxMappingParContext");
    // containers->Add(p2);

    FairContainer* p3 = new FairContainer("tttxStripCalPar", "TTTX StripCal Parameters", "TttxCalParContext");
    p3->addContext("TttxCalParContext");
    containers->Add(p3);

    FairContainer* p4 = new FairContainer("tttxHitPar", "TTTX Hit Parameters", "TttxHitParContext");
    p4->addContext("TttxHitParContext");
    containers->Add(p4);
}

FairParSet* R3BTttxContFact::createContainer(FairContainer* c)
{
    // For an actual context, which is not an empty string and not the default context
    // of this container, the name is concatinated with the context.

    const std::string name(c->GetName());
    LOG(info) << "R3BTttxContFact: Create container name: " << name.c_str();

    FairParSet* p = nullptr;
    if (name == "tttxGeoPar")
    {
        p = new R3BTGeoPar(c->getConcatName().Data(), c->GetTitle(), c->getContext());
    }
    else if (name == "tttxMappingPar")
    {
        // p = new R3BTttxMappingPar(c->getConcatName().Data(), c->GetTitle(), c->getContext());
    }
    else if (name == "tttxStripCalPar")
    {
        p = new R3BTttxStripCalPar(c->getConcatName().Data(), c->GetTitle(), c->getContext());
    }
    else if (name == "tttxHitPar")
    {
        p = new R3BTttxHitPar(c->getConcatName().Data(), c->GetTitle(), c->getContext());
    }
    return p;
}

ClassImp(R3BTttxContFact)
