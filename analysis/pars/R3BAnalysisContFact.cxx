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

#include "R3BAnalysisContFact.h"

#include "FairLogger.h"
#include "FairParAsciiFileIo.h"
#include "FairParRootFileIo.h"
#include "FairRuntimeDb.h"

#include "R3BIncomingIDPar.h"
#include "R3BTGeoPar.h"
#include "TClass.h"

static R3BAnalysisContFact gR3BAnalysisContFact;

R3BAnalysisContFact::R3BAnalysisContFact()
{
    // Constructor (called when the library is loaded)
    fName = "R3BAnalysisContFact";
    fTitle = "Factory for parameter containers needed for the online and offline analysis";
    setAllContainers();
    FairRuntimeDb::instance()->addContFactory(this);
}

void R3BAnalysisContFact::setAllContainers()
{
    // Creates the Container objects with all accepted contexts and adds them to
    // the list of containers for the STS library.

    FairContainer* p1 = new FairContainer("IncomingIDPar", "Incoming ID parameters", "IncomingIDParContext");
    p1->addContext("IncomingIDParContext");
    containers->Add(p1);
}

FairParSet* R3BAnalysisContFact::createContainer(FairContainer* c)
{
    // Trals the constructor of the corresponding parameter container.
    // For an actual context, which is not an empty string and not the default context
    // of this container, the name is concatinated with the context.

    const char* name = c->GetName();
    LOG(info) << "R3BIncomingIDContFact::Create container name: " << name;
    FairParSet* p = 0;
    if (strcmp(name, "IncomingIDPar") == 0)
    {
        p = new R3BIncomingIDPar(c->getConcatName().Data(), c->GetTitle(), c->getContext());
    }
    return p;
}

ClassImp(R3BAnalysisContFact);
