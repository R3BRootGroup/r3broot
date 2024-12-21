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

/********************************************************************************
 *    Copyright (C) 2014 GSI Helmholtzzentrum fuer Schwerionenforschung GmbH    *
 *                                                                              *
 *              This software is distributed under the terms of the             *
 *         GNU Lesser General Public Licence version 3 (LGPL) version 3,        *
 *                  copied verbatim in the file "LICENSE"                       *
 ********************************************************************************/

// -------------------------------------------------------------------------
// -----               R3BPspxContFact source file                   -----
// -----                 Created 16/05/12  by I.Syndikus               -----
// -------------------------------------------------------------------------

#include "R3BPspxContFact.h"

#include "R3BPspxCalPar.h"
#include "R3BPspxHitPar.h"
#include "R3BPspxMappedPar.h"
#include "R3BPspxPrecalPar.h"

#include "R3BTGeoPar.h"

#include "FairRuntimeDb.h"

#include <iostream>

ClassImp(R3BPspxContFact)

    static R3BPspxContFact gR3BPspxContFact;

R3BPspxContFact::R3BPspxContFact()
{
    fName = "R3BPspxContFact";
    fTitle = "Factory for parameter containers for PSPX";
    setAllContainers();
    FairRuntimeDb::instance()->addContFactory(this);
}

void R3BPspxContFact::setAllContainers()
{
    /** Creates the Container objects with all accepted
        contexts and adds them to
        the list of containers for unpacking/analysing PSPX detectors.
    */

    FairContainer* p1 = new FairContainer("R3BPspxMappedPar", "R3BPspx Mapped Parameters", "PspxDefaultContext");
    p1->addContext("PspxNonDefaultContext");

    containers->Add(p1);

    FairContainer* p2 = new FairContainer("R3BPspxPrecalPar", "R3BPspx Precal Parameters", "PspxDefaultContext");
    p2->addContext("PspxNonDefaultContext");

    containers->Add(p2);

    FairContainer* p3 = new FairContainer("R3BPspxCalPar", "R3BPspx Cal Parameters", "PspxDefaultContext");
    p3->addContext("PspxNonDefaultContext");

    containers->Add(p3);

    FairContainer* p4 = new FairContainer("R3BPspxHitPar", "R3BPspx Hit Parameters", "PspxDefaultContext");
    p4->addContext("PspxNonDefaultContext");

    containers->Add(p4);

    FairContainer* p5 = new FairContainer("PspGeoPar", "PSP geometry parameters", "TestDefaultContext");
    p5->addContext("TestNonDefaultContext");

    containers->Add(p5);
}

FairParSet* R3BPspxContFact::createContainer(FairContainer* c)
{
    /** Calls the constructor of the corresponding parameter container.
        For an actual context, which is not an empty string and not
        the default context
        of this container, the name is concatinated with the context.
    */
    const char* name = c->GetName();
    FairParSet* p = NULL;
    if (strcmp(name, "R3BPspxMappedPar") == 0)
    {
        p = new R3BPspxMappedPar(c->getConcatName().Data(), c->GetTitle(), c->getContext());
    }
    if (strcmp(name, "R3BPspxPrecalPar") == 0)
    {
        p = new R3BPspxPrecalPar(c->getConcatName().Data(), c->GetTitle(), c->getContext());
    }
    if (strcmp(name, "R3BPspxCalPar") == 0)
    {
        p = new R3BPspxCalPar(c->getConcatName().Data(), c->GetTitle(), c->getContext());
    }
    if (strcmp(name, "R3BPspxHitPar") == 0)
    {
        p = new R3BPspxHitPar(c->getConcatName().Data(), c->GetTitle(), c->getContext());
    }
    if (strcmp(name, "PspGeoPar") == 0)
    {
        p = new R3BTGeoPar(c->getConcatName().Data(), c->GetTitle(), c->getContext());
    }

    return p;
}
