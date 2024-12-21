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

/////////////////////////////////////////////////////////////
//
//  R3BPtofParFact
//
//  Factory for the parameter containers
//
/////////////////////////////////////////////////////////////
#include "R3BPtofParFact.h"

#include "R3BPtofHitPar.h"

#include "FairLogger.h"
#include "FairParSet.h"    // for FairParSet
#include "FairRuntimeDb.h" // for FairRuntimeDb

#include "TList.h"   // for TList
#include "TString.h" // for TString

// #include <string.h> // for strcmp, NULL

using namespace std;

static R3BPtofParFact gR3BPtofParFact;

R3BPtofParFact::R3BPtofParFact()
{
    // Constructor (called when the library is loaded)
    fName = "R3BPtofParFact";
    fTitle = "Tutorial factory for parameter containers";
    setAllContainers();
    FairRuntimeDb::instance()->addContFactory(this);
}

void R3BPtofParFact::setAllContainers()
{
    /** Creates the Container objects with all accepted contexts and adds them to
     *  the list of containers.*/

    FairContainer* p1 = new FairContainer("PtofHitPar", "Ptof Hit Parameters", "TestDefaultContext");
    p1->addContext("TestNonDefaultContext");
    containers->Add(p1);
}

FairParSet* R3BPtofParFact::createContainer(FairContainer* c)
{
    /** Calls the constructor of the corresponding parameter container.
     * For an actual context, which is not an empty string and not the default context
     * of this container, the name is concatinated with the context. */

    const char* name = c->GetName();
    LOG(info) << "R3BPtofParFact::createContainer : " << name;
    FairParSet* p = NULL;

    if (strcmp(name, "PtofHitPar") == 0)
    {
        p = new R3BPtofHitPar(c->getConcatName().Data(), c->GetTitle(), c->getContext());
    }

    return p;
}

ClassImp(R3BPtofParFact)
