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

#include "R3BNeulandContFact.h"
#include "FairRuntimeDb.h"
#include "R3BNeulandGeoPar.h"

static R3BNeulandContFact gR3BNeulandContFact;

R3BNeulandContFact::R3BNeulandContFact()
{
    fName = "R3BNeulandContFact";
    fTitle = "Factory for Neuland parameter containers";
    setAllContainers();
    FairRuntimeDb::instance()->addContFactory(this);
}

void R3BNeulandContFact::setAllContainers()
{
    FairContainer* p1 = new FairContainer("R3BNeulandGeoPar", "Neuland Geometry Parameters", "TestDefaultContext");
    p1->addContext("TestNonDefaultContext");
    containers->Add(p1);
}

FairParSet* R3BNeulandContFact::createContainer(FairContainer* c)
{
    const char* name = c->GetName();

    if (strcmp(name, "R3BNeulandGeoPar") == 0)
    {
        return new R3BNeulandGeoPar(c->getConcatName().Data(), c->GetTitle(), c->getContext());
    }

    return nullptr;
}

ClassImp(R3BNeulandContFact)
