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

#include "R3BCalifaContFact.h"

#include "FairLogger.h"
#include "FairParAsciiFileIo.h"
#include "FairParRootFileIo.h"
#include "FairRuntimeDb.h"

#include "R3BCalifaCrystalCalPar.h"
#include "R3BCalifaCrystalPars4Sim.h"
#include "R3BCalifaMappingPar.h"
#include "R3BCalifaTotCalPar.h"
#include "R3BTGeoPar.h"
#include "TClass.h"

static R3BCalifaContFact gR3BCalifaContFact;

R3BCalifaContFact::R3BCalifaContFact()
{
    // Constructor (called when the library is loaded)
    fName = "R3BCalifaContFact";
    fTitle = "Factory for parameter containers in libR3BCalifa";
    setAllContainers();
    FairRuntimeDb::instance()->addContFactory(this);
}

void R3BCalifaContFact::setAllContainers()
{
    // Creates the Container objects with all accepted contexts and adds them
    // to the list of containers for the CALIFA library.

    FairContainer* p1 =
        new FairContainer("califaCrystalCalPar", "Califa Calibration Parameters", "CalifaCalParContext");
    p1->addContext("CalifaCalParContext");
    containers->Add(p1);

    FairContainer* p2 = new FairContainer("califaMappingPar", "Califa Mapping Parameters", "CalifaMappingContext");
    p2->addContext("CalifaMappingContext");
    containers->Add(p2);

    FairContainer* p3 =
        new FairContainer("califaCrystalPars4Sim", "Califa Crystal Parameters for Sim", "CalifaSimParContext");
    p3->addContext("CalifaSimParContext");
    containers->Add(p3);

    FairContainer* p4 = new FairContainer("CalifaTotCalPar", "Califa Tot Cal parameters", "CalifaTotCalParContext");
    p4->addContext("CalifaTotCalParContext");
    containers->Add(p4);

    FairContainer* p5 = new FairContainer("CalifaGeoPar", "Califa geometry parameters", "GeometryParameterContext");
    p5->addContext("GeometryParameterContext");
    containers->Add(p5);
}

FairParSet* R3BCalifaContFact::createContainer(FairContainer* c)
{
    /** Calls the constructor of the corresponding parameter container.
     * For an actual context, which is not an empty string and not the default context
     * of this container, the name is concatinated with the context.
     */

    const char* name = c->GetName();
    LOG(info) << "R3BCalifaContFact: Create container name: " << name;
    FairParSet* p = 0;
    if (strcmp(name, "califaCrystalCalPar") == 0)
    {
        p = new R3BCalifaCrystalCalPar(c->getConcatName().Data(), c->GetTitle(), c->getContext());
    }
    else if (strcmp(name, "califaMappingPar") == 0)
    {
        p = new R3BCalifaMappingPar(c->getConcatName().Data(), c->GetTitle(), c->getContext());
    }
    else if (strcmp(name, "califaCrystalPars4Sim") == 0)
    {
        p = new R3BCalifaCrystalPars4Sim(c->getConcatName().Data(), c->GetTitle(), c->getContext());
    }
    else if (strcmp(name, "CalifaTotCalPar") == 0)
    {
        p = new R3BCalifaTotCalPar(c->getConcatName().Data(), c->GetTitle(), c->getContext());
    }
    else if (strcmp(name, "CalifaGeoPar") == 0)
    {
        p = new R3BTGeoPar(c->getConcatName().Data(), c->GetTitle(), c->getContext());
    }
    return p;
}

ClassImp(R3BCalifaContFact);
