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
//  R3BFiberContFact
//
//  Factory for the parameter containers in libR3BFiber
//
/////////////////////////////////////////////////////////////

#include "FairLogger.h"
#include "FairRuntimeDb.h"

#include "R3BBunchedFiberHitPar.h"
#include "R3BFiberContFact.h"
#include "R3BFiberMAPMTHitPar.h"
#include "R3BFiberMappingPar.h"
#include "R3BLogger.h"
#include "R3BTGeoPar.h"

#include "TClass.h"

static R3BFiberContFact gR3BFiberContFact;

R3BFiberContFact::R3BFiberContFact()
{
    // Constructor (called when the library is loaded)
    fName = "R3BFiberContFact";
    fTitle = "Factory for parameter containers in libR3BFiber";
    setAllContainers();
    FairRuntimeDb::instance()->addContFactory(this);
}

void R3BFiberContFact::setAllContainers()
{
    /** Creates the Container objects with all accepted contexts and adds them to
     *  the list of containers for the STS library.*/

    FairContainer* p1 = new FairContainer("Fi3aHitPar", "Fi3a Hit Parameters", "TestDefaultContext");
    p1->addContext("TestNonDefaultContext");
    containers->Add(p1);

    FairContainer* p2 = new FairContainer("Fi3bHitPar", "Fi3b Hit Parameters", "TestDefaultContext");
    p2->addContext("TestNonDefaultContext");
    containers->Add(p2);

    FairContainer* p3 = new FairContainer("Fi7HitPar", "Fi7 Hit Parameters", "TestDefaultContext");
    p3->addContext("TestNonDefaultContext");
    containers->Add(p3);

    FairContainer* p4 = new FairContainer("Fi8HitPar", "Fi8 Hit Parameters", "TestDefaultContext");
    p4->addContext("TestNonDefaultContext");
    containers->Add(p4);

    FairContainer* p5 = new FairContainer("Fi10HitPar", "Fi10 Hit Parameters", "TestDefaultContext");
    p5->addContext("TestNonDefaultContext");
    containers->Add(p5);

    FairContainer* p6 = new FairContainer("Fi11HitPar", "Fi11 Hit Parameters", "TestDefaultContext");
    p6->addContext("TestNonDefaultContext");
    containers->Add(p6);

    FairContainer* p7 = new FairContainer("Fi12HitPar", "Fi12 Hit Parameters", "TestDefaultContext");
    p7->addContext("TestNonDefaultContext");

    containers->Add(p7);
    FairContainer* p8 = new FairContainer("Fi13HitPar", "Fi13 Hit Parameters", "TestDefaultContext");
    p8->addContext("TestNonDefaultContext");
    containers->Add(p8);

    FairContainer* p9 = new FairContainer("Fi1aHitPar", "Fi1a Hit Parameters", "TestDefaultContext");
    p9->addContext("TestNonDefaultContext");
    containers->Add(p9);

    FairContainer* p10 = new FairContainer("Fi1bHitPar", "Fi1b Hit Parameters", "TestDefaultContext");
    p10->addContext("TestNonDefaultContext");
    containers->Add(p10);

    FairContainer* p11 = new FairContainer("Fi30HitPar", "Fi30 Hit Parameters", "TestDefaultContext");
    p11->addContext("TestNonDefaultContext");
    containers->Add(p11);

    FairContainer* p12 = new FairContainer("Fi31HitPar", "Fi31 Hit Parameters", "TestDefaultContext");
    p12->addContext("TestNonDefaultContext");
    containers->Add(p12);

    FairContainer* p13 = new FairContainer("Fi32HitPar", "Fi32 Hit Parameters", "TestDefaultContext");
    p13->addContext("TestNonDefaultContext");
    containers->Add(p13);

    FairContainer* p14 = new FairContainer("Fi33HitPar", "Fi33 Hit Parameters", "TestDefaultContext");
    p14->addContext("TestNonDefaultContext");
    containers->Add(p14);

    FairContainer* p15 = new FairContainer("Fi23aHitPar", "Fi23a Hit Parameters", "TestDefaultContext");
    p15->addContext("TestNonDefaultContext");
    containers->Add(p15);

    FairContainer* p16 = new FairContainer("Fi23bHitPar", "Fi23b Hit Parameters", "TestDefaultContext");
    p16->addContext("TestNonDefaultContext");
    containers->Add(p16);

    auto p17 = new FairContainer("Fi30MappingPar", "Fi30 Mapping Parameters", "TestDefaultContext");
    p17->addContext("TestNonDefaultContext");
    containers->Add(p17);

    auto p18 = new FairContainer("Fi31MappingPar", "Fi31 Mapping Parameters", "TestDefaultContext");
    p18->addContext("TestNonDefaultContext");
    containers->Add(p18);

    auto p19 = new FairContainer("Fi32MappingPar", "Fi32 Mapping Parameters", "TestDefaultContext");
    p19->addContext("TestNonDefaultContext");
    containers->Add(p19);

    auto p20 = new FairContainer("Fi33MappingPar", "Fi33 Mapping Parameters", "TestDefaultContext");
    p20->addContext("TestNonDefaultContext");
    containers->Add(p20);

    auto p21 = new FairContainer("Fi10GeoPar", "Fi10 geometry parameters", "fi10GeoContext");
    p21->addContext("fi10GeoContext");
    containers->Add(p21);

    auto p22 = new FairContainer("Fi11GeoPar", "Fi11 geometry parameters", "fi11GeoContext");
    p22->addContext("fi11GeoContext");
    containers->Add(p22);

    auto p23 = new FairContainer("Fi12GeoPar", "Fi12 geometry parameters", "fi12GeoContext");
    p23->addContext("fi12GeoContext");
    containers->Add(p23);

    auto p24 = new FairContainer("Fi13GeoPar", "Fi13 geometry parameters", "fi13GeoContext");
    p24->addContext("fi13GeoContext");
    containers->Add(p24);

    auto p25 = new FairContainer("Fi30GeoPar", "Fi30 geometry parameters", "fi30GeoContext");
    p25->addContext("fi30GeoContext");
    containers->Add(p25);

    auto p26 = new FairContainer("Fi31GeoPar", "Fi31 geometry parameters", "fi31GeoContext");
    p26->addContext("fi31GeoContext");
    containers->Add(p26);

    auto p27 = new FairContainer("Fi32GeoPar", "Fi32 geometry parameters", "fi32GeoContext");
    p27->addContext("fi32GeoContext");
    containers->Add(p27);

    auto p28 = new FairContainer("Fi33GeoPar", "Fi33 geometry parameters", "fi33GeoContext");
    p28->addContext("fi33GeoContext");
    containers->Add(p28);

    auto p29 = new FairContainer("Fi23aGeoPar", "Fi23a geometry parameters", "fi23aGeoContext");
    p29->addContext("fi23aGeoContext");
    containers->Add(p29);

    auto p30 = new FairContainer("Fi23bGeoPar", "Fi23b geometry parameters", "fi23bGeoContext");
    p30->addContext("fi23bGeoContext");
    containers->Add(p30);

    auto p31 = new FairContainer("Fi10MappingPar", "Fi10 Mapping Parameters", "TestDefaultContext");
    p31->addContext("TestNonDefaultContext");
    containers->Add(p31);

    auto p32 = new FairContainer("Fi11MappingPar", "Fi11 Mapping Parameters", "TestDefaultContext");
    p32->addContext("TestNonDefaultContext");
    containers->Add(p32);

    auto p33 = new FairContainer("Fi12MappingPar", "Fi12 Mapping Parameters", "TestDefaultContext");
    p33->addContext("TestNonDefaultContext");
    containers->Add(p33);

    auto p34 = new FairContainer("Fi13MappingPar", "Fi13 Mapping Parameters", "TestDefaultContext");
    p34->addContext("TestNonDefaultContext");
    containers->Add(p34);

    auto p35 = new FairContainer("Fi7GeoPar", "Fi7 geometry parameters", "fi7GeoContext");
    p35->addContext("fi7GeoContext");
    containers->Add(p35);

    auto p36 = new FairContainer("Fi8GeoPar", "Fi8 geometry parameters", "fi8GeoContext");
    p36->addContext("fi8GeoContext");
    containers->Add(p36);
}

FairParSet* R3BFiberContFact::createContainer(FairContainer* c)
{
    /** Fiberls the constructor of the corresponding parameter container.
     * For an actual context, which is not an empty string and not the default context
     * of this container, the name is concatinated with the context. */

    const char* name = c->GetName();
    R3BLOG(info, name);

    FairParSet* p = 0;
    if (strcmp(name, "Fi7GeoPar") == 0 || strcmp(name, "Fi8GeoPar") == 0 || strcmp(name, "Fi10GeoPar") == 0 ||
        strcmp(name, "Fi11GeoPar") == 0 || strcmp(name, "Fi12GeoPar") == 0 || strcmp(name, "Fi13GeoPar") == 0 ||
        strcmp(name, "Fi23aGeoPar") == 0 || strcmp(name, "Fi23bGeoPar") == 0 || strcmp(name, "Fi30GeoPar") == 0 ||
        strcmp(name, "Fi31GeoPar") == 0 || strcmp(name, "Fi32GeoPar") == 0 || strcmp(name, "Fi33GeoPar") == 0)
    {
        p = new R3BTGeoPar(c->getConcatName().Data(), c->GetTitle(), c->getContext());
    }
    else if (strcmp(name, "Fi10MappingPar") == 0 || strcmp(name, "Fi11MappingPar") == 0 ||
             strcmp(name, "Fi12MappingPar") == 0 || strcmp(name, "Fi13MappingPar") == 0 ||
             strcmp(name, "Fi23aMappingPar") == 0 || strcmp(name, "Fi23bMappingPar") == 0 ||
             strcmp(name, "Fi30MappingPar") == 0 || strcmp(name, "Fi31MappingPar") == 0 ||
             strcmp(name, "Fi32MappingPar") == 0 || strcmp(name, "Fi33MappingPar") == 0)
    {
        p = new R3BFiberMappingPar(c->getConcatName().Data(), c->GetTitle(), c->getContext());
    }
    else if (strcmp(name, "Fi30HitPar") == 0 || strcmp(name, "Fi31HitPar") == 0 || strcmp(name, "Fi32HitPar") == 0 ||
             strcmp(name, "Fi33HitPar") == 0 || strcmp(name, "Fi23aHitPar") == 0 || strcmp(name, "Fi23bHitPar") == 0)
    {
        p = new R3BFiberMAPMTHitPar(c->getConcatName().Data(), c->GetTitle(), c->getContext());
    }
    else
    {
        p = new R3BBunchedFiberHitPar(c->getConcatName().Data(), c->GetTitle(), c->getContext());
    }
    return p;
}

ClassImp(R3BFiberContFact);
