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

//*-- AUTHOR : D. Kresan
//*-- Created : 18/05/2015

/////////////////////////////////////////////////////////////
//
//  R3BTofdContFact
//
//  Factory for the parameter containers in libR3BTofd
//
/////////////////////////////////////////////////////////////

#include "R3BTofdContFact.h"

#include "R3BGeoTofdPar.h"
#include "R3BTGeoPar.h"

#include "FairParAsciiFileIo.h"
#include "FairParRootFileIo.h"
#include "FairRuntimeDb.h"

#include "TClass.h"

#include <iomanip>
#include <iostream>

using std::cout;
using std::endl;

static R3BTofdContFact gR3BTofdContFact;

R3BTofdContFact::R3BTofdContFact()
{
    // Constructor (called when the library is loaded)
    fName = "R3BTofdContFact";
    fTitle = "Factory for parameter containers in libR3BTofd";
    setAllContainers();
    FairRuntimeDb::instance()->addContFactory(this);
}

void R3BTofdContFact::setAllContainers()
{
    /** Creates the Container objects with all accepted contexts and adds them to
     *  the list of containers for the STS library.*/
    /*
     FairContainer* p1= new FairContainer("R3BTofdDigiPar",
                                           "Sts Digitisation Parameters",
                                           "TestDefaultContext");
     p1->addContext("TestNonDefaultContext");

     FairContainer* p2= new FairContainer("CbmGeoStsPar",
                                           "Sts Geometry Parameters",
                                           "TestDefaultContext");
     p2->addContext("TestNonDefaultContext");

     containers->Add(p1);
     containers->Add(p2);
     */

    // FairContainer* p2 = new FairContainer("tofdGeoPar", "Tofd geometry parameters", "TestDefaultContext");
    // p2->addContext("TestNonDefaultContext");
    //  containers->Add(p2);
}

FairParSet* R3BTofdContFact::createContainer(FairContainer* c)
{
    /** Tofdls the constructor of the corresponding parameter container.
     * For an actual context, which is not an empty string and not the default context
     * of this container, the name is concatinated with the context. */

    const char* name = c->GetName();
    cout << " -I container name " << name << endl;
    FairParSet* p = 0;
    /*if (strcmp(name,"R3BTofdDigiPar")==0) {
      p=new R3BTofdDigiPar(c->getConcatName().Data(),c->GetTitle(),c->getContext());
    }
    if (strcmp(name,"CbmGeoStsPar")==0) {
      p=new CbmGeoStsPar(c->getConcatName().Data(),c->GetTitle(),c->getContext());
    }*/
    /*
        if (strcmp(name, "tofdGeoPar") == 0)
        {
            p = new R3BTGeoPar(c->getConcatName().Data(), c->GetTitle(), c->getContext());
        }*/

    return p;
}

ClassImp(R3BTofdContFact);
