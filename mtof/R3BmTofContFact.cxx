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

//*-- AUTHOR : Denis Bertini
//*-- Created : 20/03/2009

/////////////////////////////////////////////////////////////
//
//  R3BmTofContFact
//
//  Factory for the parameter containers in libR3BmTof
//
/////////////////////////////////////////////////////////////

#include "R3BmTofContFact.h"

// #include "R3BmTofParRootFileIo.h"
// #include "R3BmTofParAsciiFileIo.h"
#include "R3BGeomTofPar.h"
// #include "R3BmTofDigiPar.h"
#include "R3BTGeoPar.h"

#include "FairRuntimeDb.h"
// #include "CbmParTest.h"
#include "FairParAsciiFileIo.h"
#include "FairParRootFileIo.h"

#include "TClass.h"

#include <iomanip>
#include <iostream>

using std::cout;
using std::endl;

ClassImp(R3BmTofContFact)

    static R3BmTofContFact gR3BmTofContFact;

R3BmTofContFact::R3BmTofContFact()
{
    // Constructor (called when the library is loaded)
    fName = "R3BmTofContFact";
    fTitle = "Factory for parameter containers in libR3BmTof";
    setAllContainers();
    FairRuntimeDb::instance()->addContFactory(this);
}

void R3BmTofContFact::setAllContainers()
{
    /** Creates the Container objects with all accepted contexts and adds them to
     *  the list of containers for the STS library.*/

    //    FairContainer* p1= new FairContainer("R3BmTofDigiPar",
    //                                          "Sts Digitisation Parameters",
    //                                          "TestDefaultContext");
    //    p1->addContext("TestNonDefaultContext");
    //
    //    containers->Add(p1);

    FairContainer* p2 = new FairContainer("mTofGeoPar", "mTOF geometry parameters", "TestDefaultContext");
    p2->addContext("TestNonDefaultContext");
    containers->Add(p2);
}

FairParSet* R3BmTofContFact::createContainer(FairContainer* c)
{
    /** mTofls the constructor of the corresponding parameter container.
     * For an actual context, which is not an empty string and not the default context
     * of this container, the name is concatinated with the context. */

    const char* name = c->GetName();
    cout << " -I container name " << name << endl;
    FairParSet* p = 0;

    //  if (strcmp(name,"R3BmTofDigiPar")==0) {
    //    p=new R3BmTofDigiPar(c->getConcatName().Data(),c->GetTitle(),c->getContext());
    //  }

    if (strcmp(name, "mTofGeoPar") == 0)
    {
        p = new R3BTGeoPar(c->getConcatName().Data(), c->GetTitle(), c->getContext());
    }

    return p;

    return NULL;
}

void R3BmTofContFact::activateParIo(FairParIo* io)
{
    // activates the input/output class for the parameters
    // needed by the Sts
    /*
    if (strcmp(io->IsA()->GetName(),"FairParRootFileIo")==0) {
      R3BmTofParRootFileIo* p=new R3BmTofParRootFileIo(((FairParRootFileIo*)io)->getParRootFile());
      io->setDetParIo(p);
    }
    if (strcmp(io->IsA()->GetName(),"FairParAsciiFileIo")==0) {
      R3BmTofParAsciiFileIo* p=new R3BmTofParAsciiFileIo(((FairParAsciiFileIo*)io)->getFile());
      io->setDetParIo(p);
      }
    */
}
