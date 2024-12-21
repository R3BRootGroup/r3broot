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
//  R3BTofContFact
//
//  Factory for the parameter containers in libR3BTof
//
/////////////////////////////////////////////////////////////

#include "R3BTofContFact.h"

// #include "R3BTofParRootFileIo.h"
// #include "R3BTofParAsciiFileIo.h"
#include "R3BGeoTofPar.h"
// #include "R3BTofDigiPar.h"
#include "R3BTofdHitPar.h"

#include "FairRuntimeDb.h"
// #include "CbmParTest.h"
#include "FairParAsciiFileIo.h"
#include "FairParRootFileIo.h"

#include "TClass.h"

#include <iomanip>
#include <iostream>

using std::cout;
using std::endl;

ClassImp(R3BTofContFact)

    static R3BTofContFact gR3BTofContFact;

R3BTofContFact::R3BTofContFact()
{
    // Constructor (called when the library is loaded)
    fName = "R3BTofContFact";
    fTitle = "Factory for parameter containers in libR3BTof";
    setAllContainers();
    FairRuntimeDb::instance()->addContFactory(this);
}

void R3BTofContFact::setAllContainers()
{
    /** Creates the Container objects with all accepted contexts and adds them to
     *  the list of containers for the STS library.*/
    /*
     FairContainer* p1= new FairContainer("R3BTofDigiPar",
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
    FairContainer* p1 = new FairContainer("TofdHitPar", "TOFd Hit Parameters", "TestDefaultContext");
    p1->addContext("TestNonDefaultContext");
    containers->Add(p1);
}

FairParSet* R3BTofContFact::createContainer(FairContainer* c)
{
    /** Tofls the constructor of the corresponding parameter container.
     * For an actual context, which is not an empty string and not the default context
     * of this container, the name is concatinated with the context. */

    const char* name = c->GetName();
    cout << " -I container name " << name << endl;
    FairParSet* p = 0;
    /*
    if (strcmp(name,"R3BTofDigiPar")==0) {
      p=new R3BTofDigiPar(c->getConcatName().Data(),c->GetTitle(),c->getContext());
    }
    if (strcmp(name,"CbmGeoStsPar")==0) {
      p=new CbmGeoStsPar(c->getConcatName().Data(),c->GetTitle(),c->getContext());
    }*/
    if (strcmp(name, "TofdHitPar") == 0)
    {
        p = new R3BTofdHitPar(c->getConcatName().Data(), c->GetTitle(), c->getContext());
    }
    return p;
}

void R3BTofContFact::activateParIo(FairParIo* io)
{
    // activates the input/output class for the parameters
    // needed by the Sts
    /*
    if (strcmp(io->IsA()->GetName(),"FairParRootFileIo")==0) {
      R3BTofParRootFileIo* p=new R3BTofParRootFileIo(((FairParRootFileIo*)io)->getParRootFile());
      io->setDetParIo(p);
    }
    if (strcmp(io->IsA()->GetName(),"FairParAsciiFileIo")==0) {
      R3BTofParAsciiFileIo* p=new R3BTofParAsciiFileIo(((FairParAsciiFileIo*)io)->getFile());
      io->setDetParIo(p);
      }
    */
}
