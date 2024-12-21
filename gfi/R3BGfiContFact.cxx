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
//  R3BGfiContFact
//
//  Factory for the parameter containers in libR3BGfi
//
/////////////////////////////////////////////////////////////

#include "R3BGfiContFact.h"

// #include "R3BGfiParRootFileIo.h"
// #include "R3BGfiParAsciiFileIo.h"
#include "R3BGeoGfiPar.h"
// #include "R3BGfiDigiPar.h"

#include "FairRuntimeDb.h"
// #include "CbmParTest.h"
#include "FairParAsciiFileIo.h"
#include "FairParRootFileIo.h"

#include "TClass.h"

#include <iomanip>
#include <iostream>

using std::cout;
using std::endl;

ClassImp(R3BGfiContFact)

    static R3BGfiContFact gR3BGfiContFact;

R3BGfiContFact::R3BGfiContFact()
{
    // Constructor (called when the library is loaded)
    fName = "R3BGfiContFact";
    fTitle = "Factory for parameter containers in libR3BGfi";
    setAllContainers();
    FairRuntimeDb::instance()->addContFactory(this);
}

void R3BGfiContFact::setAllContainers()
{
    /** Creates the Container objects with all accepted contexts and adds them to
     *  the list of containers for the STS library.*/
    /*
     FairContainer* p1= new FairContainer("R3BGfiDigiPar",
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
}

FairParSet* R3BGfiContFact::createContainer(FairContainer* c)
{
    /** Gfils the constructor of the corresponding parameter container.
     * For an actual context, which is not an empty string and not the default context
     * of this container, the name is concatinated with the context. */
    /*
    const char* name=c->GetName();
    cout << " -I container name " << name << endl;
    FairParSet* p=0;
    if (strcmp(name,"R3BGfiDigiPar")==0) {
      p=new R3BGfiDigiPar(c->getConcatName().Data(),c->GetTitle(),c->getContext());
    }
    if (strcmp(name,"CbmGeoStsPar")==0) {
      p=new CbmGeoStsPar(c->getConcatName().Data(),c->GetTitle(),c->getContext());
    }
    return p;
    */
    return NULL;
}

void R3BGfiContFact::activateParIo(FairParIo* io)
{
    // activates the input/output class for the parameters
    // needed by the Sts
    /*
    if (strcmp(io->IsA()->GetName(),"FairParRootFileIo")==0) {
      R3BGfiParRootFileIo* p=new R3BGfiParRootFileIo(((FairParRootFileIo*)io)->getParRootFile());
      io->setDetParIo(p);
    }
    if (strcmp(io->IsA()->GetName(),"FairParAsciiFileIo")==0) {
      R3BGfiParAsciiFileIo* p=new R3BGfiParAsciiFileIo(((FairParAsciiFileIo*)io)->getFile());
      io->setDetParIo(p);
      }
    */
}
