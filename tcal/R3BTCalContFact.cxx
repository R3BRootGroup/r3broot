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

#include "FairParSet.h"
#include "FairRuntimeDb.h"

#include "R3BLogger.h"
#include "R3BTCalContFact.h"
#include "R3BTCalPar.h"

#include "TList.h"
#include <string.h>

static R3BTCalContFact gR3BTCalContFact;

R3BTCalContFact::R3BTCalContFact()
{
    // Constructor (called when the library is loaded)
    fName = "R3BTCalContFact";
    fTitle = "Tutorial factory for parameter containers";
    setAllContainers();
    FairRuntimeDb::instance()->addContFactory(this);
}

void R3BTCalContFact::addContainer(TString name, TString description)
{
    auto container = new FairContainer(name, description, "TestDefaultContext");
    container->addContext("TestNonDefaultContext");
    containers->Add(container);
}

void R3BTCalContFact::setAllContainers()
{
    /** Creates the Container objects with all accepted contexts and adds
     *  them to the list of containers.*/
    addContainer("LandTCalPar", "NeuLAND TCAL Calibration Parameters");
    addContainer("LosTCalPar", "LOS TCAL Calibration Parameters");
    addContainer("RoluTCalPar", "ROLU TCAL Calibration Parameters");
    addContainer("Sci2TCalPar", "SCI2 TCAL Calibration Parameters");
    addContainer("Sci8TCalPar", "SCI8 TCAL Calibration Parameters");
    addContainer("TofdTCalPar", "TOFD TCAL Calibration Parameters");
    addContainer("PdcTCalPar", "PDC TCAL Calibration Parameters");
    addContainer("StrawtubesTCalPar", "Strawtubes TCAL Calibration Parameters");
    addContainer("BunchedFiberSPMTTrigTCalPar", "Bunched Fiber SPMT Trigger TCAL Calibration Parameters");
    addContainer("SfibTCalPar", "Sfib TCAL Calibration Parameters");
    addContainer("PtofTCalPar", "PTOF TCAL Calibration Parameters");
    addContainer("Sci2TCalPar", "S2 TCAL Calibration Parameters");
    addContainer("Sci8TCalPar", "S8 TCAL Calibration Parameters");
    addContainer("RpcTCalPar", "Rpc TCAL Calibration Parameters");
    ADD_FIBER(Fi0, FI0);
    ADD_FIBER(Fi1a, FI1a);
    ADD_FIBER(Fi1b, FI1b);
    ADD_FIBER(Fi2a, FI2a);
    ADD_FIBER(Fi2b, FI2b);
    ADD_FIBER(Fi3a, FI3a);
    ADD_FIBER(Fi3b, FI3b);
    ADD_FIBER(Fi4, FI4);
    ADD_FIBER(Fi5, FI5);
    ADD_FIBER(Fi6, FI6);
    ADD_FIBER(Fi7, FI7);
    ADD_FIBER(Fi8, FI8);
    ADD_FIBER(Fi9, FI9);
    ADD_FIBER(Fi10, FI10);
    ADD_FIBER(Fi11, FI11);
    ADD_FIBER(Fi12, FI12);
    ADD_FIBER(Fi13, FI13);
    ADD_FIBER(Fi30, FI30);
    ADD_FIBER(Fi31, FI31);
    ADD_FIBER(Fi32, FI32);
    ADD_FIBER(Fi33, FI33);
}

FairParSet* R3BTCalContFact::createContainer(FairContainer* c)
{
    /** Calls the constructor of the corresponding parameter container.
     * For an actual context, which is not an empty string and not the default context
     * of this container, the name is concatinated with the context. */

    const char* name = c->GetName();
    R3BLOG(info, name);

    containerNames.push_back("LandTCalPar");
    containerNames.push_back("LosTCalPar");
    containerNames.push_back("RoluTCalPar");
    containerNames.push_back("Sci2TCalPar");
    containerNames.push_back("Sci8TCalPar");
    containerNames.push_back("TofdTCalPar");
    containerNames.push_back("PdcTCalPar");
    containerNames.push_back("StrawtubesTCalPar");
    containerNames.push_back("BunchedFiberSPMTTrigTCalPar");
    containerNames.push_back("SfibTCalPar");
    containerNames.push_back("PtofTCalPar");
    containerNames.push_back("Sci2TCalPar");
    containerNames.push_back("Sci8TCalPar");
    containerNames.push_back("RpcTCalPar");
    PUSH_FIBER(Fi0);
    PUSH_FIBER(Fi1a);
    PUSH_FIBER(Fi1b);
    PUSH_FIBER(Fi2a);
    PUSH_FIBER(Fi2b);
    PUSH_FIBER(Fi3a);
    PUSH_FIBER(Fi3b);
    PUSH_FIBER(Fi4);
    PUSH_FIBER(Fi5);
    PUSH_FIBER(Fi6);
    PUSH_FIBER(Fi7);
    PUSH_FIBER(Fi8);
    PUSH_FIBER(Fi9);
    PUSH_FIBER(Fi10);
    PUSH_FIBER(Fi11);
    PUSH_FIBER(Fi12);
    PUSH_FIBER(Fi13);
    PUSH_FIBER(Fi30);
    PUSH_FIBER(Fi31);
    PUSH_FIBER(Fi32);
    PUSH_FIBER(Fi33);

    // bool found = false;
    for (auto containerName : containerNames)
    {
        if (strncmp(name, containerName, strlen(containerName)) == 0)
        {
            // found = true;
            // break;
            return new R3BTCalPar(c->getConcatName().Data(), c->GetTitle(), c->getContext());
        }
    }

    /* if (found == true)
     {
         return new R3BTCalPar(c->getConcatName().Data(), c->GetTitle(), c->getContext());
     }
     else
     {*/
    return nullptr;
    // }
}

ClassImp(R3BTCalContFact);
