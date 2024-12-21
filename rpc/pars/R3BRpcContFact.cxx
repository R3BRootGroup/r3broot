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

#include "R3BRpcContFact.h"

#include "FairLogger.h"
#include "FairParAsciiFileIo.h"
#include "FairParRootFileIo.h"
#include "FairRuntimeDb.h"

#include "R3BRpcHitPar.h"
#include "R3BRpcPars4Sim.h"
#include "R3BRpcTotCalPar.h"
#include "R3BTGeoPar.h"

#include "TClass.h"

static R3BRpcContFact gR3BRpcContFact;

R3BRpcContFact::R3BRpcContFact()
{
    // Constructor (called when the library is loaded)
    fName = "R3BRpcContFact";
    fTitle = "Factory for parameter containers in libR3BRpc";
    setAllContainers();
    FairRuntimeDb::instance()->addContFactory(this);
}

void R3BRpcContFact::setAllContainers()
{
    // Creates the Container objects with all accepted contexts and adds them
    // to the list of containers for the RPC library.

    FairContainer* p1 = new FairContainer("RpcTotCalPar", "RPC Calibration Parameters", "RPCCalParContext");
    p1->addContext("RPCCalParContext");
    containers->Add(p1);

    FairContainer* p2 = new FairContainer("RpcHitPar", "RPC Hit Parameters", "RPCHitParContext");
    p2->addContext("RPCHitParContext");
    containers->Add(p2);

    FairContainer* p3 = new FairContainer("RpcPars4Sim", "RPC Parameters for Sim", "RPCSimParContext");
    p3->addContext("RPCSimParContext");
    containers->Add(p3);

    FairContainer* p4 = new FairContainer("RpcGeoPar", "RPC geometry parameters", "GeometryParameterContext");
    p4->addContext("GeometryParameterContext");
    containers->Add(p4);
}

FairParSet* R3BRpcContFact::createContainer(FairContainer* c)
{
    /** Calls the constructor of the corresponding parameter container.
     * For an actual context, which is not an empty string and not the default context
     * of this container, the name is concatinated with the context.
     */

    const char* name = c->GetName();
    LOG(info) << "R3BRpcContFact: Create container name: " << name;
    FairParSet* p = 0;
    if (strcmp(name, "RpcTotCalPar") == 0)
    {
        p = new R3BRpcTotCalPar(c->getConcatName().Data(), c->GetTitle(), c->getContext());
    }
    else if (strcmp(name, "RpcHitPar") == 0)
    {
        p = new R3BRpcHitPar(c->getConcatName().Data(), c->GetTitle(), c->getContext());
    }
    else if (strcmp(name, "RpcPars4Sim") == 0)
    {
        p = new R3BRpcPars4Sim(c->getConcatName().Data(), c->GetTitle(), c->getContext());
    }
    else if (strcmp(name, "RpcGeoPar") == 0)
    {
        p = new R3BTGeoPar(c->getConcatName().Data(), c->GetTitle(), c->getContext());
    }
    return p;
}

ClassImp(R3BRpcContFact);
