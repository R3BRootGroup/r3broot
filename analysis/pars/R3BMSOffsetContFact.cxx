/******************************************************************************
 *   Copyright (C) 2019 GSI Helmholtzzentrum für Schwerionenforschung GmbH    *
 *   Copyright (C) 2019-2024 Members of R3B Collaboration                     *
 *                                                                            *
 *             This software is distributed under the terms of the            *
 *                 GNU General Public Licence (GPL) version 3,                *
 *                    copied verbatim in the file "LICENSE".                  *
 *                                                                            *
 * In applying this license GSI does not waive the privileges and immunities  *
 * granted to it by virtue of its status as an Intergovernmental Organization *
 * or submit itself to any jurisdiction.                                      *
 ******************************************************************************/

//
//  R3BMSOffsetContFact
//
//  Factory for the parameter containers in libR3BMSOffset
//

#include "R3BMSOffsetContFact.h"
#include "R3BMSOffsetPar.h"

#include "FairParAsciiFileIo.h"
#include "FairParRootFileIo.h"
#include "FairRuntimeDb.h"

#include "R3BLogger.h"

#include "TClass.h"

static const R3BMSOffsetContFact gR3BMSOffsetContFact;

R3BMSOffsetContFact::R3BMSOffsetContFact()
{
    // Constructor (called when the library is loaded)
    fName = "R3BMSOffsetContFact";
    fTitle = "Factory for parameter containers in libR3BMSOffset";
    setAllContainers();
    FairRuntimeDb::instance()->addContFactory(this);
}

void R3BMSOffsetContFact::setAllContainers()
{
    // Creates the Container objects with all accepted contexts and adds them to
    // the list of containers for the STS library.

    auto pOffset = std::make_unique<FairContainer>("MSOffsetPar", "MSOffsetPar Offset Parameters", "MSOffsetPar");
    pOffset->addContext("MSOffsetPar");

    containers->Add(pOffset.release());
}

FairParSet* R3BMSOffsetContFact::createContainer(FairContainer* cOffset)
{
    // Trals the constructor of the corresponding parameter container.
    // For an actual context, which is not an empty string and not the default context
    // of this container, the name is concatinated with the context.

    const auto name = std::string_view{ cOffset->GetName() };
    R3BLOG(info, name);
    FairParSet* pOff = nullptr;
    if (name == "MSOffsetPar")
    {
        pOff = new R3BMSOffsetPar(cOffset->getConcatName().Data(), cOffset->GetTitle(), cOffset->getContext());
    }
    return pOff;
}

ClassImp(R3BMSOffsetContFact);
