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

#include "R3BNeulandReconstructionContFact.h"
#include "FairRuntimeDb.h"
#include "R3BNeulandMultiplicityBayesPar.h"
#include "R3BNeulandMultiplicityCalorimetricPar.h"

static R3BNeulandReconstructionContFact gR3BNeulandReconstructionContFact;

R3BNeulandReconstructionContFact::R3BNeulandReconstructionContFact()
{
    fName = "R3BNeulandReconstructionContFact";
    fTitle = "Factory for Neuland reconstruction parameter containers";
    setAllContainers();
    FairRuntimeDb::instance()->addContFactory(this);
}

void R3BNeulandReconstructionContFact::setAllContainers()
{
    auto bayes = new FairContainer(
        "R3BNeulandMultiplicityBayesPar", "Neuland Multiplicity Bayes Parameters", "TestDefaultContext");
    bayes->addContext("TestNonDefaultContext");
    containers->Add(bayes);

    auto calorimetric = new FairContainer(
        "R3BNeulandMultiplicityCalorimetricPar", "Neuland Multiplicity Calorimetric Parameters", "TestDefaultContext");
    calorimetric->addContext("TestNonDefaultContext");
    containers->Add(calorimetric);
}

FairParSet* R3BNeulandReconstructionContFact::createContainer(FairContainer* c)
{
    const char* name = c->GetName();

    if (strcmp(name, "R3BNeulandMultiplicityBayesPar") == 0)
    {
        return new R3BNeulandMultiplicityBayesPar(c->getConcatName().Data(), c->GetTitle(), c->getContext());
    }

    if (strcmp(name, "R3BNeulandMultiplicityCalorimetricPar") == 0)
    {
        return new R3BNeulandMultiplicityCalorimetricPar(c->getConcatName().Data(), c->GetTitle(), c->getContext());
    }

    return nullptr;
}

ClassImp(R3BNeulandReconstructionContFact)
