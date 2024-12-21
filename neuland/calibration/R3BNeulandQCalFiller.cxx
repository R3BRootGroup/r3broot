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

#include "R3BNeulandQCalFiller.h"
#include "R3BEventHeader.h"
#include "R3BNeulandTacquilaMappedData.h"
#include "TClonesArray.h"
#include <FairRootManager.h>

R3BNeulandQCalFiller::R3BNeulandQCalFiller()
{
    data = new TClonesArray("R3BNeulandTacquilaMappedData", 100);
    header = new R3BEventHeader();
}

R3BNeulandQCalFiller::~R3BNeulandQCalFiller()
{
    delete data;
    delete header;
}

InitStatus R3BNeulandQCalFiller::Init()
{
    FairRootManager::Instance()->Register("NeulandTacquilaMappedData", "Neuland", data, kTRUE);
    FairRootManager::Instance()->Register("R3BEventHeader", "", header, kTRUE);
    return kSUCCESS;
}

void R3BNeulandQCalFiller::Exec(Option_t* option)
{
    Int_t value;
    i++;
    for (Int_t j = 0; j < 100; j++)
    {
        value = j + 100;
        if (i > 5000)
            value += std::rand() % 10 - 5;
        new ((*data)[j]) R3BNeulandTacquilaMappedData(0, 0, 0, value, 1, j / 2 + 1, j % 2 + 1);
    }
    header = new (header) R3BEventHeader();
    header->SetTrigger(2);
}

ClassImp(R3BNeulandQCalFiller)
