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

#include "R3BSyncCheckReader.h"
#include "R3BEventHeader.h"
#include "R3BLogger.h"
#include "R3BSyncCheckData.h"

#include <FairLogger.h>
#include <FairRootManager.h>

#include <TClonesArray.h>

extern "C"
{
#include "ext_data_client.h"
#include "ext_h101_sync_check.h"
}

R3BSyncCheckReader::R3BSyncCheckReader(EXT_STR_h101_SYNC_CHECK* data, size_t offset)
    : R3BReader("R3BSyncCheckReader")
    , fNEvent(1)
    , fData(data)
    , fOffset(offset)
    , fArray(new TClonesArray("R3BSyncCheckData"))
    , fStoreData(kFALSE)
{
}

R3BSyncCheckReader::~R3BSyncCheckReader()
{
    if (fArray)
    {
        delete fArray;
    }
}

Bool_t R3BSyncCheckReader::Init(ext_data_struct_info* a_struct_info)
{
    Int_t ok;
    R3BLOG(info, "");
    EXT_STR_h101_SYNC_CHECK_ITEMS_INFO(ok, *a_struct_info, fOffset, EXT_STR_h101_SYNC_CHECK, 0);
    if (!ok)
    {
        R3BLOG(fatal, "Failed to setup structure information");
        return kFALSE;
    }
    // Register output array in tree
    FairRootManager::Instance()->Register("SyncCheckData", "SyncCheck", fArray, fStoreData);
    Reset();
    memset(fData, 0, sizeof *fData);
    return kTRUE;
}

Bool_t R3BSyncCheckReader::R3BRead()
{
    std::vector<uint32_t> foot_input;
    foot_input.push_back(fData->SYNC_CHECK_FT1V);
    foot_input.push_back(fData->SYNC_CHECK_FT2V);
    foot_input.push_back(fData->SYNC_CHECK_FT3V);
    foot_input.push_back(fData->SYNC_CHECK_FT4V);
    foot_input.push_back(fData->SYNC_CHECK_FT5V);
    foot_input.push_back(fData->SYNC_CHECK_FT6V);
    foot_input.push_back(fData->SYNC_CHECK_FT7V);
    foot_input.push_back(fData->SYNC_CHECK_FT8V);
    foot_input.push_back(fData->SYNC_CHECK_FT9V);
    foot_input.push_back(fData->SYNC_CHECK_FT10V);
    foot_input.push_back(fData->SYNC_CHECK_FT11V);
    foot_input.push_back(fData->SYNC_CHECK_FT12V);
    foot_input.push_back(fData->SYNC_CHECK_FT13V);
    foot_input.push_back(fData->SYNC_CHECK_FT14V);
    foot_input.push_back(fData->SYNC_CHECK_FT15V);
    foot_input.push_back(fData->SYNC_CHECK_FT16V);

    new ((*fArray)[fArray->GetEntriesFast()]) R3BSyncCheckData(fData->SYNC_CHECK_MASTER,
                                                               fData->SYNC_CHECK_MASTERRR,
                                                               fData->SYNC_CHECK_MUSIC,
                                                               fData->SYNC_CHECK_RPC,
                                                               fData->SYNC_CHECK_STWO,
                                                               foot_input);
    fNEvent++;
    foot_input.clear();
    return kTRUE;
}

void R3BSyncCheckReader::Reset() { fArray->Clear(); }

ClassImp(R3BSyncCheckReader)
