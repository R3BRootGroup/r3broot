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

#include "R3BWhiterabbitReader.h"
#include "FairLogger.h"
#include "FairRootManager.h"
#include "R3BEventHeader.h"

extern "C"
{
#include "ext_data_client.h"
#include "ext_h101_whiterabbit.h"
}

R3BWhiterabbitReader::R3BWhiterabbitReader(EXT_STR_h101_whiterabbit* data, size_t offset, UInt_t whiterabbit_id)
    : R3BReader("R3BWhiterabbitReader")
    , fNEvent(0)
    , fData(data)
    , fOffset(offset)
    , fWhiterabbitId(whiterabbit_id)
    , fEventHeader(nullptr)
{
}

R3BWhiterabbitReader::~R3BWhiterabbitReader() {}

Bool_t R3BWhiterabbitReader::Init(ext_data_struct_info* a_struct_info)
{
    Int_t ok;
    LOG(info) << "R3BWhiterabbitReader::Init()";
    EXT_STR_h101_whiterabbit_ITEMS_INFO(ok, *a_struct_info, fOffset, EXT_STR_h101_whiterabbit, 0);

    if (!ok)
    {
        perror("ext_data_struct_info_item");
        LOG(error) << "Failed to setup structure information.";
        return kFALSE;
    }

    // Look for the R3BEventHeader
    FairRootManager* frm = FairRootManager::Instance();
    fEventHeader = (R3BEventHeader*)frm->GetObject("EventHeader.");
    if (!fEventHeader)
    {
        LOG(warn) << "R3BWhiterabbitReader::Init() R3BEventHeader not found";
        fEventHeader = (R3BEventHeader*)frm->GetObject("R3BEventHeader");
    }
    else
        LOG(info) << "R3BWhiterabbitReader::Init() R3BEventHeader found";

    fData->TIMESTAMP_MASTER_ID = 0;

    return kTRUE;
}

Bool_t R3BWhiterabbitReader::R3BRead()
{
    if (!fData->TIMESTAMP_MASTER_ID)
    {
        return kTRUE;
    }
    printf("%08x %08x %08x %08x %08x \n",
           fData->TIMESTAMP_MASTER_ID,
           fData->TIMESTAMP_MASTER_WR_T4,
           fData->TIMESTAMP_MASTER_WR_T3,
           fData->TIMESTAMP_MASTER_WR_T2,
           fData->TIMESTAMP_MASTER_WR_T1);
    if (fWhiterabbitId != fData->TIMESTAMP_MASTER_ID)
    {
        char strMessage[1000];
        snprintf(strMessage,
                 sizeof strMessage,
                 "Event %lu: Whiterabbit ID mismatch: expected 0x%x, got 0x%x.\n",
                 fEventHeader->GetEventno(),
                 fWhiterabbitId,
                 fData->TIMESTAMP_MASTER_ID);
        LOG(error) << strMessage;
    }

    if (fEventHeader != nullptr)
    {
        uint64_t timestamp = ((uint64_t)fData->TIMESTAMP_MASTER_WR_T4 << 48) |
                             ((uint64_t)fData->TIMESTAMP_MASTER_WR_T3 << 32) |
                             ((uint64_t)fData->TIMESTAMP_MASTER_WR_T2 << 16) | (uint64_t)fData->TIMESTAMP_MASTER_WR_T1;

        fEventHeader->SetTimeStamp(timestamp);
        fNEvent = fEventHeader->GetEventno();
    }
    else
    {
        fNEvent++;
    }

    fData->TIMESTAMP_MASTER_ID = 0;
    return kTRUE;
}

void R3BWhiterabbitReader::Reset() { fNEvent = 0; }

ClassImp(R3BWhiterabbitReader)
