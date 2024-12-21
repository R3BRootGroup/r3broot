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

#include "FairLogger.h"
#include "FairRootManager.h"

#include "R3BEventHeader.h"
#include "R3BLogger.h"
#include "R3BWRData.h"
#include "R3BWhiterabbitMasterReader.h"

#include "TClonesArray.h"

extern "C"
{
#include "ext_data_client.h"
#include "ext_h101_wrmaster.h"
}

R3BWhiterabbitMasterReader::R3BWhiterabbitMasterReader(EXT_STR_h101_WRMASTER* data,
                                                       size_t offset,
                                                       UInt_t whiterabbit_id)
    : R3BReader("R3BWhiterabbitMasterReader")
    , fNEvent(1)
    , fData(data)
    , fOffset(offset)
    , fOnline(kTRUE)
    , fWhiterabbitId(whiterabbit_id)
    , fEventHeader(nullptr)
    , fArray(new TClonesArray("R3BWRData"))
{
}

R3BWhiterabbitMasterReader::~R3BWhiterabbitMasterReader()
{
    if (fArray)
    {
        delete fArray;
    }
}

Bool_t R3BWhiterabbitMasterReader::Init(ext_data_struct_info* a_struct_info)
{
    Int_t ok;
    R3BLOG(info, "");
    EXT_STR_h101_WRMASTER_ITEMS_INFO(ok, *a_struct_info, fOffset, EXT_STR_h101_WRMASTER, 0);
    if (!ok)
    {
        R3BLOG(fatal, "Failed to setup structure information");
        return kFALSE;
    }

    // Looking for the R3BEventHeader
    FairRootManager* frm = FairRootManager::Instance();
    fEventHeader = dynamic_cast<R3BEventHeader*>(frm->GetObject("EventHeader."));
    if (!fEventHeader)
    {
        R3BLOG(warn, "EventHeader. not found");
        fEventHeader = dynamic_cast<R3BEventHeader*>(frm->GetObject("R3BEventHeader"));
    }
    else
        R3BLOG(info, "EventHeader. found");

    // Register output array in tree
    FairRootManager::Instance()->Register("WRMasterData", "WRMaster", fArray, !fOnline);
    Reset();
    memset(fData, 0, sizeof *fData);

    return kTRUE;
}

Bool_t R3BWhiterabbitMasterReader::R3BRead()
{
    if (!fData->TIMESTAMP_MASTER_ID)
    {
        return kTRUE;
    }

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
        /*printf("master %08x %08x %08x %08x %08x \n",
        fData->TIMESTAMP_MASTER_ID,
        fData->TIMESTAMP_MASTER_WR_T4, fData->TIMESTAMP_MASTER_WR_T3,
        fData->TIMESTAMP_MASTER_WR_T2, fData->TIMESTAMP_MASTER_WR_T1);*/
        uint64_t timestamp = ((uint64_t)fData->TIMESTAMP_MASTER_WR_T4 << 48) |
                             ((uint64_t)fData->TIMESTAMP_MASTER_WR_T3 << 32) |
                             ((uint64_t)fData->TIMESTAMP_MASTER_WR_T2 << 16) | (uint64_t)fData->TIMESTAMP_MASTER_WR_T1;
        fNEvent = fEventHeader->GetEventno();
        fEventHeader->SetTimeStamp(timestamp);
        new ((*fArray)[fArray->GetEntriesFast()]) R3BWRData(timestamp);
    }
    else
    {
        fNEvent++;
    }

    fData->TIMESTAMP_MASTER_ID = 0;
    return kTRUE;
}

void R3BWhiterabbitMasterReader::Reset()
{
    // Reset the output array
    fArray->Clear();
}

ClassImp(R3BWhiterabbitMasterReader);
