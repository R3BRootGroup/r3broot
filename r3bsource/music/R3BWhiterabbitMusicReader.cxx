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

#include "R3BWhiterabbitMusicReader.h"
#include "R3BEventHeader.h"
#include "R3BWRData.h"

#include <FairLogger.h>
#include <FairRootManager.h>

#include <TClonesArray.h>

extern "C"
{
#include "ext_data_client.h"
#include "ext_h101_wrmusic.h"
}

R3BWhiterabbitMusicReader::R3BWhiterabbitMusicReader(EXT_STR_h101_WRMUSIC* data, size_t offset, UInt_t whiterabbit_id)
    : R3BReader("R3BWhiterabbitMusicReader")
    , fNEvent(0)
    , fData(data)
    , fOffset(offset)
    , fOnline(kFALSE)
    , fWhiterabbitId(whiterabbit_id)
    , fEventHeader(nullptr)
    , fArray(new TClonesArray("R3BWRData"))
{
}

R3BWhiterabbitMusicReader::~R3BWhiterabbitMusicReader()
{
    if (fArray)
        delete fArray;
}

Bool_t R3BWhiterabbitMusicReader::Init(ext_data_struct_info* a_struct_info)
{
    Int_t ok;
    LOG(info) << "R3BWhiterabbitMusicReader::Init()";
    EXT_STR_h101_WRMUSIC_ITEMS_INFO(ok, *a_struct_info, fOffset, EXT_STR_h101_WRMUSIC, 0);

    if (!ok)
    {
        LOG(error) << "R3BWhiterabbitMusicReader::Failed to setup structure information.";
        return kFALSE;
    }

    // Look for the R3BEventHeader
    FairRootManager* frm = FairRootManager::Instance();
    fEventHeader = dynamic_cast<R3BEventHeader*>(frm->GetObject("EventHeader."));
    if (!fEventHeader)
    {
        LOG(warn) << "R3BWhiterabbitMusicReader::Init() EventHeader. not found";
    }
    else
        LOG(info) << "R3BWhiterabbitMusicReader::Init() R3BEventHeader found";

    // Register output array in tree
    FairRootManager::Instance()->Register("WRMusicData", "WRMusic", fArray, !fOnline);

    fData->TIMESTAMP_MUSIC_ID = 0;

    return kTRUE;
}

Bool_t R3BWhiterabbitMusicReader::R3BRead()
{
    if (!fData->TIMESTAMP_MUSIC_ID)
    {
        return kTRUE;
    }
    /*printf("ams %08x %08x %08x %08x %08x \n",
    fData->TIMESTAMP_MUSIC_ID,
    fData->TIMESTAMP_MUSIC_WR_T4, fData->TIMESTAMP_MUSIC_WR_T3,
    fData->TIMESTAMP_MUSIC_WR_T2, fData->TIMESTAMP_MUSIC_WR_T1);*/

    if (fWhiterabbitId != fData->TIMESTAMP_MUSIC_ID)
    {
        char strMessage[1000];
        snprintf(strMessage,
                 sizeof strMessage,
                 "Event %lu: Whiterabbit ID mismatch: expected 0x%x, got 0x%x.\n",
                 fEventHeader->GetEventno(),
                 fWhiterabbitId,
                 fData->TIMESTAMP_MUSIC_ID);
        LOG(error) << strMessage;
    }

    if (fEventHeader != nullptr)
    {
        uint64_t timestamp = ((uint64_t)fData->TIMESTAMP_MUSIC_WR_T4 << 48) |
                             ((uint64_t)fData->TIMESTAMP_MUSIC_WR_T3 << 32) |
                             ((uint64_t)fData->TIMESTAMP_MUSIC_WR_T2 << 16) | (uint64_t)fData->TIMESTAMP_MUSIC_WR_T1;
        // fEventHeader->SetTimeStamp(timestamp);
        fNEvent = fEventHeader->GetEventno();
        new ((*fArray)[fArray->GetEntriesFast()]) R3BWRData(timestamp);
    }
    else
    {
        fNEvent++;
    }

    fData->TIMESTAMP_MUSIC_ID = 0;
    return kTRUE;
}

void R3BWhiterabbitMusicReader::Reset()
{
    // Reset the output array
    fArray->Clear();
    fNEvent = 0;
}

ClassImp(R3BWhiterabbitMusicReader)
