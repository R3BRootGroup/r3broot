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

#include <FairLogger.h>
#include <FairRootManager.h>

#include "R3BEventHeader.h"
#include "R3BWRData.h"
#include "R3BWhiterabbitCalifaReader.h"
#include <TClonesArray.h>

extern "C"
{
#include "ext_data_client.h"
#include "ext_h101_wrcalifa.h"
}

R3BWhiterabbitCalifaReader::R3BWhiterabbitCalifaReader(EXT_STR_h101_WRCALIFA* data,
                                                       size_t offset,
                                                       UInt_t whiterabbit_id1,
                                                       UInt_t whiterabbit_id2)
    : R3BReader("R3BWhiterabbitCalifaReader")
    , fNEvent(0)
    , fData(data)
    , fOffset(offset)
    , fOnline(kFALSE)
    , fWhiterabbitId1(whiterabbit_id1)
    , fWhiterabbitId2(whiterabbit_id2)
    , fEventHeader(nullptr)
    , fArray(new TClonesArray("R3BWRData"))
{
}

R3BWhiterabbitCalifaReader::~R3BWhiterabbitCalifaReader()
{
    if (fArray)
    {
        delete fArray;
    }
}

Bool_t R3BWhiterabbitCalifaReader::Init(ext_data_struct_info* a_struct_info)
{
    Int_t ok;
    LOG(info) << "R3BWhiterabbitCalifaReader::Init()";
    EXT_STR_h101_WRCALIFA_ITEMS_INFO(ok, *a_struct_info, fOffset, EXT_STR_h101_WRCALIFA, 0);

    if (!ok)
    {
        LOG(error) << "R3BWhiterabbitCalifaReader::Failed to setup structure information.";
        return kFALSE;
    }

    // Look for the R3BEventHeader
    FairRootManager* frm = FairRootManager::Instance();
    fEventHeader = dynamic_cast<R3BEventHeader*>(frm->GetObject("EventHeader."));
    if (!fEventHeader)
    {
        LOG(warn) << "R3BWhiterabbitCalifaReader::Init() EventHeader. not found";
        fEventHeader = dynamic_cast<R3BEventHeader*>(frm->GetObject("R3BEventHeader"));
    }
    else
        LOG(info) << "R3BWhiterabbitCalifaReader::Init() R3BEventHeader found";

    // Register output array in tree
    FairRootManager::Instance()->Register("WRCalifaData", "WRCalifa", fArray, !fOnline);

    fData->TIMESTAMP_CALIFA1ID = 0;
    fData->TIMESTAMP_CALIFA2ID = 0;
    return kTRUE;
}

Bool_t R3BWhiterabbitCalifaReader::R3BRead()
{
    if (!fData->TIMESTAMP_CALIFA1ID || !fData->TIMESTAMP_CALIFA2ID)
    {
        return kTRUE;
    }
    /*printf("califa %08x %08x %08x %08x %08x \n",
    fData->TIMESTAMP_CALIFA_ID,
    fData->TIMESTAMP_CALIFA_WR_T4, fData->TIMESTAMP_CALIFA_WR_T3,
    fData->TIMESTAMP_CALIFA_WR_T2, fData->TIMESTAMP_CALIFA_WR_T1);*/

    if (fWhiterabbitId1 != fData->TIMESTAMP_CALIFA1ID)
    {
        char strMessage[1000];
        snprintf(strMessage,
                 sizeof strMessage,
                 "Event %lu: Whiterabbit ID mismatch: expected 0x%x, got 0x%x.\n",
                 fEventHeader->GetEventno(),
                 fWhiterabbitId1,
                 fData->TIMESTAMP_CALIFA1ID);
        LOG(error) << strMessage;
    }

    if (fWhiterabbitId2 != fData->TIMESTAMP_CALIFA2ID)
    {
        char strMessage[1000];
        snprintf(strMessage,
                 sizeof strMessage,
                 "Event %lu: Whiterabbit ID mismatch: expected 0x%x, got 0x%x.\n",
                 fEventHeader->GetEventno(),
                 fWhiterabbitId2,
                 fData->TIMESTAMP_CALIFA2ID);
        LOG(error) << strMessage;
    }

    if (fEventHeader != nullptr)
    {
        uint64_t timestamp1 = ((uint64_t)fData->TIMESTAMP_CALIFA1WR_T4 << 48) |
                              ((uint64_t)fData->TIMESTAMP_CALIFA1WR_T3 << 32) |
                              ((uint64_t)fData->TIMESTAMP_CALIFA1WR_T2 << 16) | (uint64_t)fData->TIMESTAMP_CALIFA1WR_T1;

        fNEvent = fEventHeader->GetEventno();
        new ((*fArray)[fArray->GetEntriesFast()]) R3BWRData(timestamp1);

        uint64_t timestamp2 = ((uint64_t)fData->TIMESTAMP_CALIFA2WR_T4 << 48) |
                              ((uint64_t)fData->TIMESTAMP_CALIFA2WR_T3 << 32) |
                              ((uint64_t)fData->TIMESTAMP_CALIFA2WR_T2 << 16) | (uint64_t)fData->TIMESTAMP_CALIFA2WR_T1;
        new ((*fArray)[fArray->GetEntriesFast()]) R3BWRData(timestamp2);
    }
    else
    {
        fNEvent++;
    }

    fData->TIMESTAMP_CALIFA1ID = 0;
    fData->TIMESTAMP_CALIFA2ID = 0;
    return kTRUE;
}

void R3BWhiterabbitCalifaReader::Reset()
{
    // Reset the output array
    fArray->Clear();
    fNEvent = 0;
}

ClassImp(R3BWhiterabbitCalifaReader)
