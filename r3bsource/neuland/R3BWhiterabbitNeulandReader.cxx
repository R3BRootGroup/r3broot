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

#include <FairRootManager.h>

#include "R3BEventHeader.h"
#include "R3BLogger.h"
#include "R3BWRData.h"
#include "R3BWhiterabbitNeulandReader.h"

#include <TClonesArray.h>

extern "C"
{
#include "ext_data_client.h"
#include "ext_h101_wrneuland.h"
}

R3BWhiterabbitNeulandReader::R3BWhiterabbitNeulandReader(EXT_STR_h101_WRNEULAND* data,
                                                         size_t offset,
                                                         UInt_t whiterabbit_id)
    : R3BReader("R3BWhiterabbitNeulandReader")
    , fNEvent(0)
    , fData(data)
    , fOffset(offset)
    , fOnline(kFALSE)
    , fWhiterabbitId(whiterabbit_id)
    , fEventHeader(nullptr)
    , fArray(new TClonesArray("R3BWRData"))
{
}

R3BWhiterabbitNeulandReader::~R3BWhiterabbitNeulandReader()
{
    if (fArray)
    {
        delete fArray;
    }
}

Bool_t R3BWhiterabbitNeulandReader::Init(ext_data_struct_info* a_struct_info)
{
    Int_t ok;
    R3BLOG(info, "");
    EXT_STR_h101_WRNEULAND_ITEMS_INFO(ok, *a_struct_info, fOffset, EXT_STR_h101_WRNEULAND, 0);
    if (!ok)
    {
        R3BLOG(error, "R3BWhiterabbitNeulandReader::Failed to setup structure information.");
        return kFALSE;
    }

    // Look for the R3BEventHeader
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
    FairRootManager::Instance()->Register("WRNeulandData", "WRNeuland", fArray, !fOnline);

    fData->NN_WR_ID = 0;

    return kTRUE;
}

Bool_t R3BWhiterabbitNeulandReader::R3BRead()
{
    if (!fData->NN_WR_ID)
    {
        return kTRUE;
    }

    if (fWhiterabbitId != fData->NN_WR_ID)
    {
        char strMessage[1000];
        snprintf(strMessage,
                 sizeof strMessage,
                 "Event %lu: Whiterabbit ID mismatch: expected 0x%x, got 0x%x.\n",
                 fEventHeader->GetEventno(),
                 fWhiterabbitId,
                 fData->NN_WR_ID);
        R3BLOG(error, strMessage);
    }

    if (fEventHeader != nullptr)
    {
        uint64_t timestamp = ((uint64_t)fData->NN_WR_T4 << 48) | ((uint64_t)fData->NN_WR_T3 << 32) |
                             ((uint64_t)fData->NN_WR_T2 << 16) | (uint64_t)fData->NN_WR_T1;
        // fEventHeader->SetTimeStamp(timestamp);
        fNEvent = fEventHeader->GetEventno();
        new ((*fArray)[fArray->GetEntriesFast()]) R3BWRData(timestamp);
    }
    else
    {
        fNEvent++;
    }

    fData->NN_WR_ID = 0;
    return kTRUE;
}

void R3BWhiterabbitNeulandReader::Reset()
{
    // Reset the output array
    fArray->Clear();
    fNEvent = 0;
}

ClassImp(R3BWhiterabbitNeulandReader);
