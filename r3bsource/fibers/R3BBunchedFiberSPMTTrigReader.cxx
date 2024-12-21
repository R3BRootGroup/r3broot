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

#include "R3BBunchedFiberSPMTTrigReader.h"
#include "FairLogger.h"
#include "FairRootManager.h"
#include "R3BFiberMappedData.h"
#include "R3BLogger.h"
#include "TClonesArray.h"

extern "C"
{
#include "ext_data_client.h"
#include "ext_h101_fib.h"
}

R3BBunchedFiberSPMTTrigReader::R3BBunchedFiberSPMTTrigReader(EXT_STR_h101_FIB* a_data, size_t a_offset)
    : R3BReader("R3BBunchedFiberSPMTTrigReader")
    , fData(reinterpret_cast<EXT_STR_h101_FIB_onion*>(a_data))
    , fOffset(a_offset)
    , fOnline(kFALSE)
    , fMappedArray(new TClonesArray("R3BFiberMappedData"))
{
}

R3BBunchedFiberSPMTTrigReader::~R3BBunchedFiberSPMTTrigReader()
{
    if (fMappedArray)
        delete fMappedArray;
}

Bool_t R3BBunchedFiberSPMTTrigReader::Init(ext_data_struct_info* a_struct_info)
{
    int ok;
    R3BLOG(info, "");
    EXT_STR_h101_FIB_ITEMS_INFO(ok, *a_struct_info, fOffset, EXT_STR_h101_FIB, 0);
    if (!ok)
    {
        R3BLOG(fatal, "Failed to setup structure information.");
        return kFALSE;
    }

    FairRootManager::Instance()->Register(
        "BunchedFiberSPMTTrigMapped", "BunchedFiberSPMTTrigMapped data", fMappedArray, !fOnline);
    memset(fData, 0, sizeof *fData);

    return kTRUE;
}

Bool_t R3BBunchedFiberSPMTTrigReader::R3BRead()
{
    auto* data = reinterpret_cast<EXT_STR_h101_FIB_onion*>(fData);

    // Leading TAMEX trigger times.
    auto numChannels = data->FIB_TRIGSLF;
    for (uint32_t i = 0; i < numChannels; i++)
    {
        uint32_t channel = data->FIB_TRIGSLFI[i];
        new ((*fMappedArray)[fMappedArray->GetEntriesFast()])
            R3BFiberMappedData(4, channel, true, data->FIB_TRIGSLCv[i], data->FIB_TRIGSLFv[i]);
    }

    return kTRUE;
}

void R3BBunchedFiberSPMTTrigReader::Reset()
{
    // Reset the output array
    fMappedArray->Clear();
}

ClassImp(R3BBunchedFiberSPMTTrigReader);
