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

#include "R3BRpcReader.h"
#include "R3BLogger.h"
#include "R3BRpcMappedData.h"

#include <FairRootManager.h>

#include <TClonesArray.h>

extern "C"
{
#include "ext_data_client.h"
#include "ext_h101_rpc.h"
}

R3BRpcReader::R3BRpcReader(EXT_STR_h101_RPC* data, size_t offset)
    : R3BReader("R3BRpcReader")
    , fData(data)
    , fArrayRpc(new TClonesArray("R3BRpcMappedData"))
    , fOnline(kFALSE)
    , fOffset(offset)
    , fCoarseReference(0)
    , fNEvent(0)
{
}

R3BRpcReader::~R3BRpcReader()
{
    if (fArrayRpc)
        delete fArrayRpc;
}

Bool_t R3BRpcReader::Init(ext_data_struct_info* a_struct_info)
{
    Int_t ok;
    R3BLOG(info, "");
    EXT_STR_h101_RPC_ITEMS_INFO(ok, *a_struct_info, fOffset, EXT_STR_h101_RPC, 0);
    if (!ok)
    {
        R3BLOG(error, "Failed to setup structure information.");
        return kFALSE;
    }

    // Register output array in tree
    FairRootManager::Instance()->Register("R3BRpcMappedData", "Rpc Mapped Data", fArrayRpc, !fOnline);
    Reset();
    memset(fData, 0, sizeof *fData);

    return kTRUE;
}

Bool_t R3BRpcReader::R3BRead()
{
    R3BLOG(debug1, "Event data.");

    uint64_t Fine_time;
    uint64_t Coarse_time;
    UShort_t DetId;
    UShort_t channelNumber;
    UShort_t Edge;
    UShort_t Side;

    for (int i = 0; i < fData->RPC_FT_RM; i++)
    {
        int begin = 0;

        if (i != 0)
        {
            begin = fData->RPC_FT_RME[i - 1];
        }

        for (int t = begin; t < fData->RPC_FT_RME[i]; t++)
        {
            DetId = 0;
            channelNumber = fData->RPC_FT_RMI[i];
            Fine_time = fData->RPC_FT_Rv[t];
            Coarse_time = fData->RPC_CT_Rv[t];
            Edge = fData->RPC_E_Rv[t];
            Side = 1;

            new ((*fArrayRpc)[fArrayRpc->GetEntriesFast()])
                R3BRpcMappedData(DetId, channelNumber, Fine_time, Coarse_time, Edge, Side);

            // R3BRpcMappedData(UShort_t DetId, UShort_t StripId, uint64_t FineTime, uint64_t CoarseTime, UShort_t Edge,
            // UShort_t Side)
        }
    }

    for (int i = 0; i < fData->RPC_FT_LM; i++)
    {
        int begin = 0;

        if (i != 0)
        {
            begin = fData->RPC_FT_LME[i - 1];
        }

        for (int t = begin; t < fData->RPC_FT_LME[i]; t++)
        {

            DetId = 0;
            channelNumber = fData->RPC_FT_RMI[i];
            channelNumber = fData->RPC_FT_LMI[i];
            Fine_time = fData->RPC_FT_Lv[t];
            Coarse_time = fData->RPC_CT_Lv[t];
            Edge = fData->RPC_E_Lv[t];
            Side = 0;

            new ((*fArrayRpc)[fArrayRpc->GetEntriesFast()])
                R3BRpcMappedData(DetId, channelNumber, Fine_time, Coarse_time, Edge, Side);

            // R3BRpcMappedData(UShort_t DetId, UShort_t StripId, uint64_t FineTime, uint64_t CoarseTime, UShort_t Edge,
            // UShort_t Side)
        }
    }

    // -------------------- Filling REF data ----------------------
    for (int i = 0; i < fData->RPC_REF_FTM; i++)
    {
        int begin = 0;

        if (i != 0)
        {
            begin = fData->RPC_REF_FTME[i - 1];
        }

        for (int t = begin; t < fData->RPC_REF_FTME[i]; t++)
        {

            DetId = 2;
            channelNumber = fData->RPC_REF_FTMI[i];
            Fine_time = fData->RPC_REF_FTv[t];
            Coarse_time = fData->RPC_REF_CTv[t];

            new ((*fArrayRpc)[fArrayRpc->GetEntriesFast()])
                R3BRpcMappedData(DetId, channelNumber, Fine_time, Coarse_time, 0, 0);

            // R3BRpcMappedData(UShort_t DetId, UShort_t StripId, uint64_t FineTime, uint64_t CoarseTime, UShort_t Edge,
            // UShort_t Side)
        }
    }

    // -------------------- Filling PMT data ----------------------
    for (int i = 0; i < fData->RPC_FT_B_PMTM; i++)
    {
        int begin = 0;

        if (i != 0)
        {
            begin = fData->RPC_FT_B_PMTME[i - 1];
        }

        for (int t = begin; t < fData->RPC_FT_B_PMTME[i]; t++)
        {

            DetId = 1;
            channelNumber = fData->RPC_FT_B_PMTMI[i];
            Fine_time = fData->RPC_FT_B_PMTv[t];
            Coarse_time = fData->RPC_CT_B_PMTv[t];

            Edge = fData->RPC_E_B_PMTv[t];
            Side = 1;

            new ((*fArrayRpc)[fArrayRpc->GetEntriesFast()])
                R3BRpcMappedData(DetId, channelNumber, Fine_time, Coarse_time, Edge, Side);

            // R3BRpcMappedData(UShort_t DetId, UShort_t StripId, uint64_t FineTime, uint64_t CoarseTime, UShort_t Edge,
            // UShort_t Side)
        }
    }

    for (int i = 0; i < fData->RPC_FT_T_PMTM; i++)
    {
        int begin = 0;

        if (i != 0)
        {
            begin = fData->RPC_FT_T_PMTME[i - 1];
        }

        for (int t = begin; t < fData->RPC_FT_T_PMTME[i]; t++)
        {

            DetId = 1;
            channelNumber = fData->RPC_FT_T_PMTMI[i];
            Fine_time = fData->RPC_FT_T_PMTv[t];
            Coarse_time = fData->RPC_CT_T_PMTv[t];

            Edge = fData->RPC_E_T_PMTv[t];
            Side = 0;

            new ((*fArrayRpc)[fArrayRpc->GetEntriesFast()])
                R3BRpcMappedData(DetId, channelNumber, Fine_time, Coarse_time, Edge, Side);

            // R3BRpcMappedData(UShort_t DetId, UShort_t StripId, uint64_t FineTime, uint64_t CoarseTime, UShort_t Edge,
            // UShort_t Side)
        }
    }

    fNEvent += 1;
    return kTRUE;
}

void R3BRpcReader::Reset()
{
    // Reset the output array
    fArrayRpc->Clear();
}

ClassImp(R3BRpcReader)
