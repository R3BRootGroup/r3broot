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

#include "R3BPdcReader.h"
#include "FairLogger.h"
#include "FairRootManager.h"
#include "R3BPdcMappedData.h"
#include "TClonesArray.h"
using namespace std;
extern "C"
{
#include "ext_data_client.h"
#include "ext_h101_pdc.h"
}

#define LENGTH(x) (sizeof x / sizeof *x)

R3BPdcReader::R3BPdcReader(EXT_STR_h101_PDC* data, UInt_t a_offset)
    : R3BReader("R3BPdcReader")
    , fData(data)
    , fOffset(a_offset)
    , fMappedArray(new TClonesArray("R3BPdcMappedData"))
{
}

R3BPdcReader::~R3BPdcReader() {}

Bool_t R3BPdcReader::Init(ext_data_struct_info* a_struct_info)
{
    int ok;

    EXT_STR_h101_PDC_ITEMS_INFO(ok, *a_struct_info, fOffset, EXT_STR_h101_PDC, 0);

    if (!ok)
    {
        perror("ext_data_struct_info_item");
        LOG(error) << "Failed to setup structure information.";
        return kFALSE;
    }

    FairRootManager::Instance()->Register("PdcMapped", "Land", fMappedArray, kTRUE);

    auto* data = reinterpret_cast<EXT_STR_h101_PDC_onion*>(fData);
    memset(data, 0, sizeof *data);
    return kTRUE;
}

Bool_t R3BPdcReader::R3BRead()
{
    // Convert plain raw data to multi-dimensional array
    auto* data = reinterpret_cast<EXT_STR_h101_PDC_onion*>(fData);

    // puts("Event");
    for (uint32_t p = 0; p < LENGTH(data->PDC_P); p++)
    {
        auto const& side = data->PDC_P[p];

        //
        // ctdc.
        //

        auto numChannels = side.TLCM;
        uint32_t curChannelStart = 0;
        for (uint32_t i = 0; i < numChannels; i++)
        {
            uint32_t channel = side.TLCMI[i];
            uint32_t nextChannelStart = side.TLCME[i];
            for (uint32_t j = curChannelStart; j < nextChannelStart; j++)
            {
                // cout << "p: " << p+1 << "  channel: " << channel << "  edge: 1   coarse: " <<   side.TLCv[j] << endl;
                new ((*fMappedArray)[fMappedArray->GetEntriesFast()])
                    R3BPdcMappedData(p + 1, channel, 1, side.TLCv[j], side.TLFv[j]);
            }
            curChannelStart = nextChannelStart;
        }

        // Trailing.
        numChannels = side.TTCM;
        curChannelStart = 0;
        for (uint32_t i = 0; i < numChannels; i++)
        {
            uint32_t channel = side.TTCMI[i];
            uint32_t nextChannelStart = side.TTCME[i];
            for (uint32_t j = curChannelStart; j < nextChannelStart; j++)
            {
                // cout << "p: " << p+1 << "  channel: " << channel << "  edge: 2   coarse: " <<   side.TTCv[j] << endl;
                new ((*fMappedArray)[fMappedArray->GetEntriesFast()])
                    R3BPdcMappedData(p + 1, channel, 2, side.TTCv[j], side.TTFv[j]);
            }
            curChannelStart = nextChannelStart;
        }
    }

    // Leading ctdc trigger times.
    {
        auto numChannels = data->PDC_TRIGLC;
        for (uint32_t i = 0; i < numChannels; i++)
        {
            uint32_t channel = data->PDC_TRIGLCI[i];
            new ((*fMappedArray)[fMappedArray->GetEntriesFast()])
                R3BPdcMappedData(LENGTH(data->PDC_P) + 1, channel, 1, data->PDC_TRIGLCv[i], data->PDC_TRIGLFv[i]);
        }
    }

    return kTRUE;
}

void R3BPdcReader::Reset() { fMappedArray->Clear(); }

ClassImp(R3BPdcReader)
