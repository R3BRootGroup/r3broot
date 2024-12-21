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

#include "R3BLogger.h"
#include "R3BRoluMappedData.h"
#include "R3BRoluReaderOnline.h"

#include <TClonesArray.h>
#include <TMath.h>
#include <iostream>

/**
 ** ext_h101_rolu.h was created by running
 ** $unpacker --ntuple=STRUCT_HH,RAW:ROLU,id=h101_ROLU,NOTRIGEVENTNO,ext_h101_rolu.h
 **/

extern "C"
{
#include "ext_data_client.h"
#include "ext_h101_rolu_online.h"
}

#define IS_NAN(x) TMath::IsNaN(x)
#define NUM_ROLU_DETECTORS (sizeof fData->ROLU / sizeof fData->ROLU[0])

using namespace std;

R3BRoluReaderOnline::R3BRoluReaderOnline(EXT_STR_h101_ROLU_NOTRIG_onion* data, size_t offset)
    : R3BReader("R3BRoluReaderOnline")
    , fData(data)
    , fNbDet(NUM_ROLU_DETECTORS)
    , fOffset(offset)
    , fOnline(kFALSE)
    , fArray(new TClonesArray("R3BRoluMappedData"))
{
}

R3BRoluReaderOnline::~R3BRoluReaderOnline()
{
    if (fArray)
    {
        delete fArray;
    }
}

bool R3BRoluReaderOnline::Init(ext_data_struct_info* a_struct_info)
{
    int ok;
    R3BLOG(info, "");
    EXT_STR_h101_ROLU_NOTRIG_ITEMS_INFO(ok, *a_struct_info, fOffset, EXT_STR_h101_ROLU_NOTRIG, 0);
    if (!ok)
    {
        R3BLOG(error, "Failed to setup ROLU structure information.");
        return kFALSE;
    }

    // Register output array in tree
    FairRootManager::Instance()->Register("RoluMapped", "Rolu Mapped", fArray, !fOnline);
    Reset();

    // clear struct_writer's output struct. Seems ucesb doesn't do that
    // for channels that are unknown to the current ucesb config.
    R3BLOG(info, "ROLU num Dets: " << fNbDet);
    for (uint32_t det = 0; det < fNbDet; det++)
    {
        fData->ROLU[det].TTFLM = 0;
        fData->ROLU[det].TTFTM = 0;
        fData->ROLU[det].TTCLM = 0;
        fData->ROLU[det].TTCTM = 0;
    }
    return kTRUE;
}

bool R3BRoluReaderOnline::R3BRead()
{
    // Convert plain raw data to multi-dimensional array
    EXT_STR_h101_ROLU_NOTRIG_onion* data = (EXT_STR_h101_ROLU_NOTRIG_onion*)fData;

    /*
     * For variable definition, see structure EXT_STR_h101_ROLU_NOTRIG_onion_t
     * in ext_str_h101_rolu.h
     * TTFL = Size of the array TFLv contaning fine leading times
     * TTFLv = Array containing the actual data on the fine leading times
     * TTFLM = No of channels having fine leading time
     * TTFLMI = Array of TFLM size containing the channel numbers of each channel with data
     * TTFLME = Array of TFLM size containing the index of the first element of the next channel in data array TFLv
     *
     * TTCL = Size of the array TCLv contaning coarse leading times
     * TTCLv = Array containing the actual data on the coarse leading times
     * TTCLM = No of channels having coarse leading time
     * TTCLMI = Array of TCLM size containing the channel numbers of each channel with data
     * TTCLME = Array of TCLM size containing the index of the first element of the next channel in data array TCLv
     *
     * The same logic is for trailing times: TTFT, TTFTv, ..., TTCTMI, TTCTME
     */

    // loop over all detectors

    for (uint32_t det = 0; det < fNbDet; det++)
    {
        // Coarse counter reset recovery:
        // Divide the course counter range 0..8191 into four pieces.
        // If we have _any_ hit in the uppermost quarter, then all hits in the
        // lowermost quarter will later be shifted up by the course counter range.

        bool do_increment = false;
        uint32_t const c_tamex3_range = 2048;

        if (data->ROLU[det].TTCL < 1)
            return kFALSE;

        uint32_t numData = data->ROLU[det].TTCL;
        uint32_t numChannels = data->ROLU[det].TTCLM;
        for (uint32_t chan = 0; chan < numData; chan++)
        {
            uint32_t coarse_leading = data->ROLU[det].TTCLv[chan];
            if (coarse_leading > 3 * c_tamex3_range / 4)
            {
                do_increment = true;
                break;
            }
        }
        if (!do_increment)
        {
            numData = data->ROLU[det].TTCT;
            numChannels = data->ROLU[det].TTCTM;
            for (uint32_t chan = 0; chan < numData; chan++)
            {
                uint32_t coarse_trailing = data->ROLU[det].TTCTv[chan];
                if (coarse_trailing > 3 * c_tamex3_range / 4)
                {
                    do_increment = true;
                    break;
                }
            }
        }

        // Leading mapping.
        numChannels = data->ROLU[det].TTCLM;
        uint32_t curChannelStart = 0;
        for (uint32_t chan = 0; chan < numChannels; chan++)
        {
            uint32_t channel = data->ROLU[det].TTCLMI[chan];
            uint32_t nextChannelStart = data->ROLU[det].TTCLME[chan];
            for (uint32_t nextchan = curChannelStart; nextchan < nextChannelStart; nextchan++)
            {
                uint32_t coarse_leading = data->ROLU[det].TTCLv[nextchan];
                if (do_increment && coarse_leading < c_tamex3_range / 4)
                {
                    coarse_leading += c_tamex3_range;
                }
                new ((*fArray)[fArray->GetEntriesFast()])
                    R3BRoluMappedData(det + 1, channel, 0, data->ROLU[det].TTFLv[nextchan], coarse_leading);
            }
            curChannelStart = nextChannelStart;
        }

        //
        // Trailing.
        // Matched against leading edges.
        //
        numChannels = data->ROLU[det].TTCTM;
        curChannelStart = 0;
        for (uint32_t chan = 0; chan < numChannels; chan++)
        {
            uint32_t channel = data->ROLU[det].TTCTMI[chan];
            uint32_t nextChannelStart = data->ROLU[det].TTCTME[chan];
            for (uint32_t nextchan = curChannelStart; nextchan < nextChannelStart; nextchan++)
            {
                uint32_t coarse_trailing = data->ROLU[det].TTCTv[nextchan];
                if (do_increment && coarse_trailing < c_tamex3_range / 4)
                {
                    coarse_trailing += c_tamex3_range;
                }

                // Find a leading edge to pair up with.
                int nHit = fArray->GetEntriesFast();
                for (int index = 0; index < nHit; index++)
                {
                    R3BRoluMappedData const* hit = dynamic_cast<R3BRoluMappedData*>(fArray->At(index));

                    unsigned int const iTypeL = hit->GetType();
                    unsigned int const iCha = hit->GetChannel();

                    // Only consider leading data.
                    if (0 == iTypeL && iCha == channel)
                    {
                        uint32_t coarse_leading = hit->GetTimeCoarse();
                        int32_t tot = coarse_trailing - coarse_leading;

                        new ((*fArray)[fArray->GetEntriesFast()])
                            R3BRoluMappedData(det + 1, channel, 1, data->ROLU[det].TTFTv[nextchan], coarse_trailing);
                    }
                }
            }
            curChannelStart = nextChannelStart;
        }
    }
    fNEvents += 1;

    return kTRUE;
}

void R3BRoluReaderOnline::Reset()
{
    // Reset the output array
    fArray->Clear();
}

ClassImp(R3BRoluReaderOnline)
