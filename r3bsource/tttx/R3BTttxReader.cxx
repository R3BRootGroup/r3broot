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
#include <TClonesArray.h>

#include "R3BLogger.h"
#include "R3BTttxMappedData.h"
#include "R3BTttxReader.h"

extern "C"
{
#include "ext_data_client.h"
#include "ext_h101_ttt10.h"
}

using namespace std;

R3BTttxReader::R3BTttxReader(EXT_STR_h101_TTTX* data, size_t offset)
    : R3BReader("R3BTttxReader")
    , fData(data)
    , fOffset(offset)
    , fOnline(kFALSE)
    , fDets(2)
    , fChannels(32)
    , fTref(1)
    , fTtrig(1)
    , fPileup(kFALSE)
    , fArray(new TClonesArray("R3BTttxMappedData"))
{
}

R3BTttxReader::~R3BTttxReader()
{
    R3BLOG(debug1, "Destructor");
    if (fArray)
    {
        delete fArray;
    }
}

Bool_t R3BTttxReader::Init(ext_data_struct_info* a_struct_info)
{
    Int_t ok;
    R3BLOG(info, "");
    EXT_STR_h101_TTTX_ITEMS_INFO(ok, *a_struct_info, fOffset, EXT_STR_h101_TTTX, 0);
    if (!ok)
    {
        R3BLOG(error, "Failed to setup structure information.");
        return kFALSE;
    }

    // Register output array in tree
    FairRootManager::Instance()->Register("tttxMappedData", "MappedTTTX", fArray, !fOnline);

    Reset();
    memset(fData, 0, sizeof *fData);

    // clear struct_writer's output struct. Seems ucesb doesn't do that
    // for channels that are unknown to the current ucesb config.
    // auto* data = (EXT_STR_h101_TTTX_onion*)fData;
    auto* data = reinterpret_cast<EXT_STR_h101_TTTX_onion*>(fData);
    if (data == nullptr)
    {
        R3BLOG(error, "tttxMappedData NOT FOUND!");
    }
    for (int i_det = 0; i_det < fDets; i_det++)
    {
        data->TTTX_ID[i_det].EM = 0;
        data->TTTX_ID[i_det].TM = 0;
        data->TTTX_ID[i_det].TREFM = 0;
        data->TTTX_ID[i_det].TTRIGM = 0;
    }

    return kTRUE;
}

Bool_t R3BTttxReader::R3BRead()
{
    // Convert plain raw data to multi-dimensional array
    auto* data = reinterpret_cast<EXT_STR_h101_TTTX_onion*>(fData);
    // auto* data = (EXT_STR_h101_TTTX_onion*)fData;

    // loop over all detectors
    for (UShort_t i_det = 0; i_det < fDets; i_det++)
    {
        ReadData(data, i_det);
    }
    return kTRUE;
}

void R3BTttxReader::Reset()
{
    // Reset the output array
    fArray->Clear();
}

Bool_t R3BTttxReader::ReadData(EXT_STR_h101_TTTX_onion* data, UShort_t i_det)
{
    Bool_t pileupFLAG = kFALSE;
    Bool_t overflowFLAG = kFALSE;
    UInt_t multPerChannel[fChannels + fTref + fTtrig];
    for (int i_ch = 0; i_ch < fChannels + fTref + fTtrig; i_ch++)
    {
        multPerChannel[i_ch] = 0;
    }
    // --- -------------------------------------- --- //
    // --- NUMBER OF CHANNELS TREF, TTRIG, Channels --- //
    // --- -------------------------------------- --- //
    UShort_t nChannelsTref = data->TTTX_ID[i_det].TREFM;
    UShort_t nChannelsTtrig = data->TTTX_ID[i_det].TTRIGM;
    UShort_t nChannelsEnergy = data->TTTX_ID[i_det].EM;
    UShort_t nChannelsTime = data->TTTX_ID[i_det].TM;

    /*
    if ((nChannelsEnergy>0)||(nChannelsTref>0)) {
      std::cout << "------------------------------" << std::endl;
      std::cout << "R3BTttxReader::ReadDatai=()" << std::endl;
      std::cout << "------------------------------" << std::endl;
      std::cout <<" * nChannelsTref = "<< nChannelsTref <<std::endl;
      std::cout <<"  * nChannelsTtrig = "<< nChannelsTtrig <<std::endl;
      std::cout << "   * nChannelsEnergy = " << nChannelsEnergy << std::endl;
      std::cout << "   * nChannelsTime = " << nChannelsTime << std::endl;
      std::cout << "------------------------------" << std::endl;
    }
    */

    // --- ----------------- --- //
    // --- TTTX MAPPED DATA --- //
    // --- ----------------- --- //

    // --- TREF --- //
    uint32_t curTref = 0;
    uint32_t nextTref = 0;
    UShort_t idChannelTref = 0;

    for (UShort_t i_tref = 0; i_tref < nChannelsTref; i_tref++)
    {
        // TREFMI gives the 1-based Tref number
        idChannelTref = data->TTTX_ID[i_det].TREFMI[i_tref] + fChannels - 1;
        nextTref = data->TTTX_ID[i_det].TREFME[i_tref];
        multPerChannel[idChannelTref] = nextTref - curTref;
        for (int i_hit = curTref; i_hit < nextTref; i_hit++)
        {
            new ((*fArray)[fArray->GetEntriesFast()])
                R3BTttxMappedData(i_det + 1, idChannelTref + 1, data->TTTX_ID[i_det].TREFv[i_hit], 0, false, false);
            curTref = nextTref;
        }
    }

    // --- TTRIG --- //
    uint32_t curTtrig = 0;
    uint32_t nextTtrig = 0;
    UShort_t idChannelTtrig = 0;

    for (UShort_t i_trig = 0; i_trig < nChannelsTtrig; i_trig++)
    {
        // again TTRIGMI is 1 based
        idChannelTtrig = data->TTTX_ID[i_det].TTRIGMI[i_trig] + fChannels - 1 + fTref;
        nextTtrig = data->TTTX_ID[i_det].TTRIGME[i_trig];
        multPerChannel[idChannelTtrig] = nextTtrig - curTtrig;
        for (int i_hit = curTtrig; i_hit < nextTtrig; i_hit++)
        {
            new ((*fArray)[fArray->GetEntriesFast()])
                R3BTttxMappedData(i_det + 1, idChannelTtrig + 1, data->TTTX_ID[i_det].TTRIGv[i_hit], 0, false, false);
            curTtrig = nextTtrig;
        }
    }

    // --- Channels --- //
    // return the number of channels in the det with data
    // 0<=nChannelsEnergy<32
    // 0<=nChannelsTime<32
    if (nChannelsEnergy != nChannelsTime)
        LOG(error) << "R3BTttxReader::ReadData error ! NOT THE SAME NUMBER OF CHANNELS HITTED IN ENERGY () AND TIME ()";

    // ENERGY AND TIME ARE SORTED
    uint32_t curChannelTimeStart = 0;
    // uint32_t curChannelEnergyStart = 0;
    for (UShort_t i_ch = 0; i_ch < nChannelsTime; i_ch++)
    {
        // EMI and TMI give the 1-based channel number
        UShort_t idChannelTime = data->TTTX_ID[i_det].TMI[i_ch] - 1;
        UShort_t idChannelEnergy = data->TTTX_ID[i_det].EMI[i_ch] - 1;

        if (idChannelEnergy != idChannelTime)
        {
            LOG(error) << "R3BTttxReader::ReadData error ! MISMATCH FOR CHANNEL ID IN ENERGY #" << idChannelEnergy
                       << " AND TIME #" << idChannelTime;
        }
        uint32_t nextChannelTimeStart = data->TTTX_ID[i_det].TME[i_ch];
        // uint32_t nextChannelEnergyStart = data->TTTX_ID[i_det].EME[i_ch];
        multPerChannel[idChannelTime] = nextChannelTimeStart - curChannelTimeStart;
        for (int i_hit = curChannelTimeStart; i_hit < nextChannelTimeStart; i_hit++)
        {
            // Attention, here the numbering is 0-based for det and channels
            pileupFLAG = (data->TTTX_ID[i_det].Ev[i_hit] & 0x00040000) >> 18;
            overflowFLAG = (data->TTTX_ID[i_det].Ev[i_hit] & 0x00080000) >> 19;
            if (pileupFLAG == kFALSE || !fPileup)
                new ((*fArray)[fArray->GetEntriesFast()]) R3BTttxMappedData(i_det + 1,
                                                                            idChannelEnergy + 1,
                                                                            data->TTTX_ID[i_det].Tv[i_hit],
                                                                            data->TTTX_ID[i_det].Ev[i_hit],
                                                                            pileupFLAG,
                                                                            overflowFLAG);
        }
        // curChannelEnergyStart = nextChannelEnergyStart;
        curChannelTimeStart = nextChannelTimeStart;
    } // end of loop over the channels from 1 to 32

    return kTRUE;
}

ClassImp(R3BTttxReader)
