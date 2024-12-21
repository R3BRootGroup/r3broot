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
#include "R3BSci8MappedData.h"
#include "R3BSci8Reader.h"
#include "TClonesArray.h"
extern "C"
{
#include "ext_data_client.h"
#include "ext_h101_sci8.h"
}
#include "TMath.h"
#define IS_NAN(x) TMath::IsNaN(x)
#define NUM_SCI8_DETECTORS 1
#define NUM_SCI8_CHANNELS 2
#include <iostream>

using namespace std;

R3BSci8Reader::R3BSci8Reader(EXT_STR_h101_SCI8* data, UInt_t offset)
    : R3BReader("R3BSci8Reader")
    , fData(data)
    , fOffset(offset)
    , fLogger(FairLogger::GetLogger())
    , fArray(new TClonesArray("R3BSci8MappedData"))
{
}

R3BSci8Reader::~R3BSci8Reader() {}

Bool_t R3BSci8Reader::Init(ext_data_struct_info* a_struct_info)
{

    int ok;

    EXT_STR_h101_SCI8_ITEMS_INFO(ok, *a_struct_info, fOffset, EXT_STR_h101_SCI8, 0);

    if (!ok)
    {
        perror("ext_data_struct_info_item");
        LOG(error) << "Failed to setup structure information.";
        return kFALSE;
    }

    // Register output array in tree
    FairRootManager::Instance()->Register("Sci8Mapped", "Land", fArray, kTRUE);
    fArray->Clear();

    // clear struct_writer's output struct. Seems ucesb doesn't do that
    // for channels that are unknown to the current ucesb config.
    auto* data = reinterpret_cast<EXT_STR_h101_SCI8_onion*>(fData);

    for (int d = 0; d < NUM_SCI8_DETECTORS; d++)
    {
        data->SCIEIGHT_VTFM = 0;
        data->SCIEIGHT_VTCM = 0;
        /* no TAMEX present
                data->SCIEIGHT_TTFLM=0;
                data->SCIEIGHT_TTFTM=0;
                data->SCIEIGHT_TTCLM=0;
                data->SCIEIGHT_TTCTM=0;
        */
    }
    return kTRUE;
}

Bool_t R3BSci8Reader::R3BRead()
{
    // Convert plain raw data to multi-dimensional array
    auto* data = reinterpret_cast<EXT_STR_h101_SCI8_onion*>(fData);

    /*
     * For variable definition, see structure EXT_STR_h101_SCI8_onion_t
     * in ext_str_h101_sci8_tamex.h
     *** VFTX DATA ***
     * VTF = Size of the array TFv contaning fine VFTX times
     * VTFv = Array containing the actual data on the fine  times
     * VTFM = No of channels having time
     * VTFMI = Array of TFM size containing the channel numbers of each channel with data
     * VTFME = Array of TFM size containing the index of the first element of the next channel in data array TFv;
     *         TFME[i]-FTME[i-1] = number of data for channel i;
     *
     * VTC = Size of the array TCv contaning coarse VFTX times
     * VTCv = Array containing the actual data on the coarse times
     * VTCM = No of channels having coarse time
     * VTCMI = Array of TCM size containing the channel numbers of each channel with data
     * VTCME = Array of TCM size containing the index of the first element of the next channel in data array TCv
     *
     *** TAMEX DATA *** Added by Aleksandra, Oct. 2016
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

    Bool_t fprint = false;

    for (int d = 0; d < NUM_SCI8_DETECTORS; d++)
    {

        Int_t Sum = data->SCIEIGHT_VTF; // no Tamex present+data->SCIEIGHT_TTFT+data->SCIEIGHT_TTFL;
        // if(data->S8TTFT != data->S8TTFL) fprint = true;
        // if(fNEvents == 9698 || fNEvents == 9701 || fNEvents == 9704) fprint = true;
        // First, we prepare time arrays for VFTX

        // VFTX first:
        uint32_t numChannels = data->SCIEIGHT_VTFM; // not necessarly number of hits! (b/c multi hit)
        // loop over channels
        uint32_t curChannelStart = 0; // index in v for first item of current channel
        Double_t mean_coarse_vftx = 0.;
        int sum_coarse_vftx = 0;
        // First get the average coarse time to shift all coarse counters in the same cycle (by calculateing, in the
        // second step, deviations from the mean value. If the coarse time is smaller than mean value by more than 200,
        // then coarse counter was reseted, and thus, to its value 8192 (in case of VFTX) will be added.
        for (int i = 0; i < numChannels; i++)
        {
            uint32_t channel = data->SCIEIGHT_VTFMI[i];          // = 1..8
            uint32_t nextChannelStart = data->SCIEIGHT_VTFME[i]; // index in v for first item of next channel

            for (int j = curChannelStart; j < nextChannelStart; j++)
            {

                int coarse_vftx = data->SCIEIGHT_VTCv[j];

                mean_coarse_vftx = mean_coarse_vftx + coarse_vftx;
                sum_coarse_vftx = sum_coarse_vftx + 1;
            }
            curChannelStart = nextChannelStart;
        }
        mean_coarse_vftx = mean_coarse_vftx / float(sum_coarse_vftx);

        curChannelStart = 0;
        for (int i = 0; i < numChannels; i++) // VFTX, now do the mapping
        {
            uint32_t channel = data->SCIEIGHT_VTFMI[i];          // = 1..8
            uint32_t nextChannelStart = data->SCIEIGHT_VTFME[i]; // index in v for first item of next channel

            for (int j = curChannelStart; j < nextChannelStart; j++)
            {

                int coarse_vftx = data->SCIEIGHT_VTCv[j];
                if ((mean_coarse_vftx - float(coarse_vftx)) > 200.)
                    coarse_vftx = coarse_vftx + 8192;

                if (fprint)
                    cout << "SCI8 READER VFTX: " << fNEvents << ", " << Sum << ", " << channel << ", "
                         << data->SCIEIGHT_VTFv[j] << ", " << data->SCIEIGHT_VTCv[j] << ", " << coarse_vftx << ", "
                         << mean_coarse_vftx << endl;

                new ((*fArray)[fArray->GetEntriesFast()])
                    R3BSci8MappedData(d + 1,                  // detector number
                                      channel,                // channel number: 1-8
                                      0,                      // VFTX (0),TAMEX leading (1), TAMEX trailing (2)
                                      data->SCIEIGHT_VTFv[j], // VFTX fine time
                                      coarse_vftx             // VFTX coarse time
                    );
            }
            curChannelStart = nextChannelStart;
        }
    }
    fNEvents += 1;

    return kTRUE;
}

void R3BSci8Reader::FinishTask() {}

void R3BSci8Reader::Reset()
{
    // Reset the output array
    fArray->Clear();
}

ClassImp(R3BSci8Reader)
