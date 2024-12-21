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

#include "R3BLosReader.h"
#include "R3BEventHeader.h"
#include "R3BLogger.h"
#include "R3BLosMappedData.h"

#include <FairLogger.h>
#include <FairRootManager.h>

#include <TClonesArray.h>
#include <TMath.h>

/**
 ** ext_h101_los.h was created by running
 ** $unpacker --ntuple=STRUCT_HH,RAW:LOS,id=h101_LOS,NOTRIGEVENTNO,ext_h101_los.h
 **/

extern "C"
{
#include "ext_data_client.h"
#include "ext_h101_los.h"
}

#define IS_NAN(x) TMath::IsNaN(x)
#define NUM_LOS_DETECTORS (sizeof data->LOS / sizeof data->LOS[0])
#define NUM_LOS_CHANNELS 8

using namespace std;

R3BLosReader::R3BLosReader(EXT_STR_h101_LOS* data, size_t offset)
    : R3BReader("R3BLosReader")
    , fNEvents(0)
    , fData(data)
    , fOffset(offset)
    , fOnline(kFALSE)
    , fSkiptriggertimes(kFALSE)
    , fArray(new TClonesArray("R3BLosMappedData"))
    , fArrayTrigger(new TClonesArray("R3BLosMappedData"))
{
}

R3BLosReader::~R3BLosReader()
{
    R3BLOG(debug1, "");
    if (fArray)
    {
        delete fArray;
    }
    if (fArrayTrigger)
    {
        delete fArrayTrigger;
    }
}

Bool_t R3BLosReader::Init(ext_data_struct_info* a_struct_info)
{
    Int_t ok;
    R3BLOG(info, "");
    EXT_STR_h101_LOS_ITEMS_INFO(ok, *a_struct_info, fOffset, EXT_STR_h101_LOS, 0);
    if (!ok)
    {
        R3BLOG(fatal, "Failed to setup structure information.");
        return kFALSE;
    }

    // try to get a handle on the EventHeader. EventHeader may not be
    // present though and hence may be null. Take care when using.
    FairRootManager* mgr = FairRootManager::Instance();
    R3BLOG_IF(fatal, !mgr, "FairRootManager not found");

    header = dynamic_cast<R3BEventHeader*>(mgr->GetObject("EventHeader."));
    if (!header)
    {
        header = dynamic_cast<R3BEventHeader*>(mgr->GetObject("R3BEventHeader"));
        R3BLOG(warn, "EventHeader. not found");
    }

    // Register output array in tree
    mgr->Register("LosMapped", "Los Mapped", fArray, !fOnline);
    if (!fSkiptriggertimes)
    {
        mgr->Register("LosTriggerMapped", "Los Trigger Mapped", fArrayTrigger, !fOnline);
    }
    else
    {
        fArrayTrigger = NULL;
    }
    Reset();

    // clear struct_writer's output struct. Seems ucesb doesn't do that
    // for channels that are unknown to the current ucesb config.
    auto* data = reinterpret_cast<EXT_STR_h101_LOS_onion*>(fData);
    for (uint32_t d = 0; d < NUM_LOS_DETECTORS; d++)
    {
        data->LOS[d].VTFM = 0;
        data->LOS[d].VTCM = 0;
        data->LOS[d].TTFLM = 0;
        data->LOS[d].TTFTM = 0;
        data->LOS[d].TTCLM = 0;
        data->LOS[d].TTCTM = 0;
        data->LOS[d].VTRIGC = 0;
        data->LOS[d].VTRIGF = 0;
        data->LOS[d].TTRIGCL = 0;
        data->LOS[d].TTRIGFL = 0;
    }

    return kTRUE;
}

Bool_t R3BLosReader::R3BRead()
{
    R3BLOG(debug1, "Event data.");
    // Convert plain raw data to multi-dimensional array
    auto* data = reinterpret_cast<EXT_STR_h101_LOS_onion*>(fData);

    /*
     * For variable definition, see structure EXT_STR_h101_LOS_onion_t
     * in ext_str_h101_los_tamex.h
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
    //  cout<<"*************** READER START ******************"<<endl;
    //  cout<<"LOS READER: "<<NUM_LOS_DETECTORS<<endl;

    for (uint32_t d = 0; d < NUM_LOS_DETECTORS; d++)
    {

        Int_t nsumv = 0, nsuml = 0, nsumt = 0;

        if (data->LOS[d].VTF > 0)
        {
            //  cout<<"FINE: iDet "<<d+1<<", "<<data->LOS[d].VTF<<", "<<data->LOS[d].TTFL<<",
            //  "<<data->LOS[d].TTFT<<endl;
            // cout<<"COARSE: iDet "<<d+1<<", "<<data->LOS[d].VTC<<", "<<data->LOS[d].TTCL<<",
            // "<<data->LOS[d].TTCT<<endl;
        }
        /*
            if(data->LOS[d].VTF != data->LOS[d].VTC) return kFALSE;
            if(data->LOS[d].TTFL != data->LOS[d].TTCL) return kFALSE;
            if(data->LOS[d].TTFT != data->LOS[d].TTCT) return kFALSE;
            if(data->LOS[d].VTF < data->LOS[d].TTFL) return kFALSE;
            if(data->LOS[d].TTFL != data->LOS[d].TTFT) return kFALSE;

            if(data->LOS[d].VTFM != data->LOS[d].VTCM) return kFALSE;
            if(data->LOS[d].TTFLM != data->LOS[d].TTCLM) return kFALSE;
            if(data->LOS[d].TTFTM != data->LOS[d].TTCTM) return kFALSE;
            if(data->LOS[d].VTFM < data->LOS[d].TTFLM) return kFALSE;
            if(data->LOS[d].TTFLM != data->LOS[d].TTFTM) return kFALSE;
        */
        //
        // VFTX2.
        //

        uint32_t const c_vftx2_range = 8192;
        uint32_t numData = data->LOS[d].VTF;

        // Coarse counter reset recovery:
        // Divide the course counter range 0..8191 into four pieces.
        // If we have _any_ hit in the uppermost quarter, then all hits in the
        // lowermost quarter will later be shifted up by the course counter range.
        bool do_increment = false;
        for (uint32_t i = 0; i < numData; i++)
        {
            uint32_t coarse_vftx = 0. / 0.;
            coarse_vftx = data->LOS[d].VTCv[i];
            if (coarse_vftx > 3 * c_vftx2_range / 4)
            {
                // FIXME:
                // This algo gives problems with vftx calibration and data processing.
                do_increment = false; // true;
                break;
            }
        }

        // VFTX2 Mapping.
        uint32_t numChannels = data->LOS[d].VTFM;
        uint32_t curChannelStart = 0;
        for (uint32_t i = 0; i < numChannels; i++)
        {
            uint32_t channel = data->LOS[d].VTFMI[i]; // = 1..8
            uint32_t nextChannelStart = data->LOS[d].VTFME[i];
            for (uint32_t j = curChannelStart; j < nextChannelStart; j++)
            {
                uint32_t coarse_vftx = data->LOS[d].VTCv[j];
                if (do_increment && coarse_vftx < c_vftx2_range / 4)
                {
                    coarse_vftx += c_vftx2_range;
                }

                //      if(header->GetTrigger() != 1) cout<< "Trigger in LosReader: "<<header->GetTrigger()<<endl;

                new ((*fArray)[fArray->GetEntriesFast()])
                    R3BLosMappedData(d + 1,                // detector number
                                     channel,              // channel number: 1-8
                                     0,                    // VFTX (0),TAMEX leading (1), TAMEX trailing (2)
                                     data->LOS[d].VTFv[j], // VFTX fine time
                                     coarse_vftx           // VFTX coarse time
                    );

                //    cout<<"Reader VFTX: "<<d+1<<", "<<channel<<"; "<<data->LOS[d].VTFv[j]<<", "<<coarse_vftx <<endl;
                nsumv += 1;
            }
            curChannelStart = nextChannelStart;
        }

        //
        // TAMEX3.
        //

        // Coarse counter recovery like VFTX2, but consider both edges.
        do_increment = false;
        uint32_t const c_tamex3_range = 2048;

        numData = data->LOS[d].TTCL;
        numChannels = data->LOS[d].TTCLM;
        for (uint32_t i = 0; i < numData; i++)
        {
            uint32_t coarse_leading = data->LOS[d].TTCLv[i];
            if (coarse_leading > 3 * c_tamex3_range / 4)
            {
                do_increment = true;
                break;
            }
        }
        if (!do_increment)
        {
            numData = data->LOS[d].TTCT;
            numChannels = data->LOS[d].TTCTM;
            for (uint32_t i = 0; i < numData; i++)
            {
                uint32_t coarse_trailing = data->LOS[d].TTCTv[i];
                if (coarse_trailing > 3 * c_tamex3_range / 4)
                {
                    do_increment = true;
                    break;
                }
            }
        }

        // TAMEX3 leading mapping.
        numChannels = data->LOS[d].TTCLM;
        curChannelStart = 0;
        for (uint32_t i = 0; i < numChannels; i++)
        {
            uint32_t channel = data->LOS[d].TTCLMI[i];
            uint32_t nextChannelStart = data->LOS[d].TTCLME[i];
            for (uint32_t j = curChannelStart; j < nextChannelStart; j++)
            {
                uint32_t coarse_leading = data->LOS[d].TTCLv[j];
                if (do_increment && coarse_leading < c_tamex3_range / 4)
                {
                    coarse_leading += c_tamex3_range;
                }
                new ((*fArray)[fArray->GetEntriesFast()])
                    R3BLosMappedData(d + 1, channel, 1, data->LOS[d].TTFLv[j], coarse_leading);

                //    cout<<"Reader TAMEX leading: "<<d + 1<<", "<<channel<<", "<< data->LOS[d].TTFLv[j]<<", "<<
                //    coarse_leading<<endl;
                nsuml += 1;
            }
            curChannelStart = nextChannelStart;
        }

        //
        // TAMEX3 trailing.
        // Matched against leading edges.
        //

        numChannels = data->LOS[d].TTCTM;
        curChannelStart = 0;
        for (uint32_t i = 0; i < numChannels; i++)
        {
            uint32_t channel = data->LOS[d].TTCTMI[i];
            uint32_t nextChannelStart = data->LOS[d].TTCTME[i];
            for (uint32_t j = curChannelStart; j < nextChannelStart; j++)
            {
                uint32_t coarse_trailing = data->LOS[d].TTCTv[j];
                if (do_increment && coarse_trailing < c_tamex3_range / 4)
                {
                    coarse_trailing += c_tamex3_range;
                }

                // Find a leading edge to pair up with.
                Int_t n = fArray->GetEntriesFast();
                for (Int_t k = 0; k < n; k++)
                {
                    R3BLosMappedData const* hit = dynamic_cast<R3BLosMappedData*>(fArray->At(k));

                    UInt_t const iTypeL = hit->GetType();
                    UInt_t const iCha = hit->GetChannel();

                    // Only consider leading data.
                    if (1 == iTypeL && iCha == channel)
                    {
                        uint32_t coarse_leading = hit->GetTimeCoarse();
                        int32_t tot = coarse_trailing - coarse_leading;
                        // 30 units -> 30 * 5 = 150 ns.
                        if ((tot >= 0))
                        { // tot <= 135) &&
                            new ((*fArray)[fArray->GetEntriesFast()])
                                R3BLosMappedData(d + 1, channel, 2, data->LOS[d].TTFTv[j], coarse_trailing);
                            nsumt += 1;
                            //     cout<<"Reader TAMEX trailing: "<<d + 1<<", "<<channel<<", "<<
                            //     data->LOS[d].TTFTv[j]<<", "<< coarse_trailing<<endl;

                            break;
                        }
                    }
                }
            }
            curChannelStart = nextChannelStart;
        }

        //
        // MTDC32
        //
#if 0
        numChannels = data->LOS[d].MTM;
        curChannelStart = 0;
        // cout<<data->LOS[d].MTM<<endl;
        for (uint32_t i = 0; i < numChannels; i++)
        {
            uint32_t channel = data->LOS[d].MTMI[i]; // = 1..8
            uint32_t nextChannelStart = data->LOS[d].MTME[i];
            /*   for (uint32_t j = curChannelStart; j < nextChannelStart; j++) {
                 new ((*fArray)[fArray->GetEntriesFast()])
                     R3BLosMappedData(
                         d + 1,               // detector number
                         channel,             // channel number: 1-8
                         3,                   // VFTX (0),TAMEX leading (1), TAMEX trailing (2), MTDC32 (3)
                         data->LOS[d].MTv[j], // MTDC32 time
                         0
                         );

                   cout<<"MTDC: "<<channel<<", "<<data->LOS[d].MTv[j] <<endl;

               }*/
            curChannelStart = nextChannelStart;
        }
#endif
        if (data->LOS[d].VTF > 0)
        {
            // cout<<"nsumv & data->LOS[d].VTF "<<d+1<<"; "<<nsumv<<"; "<< data->LOS[d].VTF<<endl;
            // cout<<"nsumt & data->LOS[d].TTFL "<<d+1<<"; "<<nsuml<<"; "<< data->LOS[d].TTFL<<endl;
            // cout<<"nsumt & data->LOS[d].TTFT "<<d+1<<"; "<<nsumt<<"; "<< data->LOS[d].TTFT<<endl;
        }

        if (fArrayTrigger)
        {
            // Trigger VFTX2 data
            uint32_t trigChannels = data->LOS[d].VTRIGF;
            for (uint32_t i = 0; i < trigChannels; i++)
            {
                uint32_t channelf = data->LOS[d].VTRIGFI[i]; // = 1..2
                uint32_t channelc = data->LOS[d].VTRIGCI[i]; // = 1..2
                if (channelf == channelc)
                    new ((*fArrayTrigger)[fArrayTrigger->GetEntriesFast()])
                        R3BLosMappedData(d + 1,                   // detector number
                                         channelf,                // channel number: 1
                                         0,                       // VFTX (0),TAMEX leading (1), TAMEX trailing (2)
                                         data->LOS[d].VTRIGFv[i], // VFTX fine time
                                         data->LOS[d].VTRIGCv[i]  // VFTX coarse time
                        );
            }

            // Trigger TAMEX3 data.
            trigChannels = data->LOS[d].TTRIGFL;
            for (uint32_t i = 0; i < trigChannels; i++)
            {
                // Leading times
                uint32_t channelfl = data->LOS[d].TTRIGFLI[i]; // = 1..2
                uint32_t channelcl = data->LOS[d].TTRIGCLI[i]; // = 1..2
                if (channelfl == channelcl)
                    new ((*fArrayTrigger)[fArrayTrigger->GetEntriesFast()])
                        R3BLosMappedData(d + 1,                    // detector number
                                         channelfl,                // channel number: 1
                                         1,                        // VFTX (0),TAMEX leading (1), TAMEX trailing (2)
                                         data->LOS[d].TTRIGFLv[i], // TAMEX3 fine time
                                         data->LOS[d].TTRIGCLv[i]  // TAMEX3 coarse time
                        );

                // Trailing times
                uint32_t channelft = data->LOS[d].TTRIGFTI[i]; // = 1..2
                uint32_t channelct = data->LOS[d].TTRIGCTI[i]; // = 1..2
                if (channelft == channelct)
                    new ((*fArrayTrigger)[fArrayTrigger->GetEntriesFast()])
                        R3BLosMappedData(d + 1,                    // detector number
                                         channelft,                // channel number: 1
                                         2,                        // VFTX (0),TAMEX leading (1), TAMEX trailing (2)
                                         data->LOS[d].TTRIGFTv[i], // TAMEX3 fine time
                                         data->LOS[d].TTRIGCTv[i]  // TAMEX3 coarse time
                        );
            }
        }
    }
    fNEvents += 1;

    return kTRUE;
}

void R3BLosReader::Reset()
{
    // Reset the output array
    fArray->Clear();
    if (fArrayTrigger)
        fArrayTrigger->Clear();
}

ClassImp(R3BLosReader)
