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

#include "R3BAmsMappedData.h"
#include "R3BCalifaJulichReader.h"
#include "R3BCalifaMappedData.h"
#include "R3BLogger.h"

#include <TClonesArray.h>
#include <ext_data_struct_info.hh>

/**
 ** ext_h101_califa.h was created by running
 ** $unpacker --ntuple=STRUCT_HH,RAW:CALIFA,id=h101_CALIFA,NOTRIGEVENTNO,ext_h101_califa.h
 **/

extern "C"
{
#include "ext_data_client.h"
#include "ext_h101_califa_julich.h"
}

R3BCalifaJulichReader::R3BCalifaJulichReader(EXT_STR_h101_CALIFA* data, size_t offset)
    : R3BReader("R3BCalifaJulichReader")
    , fNEvent(0)
    , fData(data)
    , fOffset(offset)
    , fOnline(kFALSE)
    , fArrayCalifa(new TClonesArray("R3BCalifaMappedData"))
    , fArrayAms(new TClonesArray("R3BAmsMappedData"))
{
}

R3BCalifaJulichReader::~R3BCalifaJulichReader()
{
    if (fArrayCalifa)
        delete fArrayCalifa;
    if (fArrayAms)
        delete fArrayAms;
}

Bool_t R3BCalifaJulichReader::Init(ext_data_struct_info* a_struct_info)
{
    Int_t ok;
    R3BLOG(info, "");
    EXT_STR_h101_CALIFA_ITEMS_INFO(ok, *a_struct_info, fOffset, EXT_STR_h101_CALIFA, 0);

    if (!ok)
    {
        R3BLOG(error, "Failed to setup structure information");
        return kFALSE;
    }

    // Register output array in tree
    FairRootManager::Instance()->Register("CalifaMappedData", "Califa", fArrayCalifa, !fOnline);
    fArrayCalifa->Clear();
    FairRootManager::Instance()->Register("AmsMappedData", "Si detector", fArrayAms, !fOnline);
    fArrayAms->Clear();

    return kTRUE;
}

Bool_t R3BCalifaJulichReader::R3BRead()
{
    R3BLOG(debug1, "Event data");
    // SELECT THE FOR LOOP BASED ON THE MAPPING...
    for (int crystal = 0; crystal < fData->CALIFA_ENE; ++crystal)
    {
        UShort_t channelNumber = fData->CALIFA_ENEI[crystal];
        int16_t energy = fData->CALIFA_ENEv[crystal];
        int16_t nf = fData->CALIFA_NFv[crystal];
        int16_t ns = fData->CALIFA_NSv[crystal];

        uint64_t febextime = ((uint64_t)fData->CALIFA_TSMSBv[crystal] << 32) | (uint64_t)fData->CALIFA_TSLSBv[crystal];

        uint64_t wrts = ((uint64_t)fData->CALIFA_WRTS_T4v[crystal] << 48) |
                        ((uint64_t)fData->CALIFA_WRTS_T3v[crystal] << 32) |
                        ((uint64_t)fData->CALIFA_WRTS_T2v[crystal] << 16) | (uint64_t)fData->CALIFA_WRTS_T1v[crystal];

        int32_t ov = fData->CALIFA_OVv[crystal];
        int16_t pu = fData->CALIFA_PILEUPv[crystal];
        int16_t dc = fData->CALIFA_DISCARDv[crystal];

        int16_t tot = fData->CALIFA_TOTv[crystal];

        if (channelNumber < 1 || channelNumber > 80)
            LOG(fatal) << "\033[5m\033[31m error TO CHECK!! Channel number=" << channelNumber
                       << " detected in data \033[0m";

        if (channelNumber < 17)
            new ((*fArrayCalifa)[fArrayCalifa->GetEntriesFast()])
                R3BCalifaMappedData(channelNumber, energy, nf, ns, febextime, wrts, ov, pu, dc, tot);
        else if (channelNumber < 81)
            new ((*fArrayAms)[fArrayAms->GetEntriesFast()]) R3BAmsMappedData(1, channelNumber - 16, energy);
    }
    fNEvent += 1;
    return kTRUE;
}

void R3BCalifaJulichReader::Reset()
{
    // Reset the output array
    fArrayCalifa->Clear();
    fArrayAms->Clear();
}

ClassImp(R3BCalifaJulichReader)
