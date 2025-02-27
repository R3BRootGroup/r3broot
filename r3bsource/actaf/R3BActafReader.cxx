/******************************************************************************
 *   Copyright (C) 2025 GSI Helmholtzzentrum für Schwerionenforschung GmbH    *
 *   Copyright (C) 2025 Members of R3B Collaboration                          *
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

#include "R3BActafMappedData.h"
#include "R3BActafReader.h"
#include "R3BLogger.h"

#include <TClonesArray.h>
#include <iostream>
#include <vector>

constexpr int eChn = 16;
constexpr int Bins = 2692;
constexpr int Modules = 4;

/**
 ** ext_h101_actaf.h was created by running
 ** $unpacker --ntuple=STRUCT_HH,RAW:ACTAF,id=h101_ACTAF,NOTRIGEVENTNO,ext_h101_actaf.h
 **/

extern "C"
{
#include "ext_data_client.h"
#include "ext_h101_actaf.h"
}

R3BActafReader::R3BActafReader(EXT_STR_h101_ACTAF_onion* data, size_t offset)
    : R3BReader("R3BActafReader")
    , fData(data)
    , fOffset(offset)
    , fArray(std::make_unique<TClonesArray>("R3BActafMappedData"))
{
}

Bool_t R3BActafReader::Init(ext_data_struct_info* a_struct_info)
{
    R3BLOG(info, "");
    Int_t okey = 0;
    constexpr int esf = 1;
    EXT_STR_h101_ACTAF_ITEMS_INFO(okey, *a_struct_info, fOffset, EXT_STR_h101_ACTAF, esf);
    if (!okey)
    {
        R3BLOG(error, "Failed to setup structure information");
        return kFALSE;
    }

    // Register output array in tree
    FairRootManager::Instance()->Register("ActafMappedData", "Actaf mapped data", fArray.get(), !fOnline);
    Reset();
    memset(fData, 0, sizeof *fData);

    return kTRUE;
}

Bool_t R3BActafReader::R3BRead()
{
    R3BLOG(debug1, "Event data.");

    for (int mod = 0; mod < Modules; ++mod)
    {
        std::vector<std::vector<UInt_t>> trace(eChn, std::vector<UInt_t>(Bins));

        for (int chn = 0; chn < fData->ACTAF[mod].TRACERAW; ++chn)
        {
            int row = chn / Bins;
            int col = chn % Bins;
            trace[row][col] = fData->ACTAF[mod].TRACERAWv[chn];
        }

        for (int chn = 0; chn < fData->ACTAF[mod].CH; ++chn)
        {
            new ((*fArray)[fArray->GetEntriesFast()]) R3BActafMappedData(mod * eChn + (chn + 1), trace[chn]);
        }
    }
    fNEvent += 1;
    return kTRUE;
}

void R3BActafReader::Reset()
{
    // Reset the output array
    fArray->Clear();
}

ClassImp(R3BActafReader)
