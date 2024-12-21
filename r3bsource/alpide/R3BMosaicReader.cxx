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

#include "R3BAlpideMappedData.h"
#include "R3BLogger.h"
#include "R3BMosaicReader.h"

#include <TClonesArray.h>
#include <ext_data_struct_info.hh>

/**
 ** ext_h101_mosaic.h was created by running
 ** $unpacker --ntuple=STRUCT_HH,RAW:MOSAIC,id=h101_MOSAIC,NOTRIGEVENTNO,ext_h101_mosaic.h
 **/

extern "C"
{
#include "ext_data_client.h"
#include "ext_h101_mosaic.h"
}

R3BMosaicReader::R3BMosaicReader(EXT_STR_h101_MOSAIC_onion* data, size_t offset)
    : R3BReader("R3BMosaicReader")
    , fData(data)
    , fNbMosaic(sizeof(fData->MOSAIC) / sizeof(fData->MOSAIC[0]))
    , fOffset(offset)
    , fArray(new TClonesArray("R3BAlpideMappedData"))
{
}

R3BMosaicReader::~R3BMosaicReader()
{
    if (fArray)
    {
        delete fArray;
    }
}

Bool_t R3BMosaicReader::Init(ext_data_struct_info* a_struct_info)
{
    int okay = 0;
    R3BLOG(info, "");
    EXT_STR_h101_MOSAIC_ITEMS_INFO(okay, *a_struct_info, fOffset, EXT_STR_h101_MOSAIC, 0);

    R3BLOG_IF(fatal, !okay, "Failed to setup structure information.");

    // Register output array in tree
    FairRootManager::Instance()->Register("AlpideMappedData", "ALPIDE_Map", fArray, !fOnline);
    Reset();
    memset(fData, 0, sizeof(*fData));

    return kTRUE;
}

Bool_t R3BMosaicReader::R3BRead()
{
    R3BLOG(debug1, "Event data: " << fNEvent);
    fNEvent++;

    for (int mosid = 0; mosid < fNbMosaic; mosid++)
    {
        for (int hits = 0; hits < fData->MOSAIC[mosid].CHIP; hits++)
        {
            int fChipId = fData->MOSAIC[mosid].CHIPv[hits];
            // R3BLOG_IF(error, fChipId < 1, "Wrong fChipId: " << fChipId);
            // if (fChipId == 0)
            //    continue;

            int fAlpideId = map_mosaics[mosid] * 6 + fChipId + 1; // 1-base: 1 --> 24 (6*4)
            R3BLOG_IF(error, fAlpideId < 1, "Wrong fAlpideId: " << fAlpideId);

            fAlpideId = map_sensors[fAlpideId - 1];
            new ((*fArray)[fArray->GetEntriesFast()]) R3BAlpideMappedData(
                fAlpideId, 0, mosid + 1, fChipId, fData->MOSAIC[mosid].ROWv[hits], fData->MOSAIC[mosid].COLv[hits]);
        }
    }
    return kTRUE;
}

void R3BMosaicReader::Reset()
{
    // Reset the output array
    if (fArray)
    {
        fArray->Clear();
    }
}

ClassImp(R3BMosaicReader)
