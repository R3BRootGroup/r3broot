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

#include "R3BAlpideMappedData.h"
#include "R3BAlpideReader.h"
#include "R3BLogger.h"

#include "TClonesArray.h"
#include "ext_data_struct_info.hh"

/**
 ** ext_h101_alpide.h was created by running
 ** $unpacker --ntuple=STRUCT_HH,RAW:ALPIDE,id=h101_ALPIDE,NOTRIGEVENTNO,ext_h101_alpide.h
 **/

extern "C"
{
#include "ext_data_client.h"
#include "ext_h101_alpide.h"
}

R3BAlpideReader::R3BAlpideReader(EXT_STR_h101_ALPIDE_onion* data, size_t offset)
    : R3BReader("R3BAlpideReader")
    , fNEvent(1)
    , fData(data)
    , fNbDet(sizeof(fData->ALPIDE) / sizeof(fData->ALPIDE[0]))
    , fOffset(offset)
    , fOnline(kFALSE)
    , fArray(new TClonesArray("R3BAlpideMappedData"))
{
}

R3BAlpideReader::~R3BAlpideReader()
{
    R3BLOG(debug1, "");
    if (fArray)
    {
        delete fArray;
    }
}

Bool_t R3BAlpideReader::Init(ext_data_struct_info* a_struct_info)
{
    Int_t ok;
    R3BLOG(info, "");
    EXT_STR_h101_ALPIDE_ITEMS_INFO(ok, *a_struct_info, fOffset, EXT_STR_h101_ALPIDE, 0);

    R3BLOG_IF(fatal, !ok, "Failed to setup structure information.");

    // Register output array in tree
    FairRootManager::Instance()->Register("AlpideMappedData", "ALPIDE_Map", fArray, !fOnline);
    Reset();
    memset(fData, 0, sizeof(*fData));

    return kTRUE;
}

Bool_t R3BAlpideReader::R3BRead()
{
    R3BLOG(debug1, "Event data: " << fNEvent);

    for (int d = 0; d < fNbDet; d++)
    {
        R3BLOG_IF(error,
                  fData->ALPIDE[d].REGION != fData->ALPIDE[d].ADDRESS,
                  "Region/Address sizes mismatch for detector " << d + 1 << ", Region: " << fData->ALPIDE[d].REGION
                                                                << " , Address: " << fData->ALPIDE[d].ADDRESS);

        R3BLOG_IF(error,
                  fData->ALPIDE[d].ROW != fData->ALPIDE[d].COL,
                  "Row/Col sizes mismatch for detector " << d + 1 << ", Row: " << fData->ALPIDE[d].ROW
                                                         << " , Col: " << fData->ALPIDE[d].COL);

        R3BLOG_IF(error,
                  fData->ALPIDE[d].CHIP != fData->ALPIDE[d].ROW,
                  "Chip/Row sizes mismatch for detector " << d + 1 << ", Chip: " << fData->ALPIDE[d].CHIP
                                                          << " , Row: " << fData->ALPIDE[d].ROW);

        R3BLOG_IF(error,
                  fData->ALPIDE[d].CHIP != fData->ALPIDE[d].COL,
                  "Chip/Col sizes mismatch for detector " << d + 1 << ", Chip: " << fData->ALPIDE[d].CHIP
                                                          << " , Col: " << fData->ALPIDE[d].COL);

        if (fData->ALPIDE[d].COL != fData->ALPIDE[d].ROW)
        {
            continue;
        }

        for (int r = 0; r < fData->ALPIDE[d].ROW; r++)
        {

            //  R3BLOG(debug1,"det: "<<d+1 <<", region: "<<fData->ALPIDE[d].REGIONv[fData->ALPIDE[d].REGION-1]<< " ,
            //  ads: " << fData->ALPIDE[d].ADDRESSv[fData->ALPIDE[d].ADDRESS-1]<< " , row: "<<
            //  fData->ALPIDE[d].ROWv[r]<<" , col: "<< fData->ALPIDE[d].COLv[r]);

            if (fData->ALPIDE[d].ROWv[r] > 0 && fData->ALPIDE[d].COLv[r] > 0 && fData->ALPIDE[d].COLv[r] < 1025 &&
                fData->ALPIDE[d].ROWv[r] < 513)
            {
                new ((*fArray)[fArray->GetEntriesFast()]) R3BAlpideMappedData(d + 1,
                                                                              fData->ALPIDE[d].REGIONv[r],
                                                                              fData->ALPIDE[d].ADDRESSv[r],
                                                                              fData->ALPIDE[d].CHIPv[r],
                                                                              fData->ALPIDE[d].ROWv[r],
                                                                              fData->ALPIDE[d].COLv[r]);
            }
        }
    }

    return kTRUE;
}

void R3BAlpideReader::Reset()
{
    // Reset the output array
    if (fArray)
    {
        fArray->Clear();
    }
}

ClassImp(R3BAlpideReader);
