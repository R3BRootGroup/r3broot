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

#include "R3BFootMappedData.h"
#include "R3BFootSiReader.h"
#include "R3BLogger.h"

#include "TClonesArray.h"
#include "ext_data_struct_info.hh"

/**
 ** ext_h101_foot.h was created by running
 ** $unpacker --ntuple=STRUCT_HH,RAW:FOOT,id=h101_FOOT,NOTRIGEVENTNO,ext_h101_foot.h
 **/

extern "C"
{
#include "ext_data_client.h"
#include "ext_h101_foot.h"
}

R3BFootSiReader::R3BFootSiReader(EXT_STR_h101_FOOT_onion* data, size_t offset)
    : R3BReader("R3BFootSiReader")
    , fNEvent(0)
    , fData(data)
    , fOffset(offset)
    , fOnline(kFALSE)
    , fNbDet(16)
    //, fNbDet(sizeof(EXT_STR_h101_FOOT_onion) / sizeof(EXT_STR_h101_FOOT_onion.FOOT[0])) // Auto-gets # FEET from
    // struct!
    , fArray(new TClonesArray("R3BFootMappedData"))
{
}

R3BFootSiReader::~R3BFootSiReader() { delete fArray; }

Bool_t R3BFootSiReader::Init(ext_data_struct_info* a_struct_info)
{
    Int_t ok;
    R3BLOG(info, "");
    EXT_STR_h101_FOOT_ITEMS_INFO(ok, *a_struct_info, fOffset, EXT_STR_h101_FOOT, 0);

    if (!ok)
    {
        R3BLOG(error, "Failed to setup structure information");
        return kFALSE;
    }

    // Register output array in tree
    FairRootManager::Instance()->Register("FootMappedData", "Foot mapped data", fArray, !fOnline);
    fArray->Clear();

    memset(fData, 0, sizeof *fData);

    return kTRUE;
}

Bool_t R3BFootSiReader::R3BRead()
{
    R3BLOG(debug1, "Event data");
    // Read FOOT detectors
    for (Int_t d = 0; d < fNbDet; d++)
    {
        if (fData->FOOT[d]._ == 640)
        {
            for (Int_t strip = 0; strip < fData->FOOT[d]._; ++strip)
            {
                new ((*fArray)[fArray->GetEntriesFast()]) R3BFootMappedData(d + 1, strip + 1, fData->FOOT[d].E[strip]);
            }
        }
        else if (fData->FOOT[d]._ == 0)
        {
        }
        else
        {
            if (fNEvent > 0)
                LOG(fatal) << "\033[5m\033[31m R3BFootSiReader::Failed number of strips per detector. \033[0m";
        }
    }
    fNEvent += 1;
    return kTRUE;
}

void R3BFootSiReader::Reset()
{
    // Reset the output array
    fArray->Clear();
}

ClassImp(R3BFootSiReader);
