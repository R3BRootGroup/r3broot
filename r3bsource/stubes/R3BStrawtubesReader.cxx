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
#include "R3BStrawtubesMappedData.h"
#include "R3BStrawtubesReader.h"
#include "TClonesArray.h"

extern "C"
{
#include "ext_data_client.h"
#include "ext_h101_strawtubes.h"
}

#define LENGTH(x) (sizeof x / sizeof *x)

R3BStrawtubesReader::R3BStrawtubesReader(EXT_STR_h101_STRAWTUBES* data, UInt_t offset)
    : R3BReader("R3BStrawtubesReader")
    , fOnion(reinterpret_cast<EXT_STR_h101_STRAWTUBES_onion*>(data))
    , fOffset(offset)
    , fLogger(FairLogger::GetLogger())
    , fArray(new TClonesArray("R3BStrawtubesMappedData"))
{
}

R3BStrawtubesReader::~R3BStrawtubesReader() {}

Bool_t R3BStrawtubesReader::Init(ext_data_struct_info* a_struct_info)
{
    int ok;
    EXT_STR_h101_STRAWTUBES_ITEMS_INFO(ok, *a_struct_info, fOffset, EXT_STR_h101_STRAWTUBES, 1);
    if (!ok)
    {
        perror("ext_data_struct_info_item");
        LOG(error) << "Failed to setup structure information.";
        return kFALSE;
    }

    FairRootManager::Instance()->Register("StrawtubesMapped", "Land", fArray, kTRUE);

    memset(fOnion, 0, sizeof *fOnion);

    return kTRUE;
}

Bool_t R3BStrawtubesReader::R3BRead()
{
    for (int plane_i = 0; LENGTH(fOnion->STRAW) > plane_i; ++plane_i)
    {
        for (int side_i = 0; LENGTH(fOnion->STRAW[plane_i].S) > side_i; ++side_i)
        {
            uint32_t vector_i = 0;
            uint32_t const c_channel_nonzero_num = fOnion->STRAW[plane_i].S[side_i].TFM;
            for (uint32_t i = 0; c_channel_nonzero_num > i; ++i)
            {
                // UCESB channels are 1-based!
                uint32_t channel_i = fOnion->STRAW[plane_i].S[side_i].TFMI[i];
                uint32_t channel_end = fOnion->STRAW[plane_i].S[side_i].TFME[i];
                for (; channel_end > vector_i; ++vector_i)
                {
                    new ((*fArray)[fArray->GetEntriesFast()])
                        R3BStrawtubesMappedData(1 + plane_i,
                                                channel_i,
                                                1 + side_i,
                                                fOnion->STRAW[plane_i].S[side_i].TCv[vector_i],
                                                fOnion->STRAW[plane_i].S[side_i].TFv[vector_i]);
                }
            }
            if (fOnion->STRAW[plane_i].S[side_i].TF != vector_i)
            {
                LOG(error) << "Straw plane=" << plane_i
                           << ": UCESB multi-hit "
                              "structure corrupt (TF="
                           << fOnion->STRAW[plane_i].S[side_i].TF << ", vector_i=" << vector_i << ")!";
            }
        }
    }
    return kTRUE;
}

void R3BStrawtubesReader::Reset() { fArray->Clear(); }

ClassImp(R3BStrawtubesReader)
