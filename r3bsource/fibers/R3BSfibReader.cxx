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

#include "R3BSfibReader.h"
#include "FairLogger.h"
#include "FairRootManager.h"
#include "R3BEventHeader.h"
#include "R3BSfibMappedData.h"
#include "TClonesArray.h"

extern "C"
{
#include "ext_data_client.h"
#include "ext_h101_sfib.h"
}
#include <iostream>

using namespace std;

R3BSfibReader::R3BSfibReader(EXT_STR_h101_sfib* data, UInt_t offset)
    : R3BReader("R3BSfibReader")
    , fData(data)
    , fOffset(offset)
    , fLogger(FairLogger::GetLogger())
    , fMappedArray(new TClonesArray("R3BSfibMappedData"))
{
}

Bool_t R3BSfibReader::Init(ext_data_struct_info* a_struct_info)
{
    int ok;

    EXT_STR_h101_sfib_ITEMS_INFO(ok, *a_struct_info, fOffset, EXT_STR_h101_sfib, 0);

    if (!ok)
    {
        perror("ext_data_struct_info_item");
        LOG(error) << "Failed to setup structure information.";
        return kFALSE;
    }

    FairRootManager* mgr = FairRootManager::Instance();

    FairRootManager::Instance()->Register("SfibMapped", "Land", fMappedArray, kTRUE);

    return kTRUE;
}

Bool_t R3BSfibReader::R3BRead()
{
#define PERM(side, side_i, edge, edge_i)                                   \
    do                                                                     \
    {                                                                      \
        uint32_t c_M = fData->SFIB_##side##T##edge##CM;                    \
        uint32_t cur_entry = 0;                                            \
        for (uint32_t i = 0; i < c_M; i++)                                 \
        {                                                                  \
            uint32_t c_MI = fData->SFIB_##side##T##edge##CMI[i];           \
            uint32_t c_ME = fData->SFIB_##side##T##edge##CME[i];           \
            for (; cur_entry < c_ME; cur_entry++)                          \
            {                                                              \
                new ((*fMappedArray)[fMappedArray->GetEntriesFast()])      \
                    R3BSfibMappedData(side_i,                              \
                                      c_MI,                                \
                                      0 == edge_i,                         \
                                      fData->SFIB_##side##T##edge##Cv[i],  \
                                      fData->SFIB_##side##T##edge##Fv[i]); \
            }                                                              \
        }                                                                  \
    } while (0)
    PERM(B, 0, L, 0);
    PERM(B, 0, T, 1);
    PERM(T, 1, L, 0);
    PERM(T, 1, T, 1);

    return kTRUE;
}

void R3BSfibReader::Reset() { fMappedArray->Clear(); }

ClassImp(R3BSfibReader)
