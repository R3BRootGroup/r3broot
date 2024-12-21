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
#include "R3BMusliMappedData.h"
#include "R3BMusliReader.h"

extern "C"
{
#include "ext_data_client.h"
#include "ext_h101_musli.h"
}

using namespace std;

R3BMusliReader::R3BMusliReader(EXT_STR_h101_MUSLI* data, size_t offset)
    : R3BReader("R3BMusliReader")
    , fData(data)
    , fOffset(offset)
    , fOnline(kFALSE)
    , fArray(new TClonesArray("R3BMusliMappedData"))
{
}

R3BMusliReader::~R3BMusliReader()
{
    R3BLOG(debug1, "Destructor");
    if (fArray)
    {
        delete fArray;
    }
}

Bool_t R3BMusliReader::Init(ext_data_struct_info* a_struct_info)
{
    Int_t ok;
    R3BLOG(info, "");
    EXT_STR_h101_MUSLI_ITEMS_INFO(ok, *a_struct_info, fOffset, EXT_STR_h101_MUSLI, 0);
    if (!ok)
    {
        R3BLOG(error, "Failed to setup structure information.");
        return kFALSE;
    }

    // Register output array in tree
    FairRootManager::Instance()->Register("MusliMappedData", "MappedMusli", fArray, !fOnline);

    Reset();
    memset(fData, 0, sizeof *fData);

    // clear struct_writer's output struct. Seems ucesb doesn't do that
    // for channels that are unknown to the current ucesb config.
    auto* data = reinterpret_cast<EXT_STR_h101_MUSLI_onion*>(fData);
    data->MUSLI_EM = 0;
    data->MUSLI_TM = 0;
    data->MUSLI_TREFM = 0;
    data->MUSLI_TTRIGM = 0;

    return kTRUE;
}

Bool_t R3BMusliReader::R3BRead()
{
    // Convert plain raw data to multi-dimensional array
    auto* data = reinterpret_cast<EXT_STR_h101_MUSLI_onion*>(fData);

    ReadData(data);
    return kTRUE;
}

void R3BMusliReader::Reset()
{
    // Reset the output array
    if (fArray)
        fArray->Clear();
}

Bool_t R3BMusliReader::ReadData(EXT_STR_h101_MUSLI_onion* data)
{
    Bool_t pileupFLAG;
    Bool_t overflowFLAG;

    // --- ---------------------------- --- //
    // --- NUMBER OF CHANNELS WITH HITS --- //
    // --- ---------------------------- --- //
    UShort_t nTref = data->MUSLI_TREFM;   // should always be 1
    UShort_t nTtrig = data->MUSLI_TTRIGM; // should always be 0 or 1
    UShort_t nEnergy = data->MUSLI_EM;
    UShort_t nTime = data->MUSLI_TM;

    // --- ----------------- --- //
    // --- MUSLI MAPPED DATA --- //
    // --- ----------------- --- //

    // --- TREF --- //
    uint32_t curTref = 0;
    uint32_t nextTref = 0;

    for (UShort_t a = 0; a < nTref; a++)
    {
        // TREFMI gives the 1-based Tref number
        nextTref = data->MUSLI_TREFME[a];
        for (int hit = curTref; hit < nextTref; hit++)
        {
            pileupFLAG = (data->MUSLI_TREFv[hit] & 0x00040000) >> 18;
            overflowFLAG = (data->MUSLI_TREFv[hit] & 0x00080000) >> 19;
            new ((*fArray)[fArray->GetEntriesFast()])
                R3BMusliMappedData(17, 1, 17, data->MUSLI_TREFv[hit], 0, pileupFLAG, overflowFLAG);
        }
        curTref = nextTref;
    }

    // --- TTRIG --- //
    uint32_t curTtrig = 0;
    uint32_t nextTtrig = 0;

    for (UShort_t a = 0; a < nTtrig; a++)
    {
        // again TTRIGMI is 1 based
        nextTtrig = data->MUSLI_TTRIGME[a];
        for (int hit = curTtrig; hit < nextTtrig; hit++)
        {
            pileupFLAG = (data->MUSLI_TTRIGv[hit] & 0x00040000) >> 18;
            overflowFLAG = (data->MUSLI_TTRIGv[hit] & 0x00080000) >> 19;
            new ((*fArray)[fArray->GetEntriesFast()])
                R3BMusliMappedData(17, 2, 18, data->MUSLI_TTRIGv[hit], 0, pileupFLAG, overflowFLAG);
        }
        curTtrig = nextTtrig;
    }

    // --- ENERGY AND TIME SIGNALS --- //
    if (nEnergy != nTime)
        LOG(error) << "R3BMusliReader::ReadData error ! NOT THE SAME NUMBER OF ANODES HITTED IN ENERGY () AND TIME ()";

    // ENERGY AND TIME ARE SORTED
    uint32_t curAnodeTimeStart = 0;
    uint32_t curAnodeEnergyStart = 0;
    for (UShort_t a = 0; a < nTime; a++)
    {
        // EMI and TMI give the 1-based signal number
        UShort_t idTime = data->MUSLI_TMI[a];
        UShort_t idEnergy = data->MUSLI_EMI[a];
        if (idTime != idEnergy)
            LOG(error) << "R3BMusliReader::ReadData error ! Mismatch between numbering of Time and Energy signals";

        uint32_t nextAnodeTimeStart = data->MUSLI_TME[a];
        uint32_t nextAnodeEnergyStart = data->MUSLI_EME[a];
        uint32_t index, type;
        for (int hit = curAnodeTimeStart; hit < nextAnodeTimeStart; hit++)
        {
            pileupFLAG = (data->MUSLI_Ev[hit] & 0x00040000) >> 18;
            overflowFLAG = (data->MUSLI_Ev[hit] & 0x00080000) >> 19;
            type = 2 * 2 * pow(2, uint32_t(idTime / 9) + uint32_t(idTime / 13) + uint32_t(idTime / 15));
            index = idTime - uint32_t(idTime / 9) * 8 - uint32_t(idTime / 13) * 4 - uint32_t(idTime / 15) * 2;
            new ((*fArray)[fArray->GetEntriesFast()]) R3BMusliMappedData(
                type, index, idTime, data->MUSLI_Tv[hit], data->MUSLI_Ev[hit], pileupFLAG, overflowFLAG);
        }
        curAnodeEnergyStart = nextAnodeEnergyStart;
        curAnodeTimeStart = nextAnodeTimeStart;
    } // end of loop over the e and t signals from 1 to 15
    return kTRUE;
}

ClassImp(R3BMusliReader)
