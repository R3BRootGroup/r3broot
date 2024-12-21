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

#include "R3BMusicReader.h"
#include "R3BMusicMappedData.h"

#include <FairLogger.h>
#include <FairRootManager.h>

#include <TClonesArray.h>
#include <iostream>

extern "C"
{
#include "ext_data_client.h"
#include "ext_h101_music.h"
}

R3BMusicReader::R3BMusicReader(EXT_STR_h101_MUSIC* data, size_t offset)
    : R3BReader("R3BMusicReader")
    , fData(data)
    , fOffset(offset)
    , fOnline(kFALSE)
    , fArray(new TClonesArray("R3BMusicMappedData"))
{
}

R3BMusicReader::~R3BMusicReader()
{
    if (fArray)
    {
        delete fArray;
    }
}

Bool_t R3BMusicReader::Init(ext_data_struct_info* a_struct_info)
{
    Int_t ok;
    LOG(info) << "R3BMusicReader::Init()";
    EXT_STR_h101_MUSIC_ITEMS_INFO(ok, *a_struct_info, fOffset, EXT_STR_h101_MUSIC, 0);
    if (!ok)
    {
        perror("ext_data_struct_info_item");
        LOG(error) << "Failed to setup structure information.";
        return kFALSE;
    }

    // Register output array in tree
    FairRootManager::Instance()->Register("MusicMappedData", "MUSIC", fArray, !fOnline);
    Reset();

    // clear struct_writer's output struct. Seems ucesb doesn't do that
    // for channels that are unknown to the current ucesb config.
    fData->MUSIC_EM = 0;
    fData->MUSIC_TM = 0;
    fData->MUSIC_TREFM = 0;
    for (int mult = 0; mult < NUM_MUSIC_ANODES; mult++)
        multPerAnode[mult] = 0;

    return kTRUE;
}

Bool_t R3BMusicReader::R3BRead()
{
    // Read data for different dates, if needed!
    ReadData201911();
    return kTRUE;
}

void R3BMusicReader::Reset()
{
    // Reset the output array
    fArray->Clear();
}

Bool_t R3BMusicReader::ReadData201911()
{
    // --- ---------------- --- //
    // --- TREF : TIME ONLY --- //
    // --- ---------------- --- //

    // mail from R. Schneider from June 16, 2016: "if you have several hits of one
    // channel in the window of interest, you get an energy and a timing for each
    // hit."
    // --> multiplicity per anode (for anode 1 to anode 8) should also be the same
    // in energy and time mail from R. Schneider from May 21st 2019 : "the hits
    // from one channel are kept in the chronological order."
    // --> for one anode with multi-hit, the first hit in energy correspond to the
    // first hit in time

    // return the number of Tref
    // good event : nTref=1 (one channel only, which can have multi hits)
    // uint32_t nTref = fData->MUSIC_TREF;
    // if (nTref > 1)
    //  LOG(warning) << "R3BMusicReader::ReadData error ! ONLY ONE Tref SIGNAL BUT
    //  HERE MULTIPLE SIGNALS IS FOUND!";

    uint32_t curTref = 0;
    uint32_t nextTref = 0;
    UShort_t idAnodeTref = 0;
    UShort_t nAnodesTref = fData->MUSIC_TREFM;

    for (UShort_t a = 0; a < nAnodesTref; a++)
    {
        idAnodeTref = fData->MUSIC_TREFMI[a] + 7;
        nextTref = fData->MUSIC_TREFME[a];
        multPerAnode[idAnodeTref] = nextTref - curTref;
        for (int hit = curTref; hit < nextTref; hit++)
            new ((*fArray)[fArray->GetEntriesFast()]) R3BMusicMappedData(idAnodeTref, fData->MUSIC_TREFv[hit], 0);
        curTref = nextTref;
    }

    // --- ---------------------- --- //
    // --- ANODES : ENERGY + TIME --- //
    // --- ---------------------- --- //
    // 1<=nAnodesEnergy<=8
    // 1<=nAnodesTime<=8
    UShort_t nAnodesEnergy = fData->MUSIC_EM;
    UShort_t nAnodesTime = fData->MUSIC_TM;

    // mail from R. Schneider from June 16, 2016: "if you have several hits of one
    // channel in the window of interest, you get an energy and a timing for each
    // hit."
    // --> multiplicity per anode (for anode 1 to anode 8) should also be the same
    // in energy and time mail from R. Schneider from May 21st 2019 : "the hits
    // from one channel are kept in the chronological order."
    // --> for one anode with multi-hit, the first hit in energy correspond to the
    // first hit in time
    if (nAnodesEnergy != nAnodesTime)
        LOG(warning) << "R3BMusicReader::ReadData error ! NOT THE SAME NUMBER OF "
                        "ANODES HITTED IN ENERGY AND TIME";

    // ENERGY AND TIME ARE SORTED
    uint32_t curAnodeTimeStart = 0;
    uint32_t curAnodeEnergyStart = 0;
    for (UShort_t a = 0; a < nAnodesTime; a++)
    {
        // EMI and TMI give the 1-based anode number
        UShort_t idAnodeTime = fData->MUSIC_TMI[a] - 1;
        UShort_t idAnodeEnergy = fData->MUSIC_EMI[a] - 1;
        if (idAnodeEnergy != idAnodeTime)
            LOG(error) << "R3BMusicReader::ReadData error ! MISMATCH FOR ANODE ID IN "
                          "ENERGY #"
                       << idAnodeEnergy << " AND TIME #" << idAnodeTime;
        uint32_t nextAnodeTimeStart = fData->MUSIC_TME[a];
        uint32_t nextAnodeEnergyStart = fData->MUSIC_EME[a];
        multPerAnode[idAnodeTime] = nextAnodeTimeStart - curAnodeTimeStart;
        if (multPerAnode[idAnodeTime] != (nextAnodeEnergyStart - curAnodeEnergyStart))
            LOG(error) << "R3BMusicReader::ReadData error ! MISMATCH FOR "
                          "MULTIPLICITY PER ANODE IN ENERGY AND TIME";
        for (int hit = curAnodeTimeStart; hit < nextAnodeTimeStart; hit++)
        { // idAnodeEnergy from 1 to 8, so -1 to get from 0 to 7
            new ((*fArray)[fArray->GetEntriesFast()])
                R3BMusicMappedData(idAnodeEnergy, fData->MUSIC_Tv[hit], fData->MUSIC_Ev[hit]);
        }
        curAnodeEnergyStart = nextAnodeEnergyStart;
        curAnodeTimeStart = nextAnodeTimeStart;
    } // end of loop over the anodes from 1 to 8

    return kTRUE;
}

ClassImp(R3BMusicReader)
