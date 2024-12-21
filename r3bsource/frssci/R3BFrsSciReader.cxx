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

#include "R3BFrsSciReader.h"
#include "R3BFrsSciMappedData.h"

#include <FairLogger.h>
#include <FairRootManager.h>

extern "C"
{
#include "ext_data_client.h"
#include "ext_h101_frssci.h"
}

#include <iostream>

using namespace std;

R3BFrsSciReader::R3BFrsSciReader(EXT_STR_h101_FRSSCI* data, size_t offset)
    : R3BFrsSciReader(data, offset, 1)
{
}

R3BFrsSciReader::R3BFrsSciReader(EXT_STR_h101_FRSSCI* data, size_t offset, UShort_t NumSci)
    : R3BReader("R3BFrsSciReader")
    , fData(data)
    , fOffset(offset)
    , fOnline(kFALSE)
    , fArray(new TClonesArray("R3BFrsSciMappedData")) // class name
    , fNumSci(NumSci)
{
}

R3BFrsSciReader::~R3BFrsSciReader()
{
    LOG(debug) << "R3BFrsSciReader: Delete instance";
    if (fArray)
    {
        delete fArray;
    }
}

Bool_t R3BFrsSciReader::Init(ext_data_struct_info* a_struct_info)
{
    Int_t ok;
    LOG(info) << "R3BFrsSciReader::Init()";
    EXT_STR_h101_FRSSCI_ITEMS_INFO(ok, *a_struct_info, fOffset, EXT_STR_h101_FRSSCI, 0);
    if (!ok)
    {
        perror("ext_data_struct_info_item");
        LOG(error) << "R3BFrsSciReader::Failed to setup structure information.";
        return kFALSE;
    }

    // Register output array in tree
    FairRootManager::Instance()->Register("FrsSciMappedData", "FrsSci", fArray, !fOnline);
    fArray->Clear();

    // clear struct_writer's output struct. Seems ucesb doesn't do that
    // for channels that are unknown to the current ucesb config.
    EXT_STR_h101_FRSSCI_onion* data = (EXT_STR_h101_FRSSCI_onion*)fData;
    for (UShort_t d = 0; d < fNumSci; d++)
        data->FRSSCI[d].VTFM = 0;

    return kTRUE;
}

Bool_t R3BFrsSciReader::R3BRead()
{
    // Convert plain raw data to multi-dimensional array
    EXT_STR_h101_FRSSCI_onion* data = (EXT_STR_h101_FRSSCI_onion*)fData;

    // loop over all detectors
    for (UShort_t d = 0; d < fNumSci; d++)
    {
        uint32_t numberOfPMTsWithHits_TC = (uint32_t)data->FRSSCI[d].VTCM;
        uint32_t numberOfPMTsWithHits_TF = (uint32_t)data->FRSSCI[d].VTFM;
        if (numberOfPMTsWithHits_TF != numberOfPMTsWithHits_TC)
        {
            LOG(error) << "R3BFrsSciReader::R3BRead() Error in unpacking, unconsistency between TF and TC for FrsSci !";
        }
        else
        {
            // loop over channels with hits
            uint32_t curChannelStart = 0;
            for (Int_t pmmult = 0; pmmult < numberOfPMTsWithHits_TF; pmmult++)
            {
                uint32_t pmtid_TC = (uint32_t)data->FRSSCI[d].VTCMI[pmmult];
                uint32_t pmtid_TF = (uint32_t)data->FRSSCI[d].VTFMI[pmmult];
                if (pmtid_TF != pmtid_TC)
                {
                    LOG(error)
                        << "R3BFrsSciReader::R3BRead() Error in unpacking, unconsistency between the PMT id for TF "
                           "and TC for FrsSci !";
                }
                uint32_t nextChannelStart = (uint32_t)data->FRSSCI[d].VTFME[pmmult];
                // put the mapped items {det,pmt,finetime, coarsetime} one after the other in the fArray
                for (Int_t hit = curChannelStart; hit < nextChannelStart; hit++)
                {
                    new ((*fArray)[fArray->GetEntriesFast()]) R3BFrsSciMappedData(d + 1, // 1-based numbering
                                                                                  pmtid_TF,
                                                                                  data->FRSSCI[d].VTCv[hit],
                                                                                  data->FRSSCI[d].VTFv[hit]);
                }
                curChannelStart = nextChannelStart;
            }
        }
    } // end of for(d)
    return kTRUE;
}

void R3BFrsSciReader::Reset()
{
    // Reset the output array
    fArray->Clear();
}

ClassImp(R3BFrsSciReader)
