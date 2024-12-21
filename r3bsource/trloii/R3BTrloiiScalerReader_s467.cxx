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

#include "R3BTrloiiData.h"
#include "R3BTrloiiScalerReader_s467.h"

#include "TClonesArray.h"
#include "ext_data_struct_info.hh"

/**
 ** ext_h101_trlo.h was created by running
 ** ./202002_s467
 *--ntuple=STRUCT_HH,RAW:TRLORAW,RAW:TRLOADT,RAW:TRLOBDT,RAW:TRLOARD,id=h101_TRLO,NOTRIGEVENTNO,ext_h101_trlo_s467.h
 **/

extern "C"
{
#include "ext_data_client.h"
#include "ext_h101_trlo_s467.h"
}

R3BTrloiiScalerReader_s467::R3BTrloiiScalerReader_s467(EXT_STR_h101_TRLO_s467_onion* data, size_t offset)
    : R3BReader("R3BTrloiiScalerReader_s467")
    , fNEvent(0)
    , fData(data)
    , fOffset(offset)
    , fOnline(kFALSE)
    , fArray(new TClonesArray("R3BTrloiiData"))
{
}

R3BTrloiiScalerReader_s467::~R3BTrloiiScalerReader_s467()
{
    if (fArray)
        delete fArray;
}

Bool_t R3BTrloiiScalerReader_s467::Init(ext_data_struct_info* a_struct_info)
{
    Int_t ok;
    LOG(info) << "R3BTrloiiScalerReader_s467::Init()";
    EXT_STR_h101_TRLO_s467_ITEMS_INFO(ok, *a_struct_info, fOffset, EXT_STR_h101_TRLO_s467, 0);

    if (!ok)
    {
        LOG(error) << "R3BTrloiiScalerReader_s467::Failed to setup structure information.";
        return kFALSE;
    }

    // Register output array in tree
    FairRootManager::Instance()->Register("TrloiiData", "TRLOII_s467 info", fArray, !fOnline);
    fArray->Clear();

    return kTRUE;
}

Bool_t R3BTrloiiScalerReader_s467::R3BRead()
{
    LOG(debug) << "R3BTrloiiScalerReader_s467::R3BRead() Event data.";

    for (int ch = 0; ch < 16; ++ch)
    {
        new ((*fArray)[fArray->GetEntriesFast()]) R3BTrloiiData(1, ch + 1, fData->TRLORAW_MAIN[ch]);
        new ((*fArray)[fArray->GetEntriesFast()]) R3BTrloiiData(2, ch + 1, fData->TRLOBDT_MAIN[ch]);
        new ((*fArray)[fArray->GetEntriesFast()]) R3BTrloiiData(3, ch + 1, fData->TRLOADT_MAIN[ch]);
        new ((*fArray)[fArray->GetEntriesFast()]) R3BTrloiiData(4, ch + 1, fData->TRLOARD_MAIN[ch]);

        new ((*fArray)[fArray->GetEntriesFast()]) R3BTrloiiData(5, ch + 1, fData->TRLORAW_SCITWO[ch]);
        new ((*fArray)[fArray->GetEntriesFast()]) R3BTrloiiData(6, ch + 1, fData->TRLOBDT_SCITWO[ch]);
        new ((*fArray)[fArray->GetEntriesFast()]) R3BTrloiiData(7, ch + 1, fData->TRLOADT_SCITWO[ch]);
        new ((*fArray)[fArray->GetEntriesFast()]) R3BTrloiiData(8, ch + 1, fData->TRLOARD_SCITWO[ch]);

        new ((*fArray)[fArray->GetEntriesFast()]) R3BTrloiiData(9, ch + 1, fData->TRLORAW_SCIEIGHT[ch]);
        new ((*fArray)[fArray->GetEntriesFast()]) R3BTrloiiData(10, ch + 1, fData->TRLOBDT_SCIEIGHT[ch]);
        new ((*fArray)[fArray->GetEntriesFast()]) R3BTrloiiData(11, ch + 1, fData->TRLOADT_SCIEIGHT[ch]);
        new ((*fArray)[fArray->GetEntriesFast()]) R3BTrloiiData(12, ch + 1, fData->TRLOARD_SCIEIGHT[ch]);
    }

    fNEvent += 1;
    return kTRUE;
}

void R3BTrloiiScalerReader_s467::Reset()
{
    // Reset the output array
    fArray->Clear();
}

ClassImp(R3BTrloiiScalerReader_s467);
