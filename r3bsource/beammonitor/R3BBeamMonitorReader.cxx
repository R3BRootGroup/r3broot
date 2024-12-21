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

#include "R3BBeamMonitorReader.h"
#include "FairRootManager.h"
#include "R3BBeamMonitorMappedData.h"
#include "R3BLogger.h"
#include "TClonesArray.h"

extern "C"
{
#include "ext_data_client.h"
#include "ext_h101_bmon.h"
}
#include "TMath.h"
#define IS_NAN(x) TMath::IsNaN(x)

using namespace std;

R3BBeamMonitorReader::R3BBeamMonitorReader(EXT_STR_h101_BMON* data, size_t offset)
    : R3BReader("R3BBeamMonitorReader")
    , fData(reinterpret_cast<EXT_STR_h101_BMON_onion*>(data))
    , fOffset(offset)
    , fOnline(kFALSE)
    , fArray(new TClonesArray("R3BBeamMonitorMappedData"))
{
}

R3BBeamMonitorReader::~R3BBeamMonitorReader()
{
    if (fArray)
        delete fArray;
}

Bool_t R3BBeamMonitorReader::Init(ext_data_struct_info* a_struct_info)
{
    Int_t ok;
    R3BLOG(info, "");
    EXT_STR_h101_BMON_ITEMS_INFO(ok, *a_struct_info, fOffset, EXT_STR_h101_BMON, 0);
    if (!ok)
    {
        R3BLOG(error, "Failed to setup structure information.");
        return kFALSE;
    }

    // Register output array in tree
    FairRootManager::Instance()->Register("BeamMonitorMapped", "BeamMonitorMapped data", fArray, !fOnline);
    Reset();
    memset(fData, 0, sizeof *fData);

    // clear struct_writer's output struct. Seems ucesb doesn't do that
    // for channels that are unknown to the current ucesb config.
    auto* data = fData;
    data->IC = 0;
    data->SEETRAM = 0;
    data->TOFDOR = 0;

    return kTRUE;
}

Bool_t R3BBeamMonitorReader::R3BRead()
{
    // Convert plain raw data to multi-dimensional array
    EXT_STR_h101_BMON_onion* data = fData;

    new ((*fArray)[fArray->GetEntriesFast()]) R3BBeamMonitorMappedData(data->IC, data->SEETRAM, data->TOFDOR);

    fNEvents += 1;

    return kTRUE;
}

void R3BBeamMonitorReader::Reset()
{
    // Reset the output array
    fArray->Clear();
}

ClassImp(R3BBeamMonitorReader)
