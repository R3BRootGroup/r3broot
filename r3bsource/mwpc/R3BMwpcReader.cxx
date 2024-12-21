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
#include "R3BMwpcMappedData.h"
#include "R3BMwpcReader.h"

extern "C"
{
#include "ext_data_client.h"
#include "ext_h101_mwpc.h"
}

#define NUM_SOFMWPC_DETECTORS 4
#define NUM_SOFMWPC_PLANES_MAX 3

using namespace std;

R3BMwpcReader::R3BMwpcReader(EXT_STR_h101_SOFMWPC* data, size_t offset)
    : R3BReader("R3BMwpcReader")
    , fMaxDet(4)
    , fData(data)
    , fOffset(offset)
    , fOnline(kFALSE)
    , fArrayMwpc0(new TClonesArray("R3BMwpcMappedData")) // class name
    , fArrayMwpc1(new TClonesArray("R3BMwpcMappedData")) // class name
    , fArrayMwpc2(new TClonesArray("R3BMwpcMappedData")) // class name
    , fArrayMwpc3(new TClonesArray("R3BMwpcMappedData")) // class name
{
}

R3BMwpcReader::~R3BMwpcReader()
{
    R3BLOG(debug1, "Destructor");
    if (fArrayMwpc0)
    {
        delete fArrayMwpc0;
    }
    if (fArrayMwpc1)
    {
        delete fArrayMwpc1;
    }
    if (fArrayMwpc2)
    {
        delete fArrayMwpc2;
    }
    if (fArrayMwpc3)
    {
        delete fArrayMwpc3;
    }
}

Bool_t R3BMwpcReader::Init(ext_data_struct_info* a_struct_info)
{
    Int_t ok;
    R3BLOG(info, "");
    EXT_STR_h101_SOFMWPC_ITEMS_INFO(ok, *a_struct_info, fOffset, EXT_STR_h101_SOFMWPC, 0);
    if (!ok)
    {
        R3BLOG(error, "Failed to setup structure information.");
        return kFALSE;
    }

    // Register output array in tree
    FairRootManager::Instance()->Register("Mwpc0MappedData", "MWPC0", fArrayMwpc0, !fOnline);
    if (fMaxDet > 1)
        FairRootManager::Instance()->Register("Mwpc1MappedData", "MWPC1", fArrayMwpc1, !fOnline);
    if (fMaxDet > 2)
        FairRootManager::Instance()->Register("Mwpc2MappedData", "MWPC2", fArrayMwpc2, !fOnline);
    if (fMaxDet > 3)
        FairRootManager::Instance()->Register("Mwpc3MappedData", "MWPC3", fArrayMwpc3, !fOnline);
    R3BLOG(info, "Max. number of MWPCs: " << fMaxDet);
    Reset();
    memset(fData, 0, sizeof *fData);

    // clear struct_writer's output struct. Seems ucesb doesn't do that
    // for channels that are unknown to the current ucesb config.
    // no multi-hit capabiity
    auto* data = reinterpret_cast<EXT_STR_h101_SOFMWPC_onion*>(fData);
    for (int d = 0; d < NUM_SOFMWPC_DETECTORS; d++)
        for (int p = 0; p < NUM_SOFMWPC_PLANES_MAX; p++)
            data->SOFMWPC[d].Plane[p].Q = 0;
    return kTRUE;
}

Bool_t R3BMwpcReader::R3BRead()
{
    // Convert plain raw data to multi-dimensional array
    auto* data = reinterpret_cast<EXT_STR_h101_SOFMWPC_onion*>(fData);

    // loop over all detectors
    for (int d = 0; d < NUM_SOFMWPC_DETECTORS; d++)
    {
        for (int p = 0; p < NUM_SOFMWPC_PLANES_MAX; p++)
        {
            uint32_t numberOfPadsPerPlane = data->SOFMWPC[d].Plane[p].Q;
            for (int mult = 0; mult < numberOfPadsPerPlane; mult++)
            {
                uint16_t pad = data->SOFMWPC[d].Plane[p].QI[mult];
                uint16_t qval = data->SOFMWPC[d].Plane[p].Qv[mult];
                switch (d)
                {
                    case 0:
                        new ((*fArrayMwpc0)[fArrayMwpc0->GetEntriesFast()]) R3BMwpcMappedData(p + 1, pad, qval);
                        break;
                    case 1:
                        new ((*fArrayMwpc1)[fArrayMwpc1->GetEntriesFast()]) R3BMwpcMappedData(p + 1, pad, qval);
                        break;
                    case 2:
                        new ((*fArrayMwpc2)[fArrayMwpc2->GetEntriesFast()]) R3BMwpcMappedData(p + 1, pad, qval);
                        break;
                    case 3:
                        new ((*fArrayMwpc3)[fArrayMwpc3->GetEntriesFast()]) R3BMwpcMappedData(p + 1, pad, qval);
                        break;
                    default:
                        break;
                }
            }
        }
    }

    return kTRUE;
}

void R3BMwpcReader::Reset()
{
    // Reset the output array
    fArrayMwpc0->Clear();
    fArrayMwpc1->Clear();
    fArrayMwpc2->Clear();
    fArrayMwpc3->Clear();
}

ClassImp(R3BMwpcReader)
