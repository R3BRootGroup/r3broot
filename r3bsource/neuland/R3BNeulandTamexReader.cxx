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

#include "R3BNeulandTamexReader.h"
#include "R3BLogger.h"
#include "R3BPaddleTamexMappedData.h"

#include <FairLogger.h>
#include <FairRootManager.h>

#include <TClonesArray.h>
#include <ext_data_struct_info.hh>

/**
 ** ext_h101_raw_nnp_tamex.h was created by running
 ** $unpacker --ntuple=STRUCT_HH,RAW:NN,id=h101_raw_nnp_tamex,NOTRIGEVENTNO,ext_h101_raw_nnp_tamex.h
 **/

extern "C"
{
#include "ext_data_client.h"
#include "ext_h101_raw_nnp_tamex.h"
}

R3BNeulandTamexReader::R3BNeulandTamexReader(EXT_STR_h101_raw_nnp_tamex_onion* data, size_t offset)
    : R3BReader("R3BNeulandTamexReader")
    , fData(data)
    , fOffset(offset)
    , fOnline(kFALSE)
    , fSkiptriggertimes(kFALSE)
    , fArray(new TClonesArray("R3BPaddleTamexMappedData"))
    , fArrayTrigger(new TClonesArray("R3BPaddleTamexMappedData"))
    , fNofPlanes(sizeof(fData->NN_P) / sizeof(fData->NN_P[0]))
{
}

R3BNeulandTamexReader::~R3BNeulandTamexReader()
{
    R3BLOG(debug1, "");
    if (fArray)
    {
        delete fArray;
    }
    if (fArrayTrigger)
    {
        delete fArrayTrigger;
    }
}

Bool_t R3BNeulandTamexReader::Init(ext_data_struct_info* a_struct_info)
{
    Int_t ok;
    R3BLOG(info, "");
    EXT_STR_h101_raw_nnp_tamex_ITEMS_INFO(ok, *a_struct_info, fOffset, EXT_STR_h101_raw_nnp_tamex, 0);

    if (!ok)
    {
        R3BLOG(error, "Failed to setup structure information.");
        return kFALSE;
    }

    R3BLOG(info, "Number of planes " << fNofPlanes);

    // Register output arrays in tree
    FairRootManager::Instance()->Register("NeulandMappedData", "Neuland", fArray, !fOnline);
    if (!fSkiptriggertimes)
    {
        FairRootManager::Instance()->Register("NeulandTrigMappedData", "Neuland", fArrayTrigger, !fOnline);
    }
    else
    {
        fArrayTrigger = NULL;
    }
    Reset();
    memset(fData, 0, sizeof *fData);

    return kTRUE;
}

Bool_t R3BNeulandTamexReader::R3BRead()
{
    const auto data = fData;

    for (int plane = 0; plane < fNofPlanes; ++plane)
    {
        for (int pm = 0; pm < 2; ++pm)
        {
            int idx = 0;
            int start = 0;
            int stop = 0;
            int bar;
            int cLE, fLE, cTE, fTE;

            // the counter for coarse time and fine time should be always the same:
            if (data->NN_P[plane].tcl_T[pm].BM != data->NN_P[plane].tfl_T[pm].BM ||
                data->NN_P[plane].tcl_T[pm].B != data->NN_P[plane].tfl_T[pm].B ||
                data->NN_P[plane].tct_T[pm].BM != data->NN_P[plane].tft_T[pm].BM ||
                data->NN_P[plane].tct_T[pm].B != data->NN_P[plane].tft_T[pm].B)
            {
                // fLogger->Info(MESSAGE_ORIGIN, "  Bad event, counter of coarse times and fine times do not match \n");
                return kFALSE;
            }

            // the counter for leading and trailing edge should be always the same:
            if (data->NN_P[plane].tcl_T[pm].B != data->NN_P[plane].tct_T[pm].B ||
                data->NN_P[plane].tfl_T[pm].B != data->NN_P[plane].tft_T[pm].B)
            {
                // fLogger->Info(MESSAGE_ORIGIN, "  Bad event, mismatch of trailing and leading edges \n");
                return kFALSE;
            }

            for (int hit = 0; hit < data->NN_P[plane].tcl_T[pm].BM; hit++)
            {
                bar = data->NN_P[plane].tcl_T[pm].BMI[hit];
                stop = data->NN_P[plane].tcl_T[pm].BME[hit];
                // fLogger->Info(MESSAGE_ORIGIN, "  bar %d, multihit %d \n", bar, stop-start);

                for (int multi = start; multi < stop; multi++)
                {
                    cLE = data->NN_P[plane].tcl_T[pm].Bv[multi];
                    fLE = data->NN_P[plane].tfl_T[pm].Bv[multi];
                    cTE = data->NN_P[plane].tct_T[pm].Bv[multi];
                    fTE = data->NN_P[plane].tft_T[pm].Bv[multi];

                    // fLogger->Info(MESSAGE_ORIGIN, " leading coarse time %d, fine time %d \n", cLE, fLE);
                    // fLogger->Info(MESSAGE_ORIGIN, " trailing coarse time %d, fine time %d \n",cTE, fTE);

                    auto mapped = new ((*fArray)[fArray->GetEntriesFast()]) R3BPaddleTamexMappedData(plane + 1, bar);
                    if (0 == pm)
                    {
                        mapped->fCoarseTime1LE = cLE;
                        mapped->fFineTime1LE = fLE;
                        mapped->fCoarseTime1TE = cTE;
                        mapped->fFineTime1TE = fTE;
                    }
                    else
                    {
                        mapped->fCoarseTime2LE = cLE;
                        mapped->fFineTime2LE = fLE;
                        mapped->fCoarseTime2TE = cTE;
                        mapped->fFineTime2TE = fTE;
                    }
                }
                start = stop;
            }
        }
    }

    // Cards' trigger.
    // NOLINTBEGIN
    if (fArrayTrigger != nullptr)
    {
        int v_i = 0;
        for (int i_i = 0; i_i < fData->NN_TRIGCM; ++i_i)
        {
            auto channel_id = fData->NN_TRIGCMI[i_i];
            for (; v_i < fData->NN_TRIGCME[i_i]; ++v_i)
            {
                auto coarse = fData->NN_TRIGCv[v_i];
                auto fine = fData->NN_TRIGFv[v_i];
                auto mapped =
                    new ((*fArrayTrigger)[fArrayTrigger->GetEntriesFast()]) R3BPaddleTamexMappedData(0, channel_id);
                mapped->fCoarseTime1LE = coarse;
                mapped->fFineTime1LE = fine;
            }
        }
    }
    // NOLINTEND

    return kTRUE;
}

void R3BNeulandTamexReader::Reset()
{
    fArray->Clear();
    if (fArrayTrigger != nullptr)
    {
        fArrayTrigger->Clear();
    }
}

ClassImp(R3BNeulandTamexReader)
