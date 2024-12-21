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
#include <FairRunOnline.h>

#include "R3BEventHeader.h"
#include "R3BLogger.h"
#include "R3BTrloiiTpatReader.h"

/**
 ** ext_h101_tpat.h was created by running
 ** $unpacker --ntuple=STRUCT_HH,RAW:TPAT,RAW:TPATFAIL,id=h101_TPAT,NOTRIGEVENTNO,ext_h101_tpat.h
 **/

extern "C"
{
#include "ext_data_client.h"
#include "ext_h101_tpat.h"
}

R3BTrloiiTpatReader::R3BTrloiiTpatReader(EXT_STR_h101_TPAT* data, size_t offset)
    : R3BReader("R3BTrloiiTpatReader")
    , fNEvent(1)
    , fData(data)
    , fOffset(offset)
    , fTrigger(-1)
    , fTpatmin(-1)
    , fTpatmax(-1)
    , fSkipTpatZero(kFALSE)
    , fEventHeader(nullptr)
{
}

R3BTrloiiTpatReader::~R3BTrloiiTpatReader() { R3BLOG(debug, ""); }

Bool_t R3BTrloiiTpatReader::Init(ext_data_struct_info* a_struct_info)
{
    Int_t ok;
    R3BLOG(info, "");
    EXT_STR_h101_TPAT_ITEMS_INFO(ok, *a_struct_info, fOffset, EXT_STR_h101_TPAT, 0);
    if (!ok)
    {
        R3BLOG(fatal, "Failed to setup structure information");
        return kFALSE;
    }

    FairRootManager* frm = FairRootManager::Instance();
    fEventHeader = dynamic_cast<R3BEventHeader*>(frm->GetObject("EventHeader."));
    if (!fEventHeader)
    {
        R3BLOG(fatal, "EventHeader. not found");
        return kFALSE;
    }
    else
    {
        R3BLOG(info, "EventHeader. found");
    }
    Reset();
    memset(fData, 0, sizeof *fData);

    return kTRUE;
}

Bool_t R3BTrloiiTpatReader::R3BRead()
{
    // R3BLOG(info, "Event : " << fNEvent);

    R3BLOG(debug1, "Event data.");
    if (fEventHeader)
    {
        fEventHeader->SetTpat(0);
        if (fData->TPAT > 0)
        {
            fEventHeader->SetTpat(fData->TPATv[0]);
            fNEvent = fEventHeader->GetEventno();
        }
    }
    else if (fData->TPAT <= 0)
        fNEvent++;

    if (0 == (fNEvent % 1000000))
    {
        R3BLOG(debug1, "Event : " << fNEvent);
    }

    if (fTrigger > 0 && fEventHeader && fTpatmin < 0 && fTpatmax < 0)
    {
        if (fEventHeader->GetTrigger() == fTrigger)
        {
            FairRunOnline::Instance()->MarkFill(kTRUE);
        }
        else
        {
            FairRunOnline::Instance()->MarkFill(kFALSE);
        }
    }

    if (fSkipTpatZero && fTrigger < 0 && fTpatmin < 0 && fTpatmax < 0)
    {
        if (fEventHeader->GetTpat() > 0)
        {
            FairRunOnline::Instance()->MarkFill(kTRUE);
        }
        else
        {
            FairRunOnline::Instance()->MarkFill(kFALSE);
        }
    }

    if (fTpatmin >= 0 && fTpatmax >= 0 && fEventHeader && fTrigger < 0)
    {
        Bool_t valid = false;
        Int_t tpatbin = 0;
        for (Int_t i = 0; i < 16; i++)
        {
            tpatbin = false;
            tpatbin = (fEventHeader->GetTpat() & (1 << i));
            if (tpatbin != 0 && (i + 1 < fTpatmin))
            {
                valid = false;
            }
            else
            {
                if (tpatbin != 0 && (i + 1 > fTpatmax))
                {
                    valid = false;
                }
                else if (tpatbin != 0)
                {
                    R3BLOG(debug1, "Accepted Tpat" << i + 1);
                    valid = true;
                }
            }
        }
        if (valid)
        {
            FairRunOnline::Instance()->MarkFill(kTRUE);
        }
        else
        {
            FairRunOnline::Instance()->MarkFill(kFALSE);
        }
    }

    if (fTpatmin >= 0 && fTpatmax >= 0 && fEventHeader && fTrigger > 0)
    {
        if (fEventHeader->GetTrigger() == fTrigger)
        {
            Bool_t valid = false;
            Int_t tpatbin = 0;
            for (Int_t i = 0; i < 16; i++)
            {
                tpatbin = false;
                tpatbin = (fEventHeader->GetTpat() & (1 << i));
                if (tpatbin != 0 && (i + 1 < fTpatmin))
                {
                    valid = false;
                }
                else
                {
                    if (tpatbin != 0 && (i + 1 > fTpatmax))
                    {
                        valid = false;
                    }
                    else if (tpatbin != 0)
                    {
                        R3BLOG(debug1, "Accepted Tpat" << i + 1);
                        valid = true;
                    }
                }
            }
            if (valid)
            {
                FairRunOnline::Instance()->MarkFill(kTRUE);
            }
            else
            {
                FairRunOnline::Instance()->MarkFill(kFALSE);
            }
        }
        else
        {
            FairRunOnline::Instance()->MarkFill(kFALSE);
        }
    }

    return kTRUE;
}

void R3BTrloiiTpatReader::Reset() {}

ClassImp(R3BTrloiiTpatReader)
