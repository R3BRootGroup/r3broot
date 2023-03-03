/******************************************************************************
 *   Copyright (C) 2019 GSI Helmholtzzentrum für Schwerionenforschung GmbH    *
 *   Copyright (C) 2019 Members of R3B Collaboration                          *
 *                                                                            *
 *             This software is distributed under the terms of the            *
 *                 GNU General Public Licence (GPL) version 3,                *
 *                    copied verbatim in the file "LICENSE".                  *
 *                                                                            *
 * In applying this license GSI does not waive the privileges and immunities  *
 * granted to it by virtue of its status as an Intergovernmental Organization *
 * or submit itself to any jurisdiction.                                      *
 ******************************************************************************/

#include "R3BTrloiiTpatReader.h"
#include "FairLogger.h"
#include "FairRootManager.h"
#include "R3BEventHeader.h"
#include "R3BLogger.h"

extern "C"
{
#include "ext_data_client.h"
#include "ext_h101_tpat.h"
}

using namespace std;

R3BTrloiiTpatReader::R3BTrloiiTpatReader(EXT_STR_h101_TPAT* data, UInt_t offset)
    : R3BReader("R3BTrloiiTpatReader")
    , fNEvent(0)
    , fData(data)
    , fOffset(offset)
    , fEventHeader(nullptr)
{
}

R3BTrloiiTpatReader::~R3BTrloiiTpatReader() {}

Bool_t R3BTrloiiTpatReader::Init(ext_data_struct_info* a_struct_info)
{
    int ok;

    EXT_STR_h101_TPAT_ITEMS_INFO(ok, *a_struct_info, fOffset, EXT_STR_h101_TPAT, 0);

    if (!ok)
    {
        perror("ext_data_struct_info_item");
        LOG(error) << "Failed to setup structure information.";
        return kFALSE;
    }

    auto mgr = FairRootManager::Instance();
    fEventHeader = dynamic_cast<R3BEventHeader*>(mgr->GetObject("EventHeader."));
    if (fEventHeader)
        R3BLOG(info, "EventHeader. was found");
    else
        R3BLOG(info, "EventHeader. was not found");

    return kTRUE;
}

Bool_t R3BTrloiiTpatReader::Read()
{

    // LOG(info) << "TrloiiTpatReader::Read BEGIN";

    if (fEventHeader)
    {
        fEventHeader->SetTpat(0);
        if (fData->TPAT > 0)
        {
            fEventHeader->SetTpat(fData->TPATv[0]);
            fNEvent = fEventHeader->GetEventno();
        }
    }
    if (!fEventHeader || fData->TPAT <= 0)
        fNEvent++;

    if (0 == (fNEvent % 1000000))
    {

        LOG(debug1) << "R3BTrloiiTpatReader : event : " << fNEvent;
    }

    /* Display data */
    // char str[256];
    // sprintf(str, "  Trlo II Tpat = 0x%04x.", fData->TPATv[0]);
    // LOG(info) << str;

    // LOG(info) << "TrloiiTpatReader::Read END";

    return kTRUE;
}

void R3BTrloiiTpatReader::Reset() { fNEvent = 0; }

ClassImp(R3BTrloiiTpatReader)
