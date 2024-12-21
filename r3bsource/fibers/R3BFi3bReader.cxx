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

#include "R3BFi3bReader.h"
#include "FairLogger.h"

extern "C"
{
#include "ext_data_client.h"
#include "ext_h101_fibthreeb.h"
}

R3BFi3bReader::R3BFi3bReader(EXT_STR_h101_FIBTHREEB* a_data, UInt_t a_offset)
    : R3BBunchedFiberReader("Fi3b", a_offset, 1, 256, 2)
    , fData(reinterpret_cast<EXT_STR_h101_FIBTHREEB_onion*>(a_data))
{
}

Bool_t R3BFi3bReader::Init(ext_data_struct_info* a_struct_info) { R3B_BUNCHED_FIBER_INIT_MAPMT_TRIG(FIBTHREEB, fData); }

ClassImp(R3BFi3bReader)
