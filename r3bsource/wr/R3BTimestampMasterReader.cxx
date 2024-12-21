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

#include "R3BTimestampMasterReader.h"
#include "FairLogger.h"
#include "FairRootManager.h"
#include "R3BEventHeader.h"
#include "R3BWhiterabbitReaderImpl.h"

extern "C"
{
#include "ext_data_client.h"
#include "ext_h101_timestamp_master.h"
}

R3B_WHITERABBIT_READER_IMPL(TimestampMaster, timestamp_master, 0x1000);
// 0x300 for dec2019
// 0x100 for data before dec2019
// 0x1000 for S455 forwards.
