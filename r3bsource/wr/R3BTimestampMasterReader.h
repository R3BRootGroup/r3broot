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

#ifndef R3BTIMESTAMPMASTERREADER_H
#define R3BTIMESTAMPMASTERREADER_H

#include "R3BReader.h"
#include "R3BWhiterabbitReaderImpl.h"

struct EXT_STR_h101_timestamp_master_t;
typedef struct EXT_STR_h101_timestamp_master_t EXT_STR_h101_timestamp_master;

class R3BEventHeader;

class R3BTimestampMasterReader : public R3BReader
{
  public:
    R3BTimestampMasterReader(EXT_STR_h101_timestamp_master*, size_t);

    Bool_t Init(ext_data_struct_info*);
    Bool_t R3BRead();
    void Reset();

  private:
    UInt_t fNEvent;
    EXT_STR_h101_timestamp_master* fData;
    size_t fOffset;
    UInt_t fWhiterabbitId;
    R3BEventHeader* fEventHeader;

  public:
    ClassDef(R3BTimestampMasterReader, 0);
};

#endif
