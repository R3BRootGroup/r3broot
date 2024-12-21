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

#ifndef R3BPDCREADER_H
#define R3BPDCREADER_H

#include "R3BReader.h"
#include <TClonesArray.h>
struct EXT_STR_h101_PDC_t;
typedef struct EXT_STR_h101_PDC_t EXT_STR_h101_PDC;

class R3BPdcReader : public R3BReader
{
  public:
    R3BPdcReader(EXT_STR_h101_PDC*, UInt_t);
    ~R3BPdcReader();

    Bool_t Init(ext_data_struct_info*);
    Bool_t R3BRead();
    void Reset();

  private:
    /* Reader specific data structure from ucesb */
    EXT_STR_h101_PDC* fData;
    /* Data offset */
    UInt_t fOffset;

    TClonesArray* fMappedArray;

  public:
    ClassDef(R3BPdcReader, 1);
};

#endif
