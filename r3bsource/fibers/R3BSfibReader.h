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

#ifndef R3BSFIBREADER_H
#define R3BSFIBREADER_H

#include "R3BReader.h"
#include <TClonesArray.h>

struct EXT_STR_h101_sfib_t;
typedef struct EXT_STR_h101_sfib_t EXT_STR_h101_sfib;

class FairLogger;
class R3BEventHeader;

class R3BSfibReader : public R3BReader
{
  public:
    R3BSfibReader(EXT_STR_h101_sfib*, UInt_t);

    Bool_t Init(ext_data_struct_info*);
    Bool_t R3BRead();
    void Reset();

  private:
    /* An event counter */
    UInt_t fNEvent;
    /* Reader specific data structure from ucesb */
    EXT_STR_h101_sfib* fData;
    /* Offset of detector specific data in full data structure */
    UInt_t fOffset;
    /* FairLogger */
    FairLogger* fLogger;
    TClonesArray* fMappedArray;

  public:
    ClassDef(R3BSfibReader, 0);
};

#endif
