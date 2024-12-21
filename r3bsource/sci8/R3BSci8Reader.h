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

#ifndef R3BSCI8READER_H
#define R3BSCI8READER_H

#include "FairTask.h"
#include "R3BReader.h"
#include <TClonesArray.h>
struct EXT_STR_h101_SCI8_t;
typedef struct EXT_STR_h101_SCI8_t EXT_STR_h101_SCI8;
class FairLogger;
class TH1F;
class TH2F;

class R3BSci8Reader : public R3BReader
{
  public:
    R3BSci8Reader(EXT_STR_h101_SCI8*, UInt_t);
    ~R3BSci8Reader();

    Bool_t Init(ext_data_struct_info*);
    Bool_t R3BRead();
    void Reset();
    virtual void FinishTask();

  private:
    /* Reader specific data structure from ucesb */
    EXT_STR_h101_SCI8* fData;
    /* Data offset */
    UInt_t fOffset;
    /* FairLogger */
    FairLogger* fLogger;
    /* the structs of type R3BSci8xMappedItem */
    TClonesArray* fArray; /**< Output array. */
    Int_t fNEvents = 0;

  public:
    ClassDef(R3BSci8Reader, 0);
};

#endif
