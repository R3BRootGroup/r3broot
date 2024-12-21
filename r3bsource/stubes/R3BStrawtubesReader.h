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

#ifndef R3BSTRAWTUBESREADER_H
#define R3BSTRAWTUBESREADER_H

#include "R3BReader.h"
#include <TClonesArray.h>

typedef struct EXT_STR_h101_STRAWTUBES_t EXT_STR_h101_STRAWTUBES;
typedef struct EXT_STR_h101_STRAWTUBES_onion_t EXT_STR_h101_STRAWTUBES_onion;
class FairLogger;

class R3BStrawtubesReader : public R3BReader
{
  public:
    explicit R3BStrawtubesReader(EXT_STR_h101_STRAWTUBES*, UInt_t);
    virtual ~R3BStrawtubesReader();

    Bool_t Init(ext_data_struct_info*);
    Bool_t R3BRead();
    void Reset();

  private:
    EXT_STR_h101_STRAWTUBES_onion* fOnion;
    UInt_t fOffset;
    FairLogger* fLogger;
    TClonesArray* fArray;

  public:
    ClassDef(R3BStrawtubesReader, 0);
};

#endif
