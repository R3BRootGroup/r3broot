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

#ifndef R3BSCI2READER_H
#define R3BSCI2READER_H

#include "FairTask.h"
#include "R3BReader.h"
#include <Rtypes.h>
#include <TClonesArray.h>

struct EXT_STR_h101_SCI2_t;
typedef struct EXT_STR_h101_SCI2_t EXT_STR_h101_SCI2;
class ext_data_struct_info;

class R3BSci2Reader : public R3BReader
{
  public:
    // Standard constructor
    R3BSci2Reader(EXT_STR_h101_SCI2*, size_t);

    // Destructor
    virtual ~R3BSci2Reader();

    // Setup structure information
    virtual Bool_t Init(ext_data_struct_info*) override;

    // Read data from full event structure
    virtual Bool_t R3BRead() override;

    // Reset
    virtual void Reset() override;

    // Accessor to select online mode
    void SetOnline(Bool_t option) { fOnline = option; }

  private:
    // An event counter
    unsigned int fNEvent;
    // Reader specific data structure from ucesb
    EXT_STR_h101_SCI2* fData;
    // Data offset
    size_t fOffset;
    // Don't store data for online
    Bool_t fOnline;
    // Output array of type R3BSci2xMappedItem
    TClonesArray* fArray;

  public:
    ClassDefOverride(R3BSci2Reader, 0);
};

#endif /* R3BSCI2READER_H */
