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

#pragma once

#include "R3BReader.h"
#include <TClonesArray.h>

#include <Rtypes.h>

struct EXT_STR_h101_FRSSCI_t;
typedef struct EXT_STR_h101_FRSSCI_t EXT_STR_h101_FRSSCI;

class R3BFrsSciReader : public R3BReader
{
  public:
    // Standard constructor
    R3BFrsSciReader(EXT_STR_h101_FRSSCI*, size_t);
    R3BFrsSciReader(EXT_STR_h101_FRSSCI*, size_t, UShort_t);

    // Destructor
    virtual ~R3BFrsSciReader();

    // Setup structure information
    virtual Bool_t Init(ext_data_struct_info*) override;

    // Read data from full event structure
    virtual Bool_t R3BRead() override;

    // Reset
    virtual void Reset() override;

    // Accessor to select online mode
    void SetOnline(Bool_t option) { fOnline = option; }

  private:
    // Reader specific data structure from ucesb
    EXT_STR_h101_FRSSCI* fData;
    // Data offset
    size_t fOffset;
    // Don't store data for online
    Bool_t fOnline;
    // Output array of type R3BFrsSciMapped
    TClonesArray* fArray;

    UInt_t fNumEntries;
    UShort_t fNumSci;

  public:
    ClassDefOverride(R3BFrsSciReader, 0);
};
