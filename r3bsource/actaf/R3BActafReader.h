/******************************************************************************
 *   Copyright (C) 2025 GSI Helmholtzzentrum für Schwerionenforschung GmbH    *
 *   Copyright (C) 2025 Members of R3B Collaboration                          *
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
#include <Rtypes.h>
#include <memory>

class TClonesArray;

struct EXT_STR_h101_ACTAF_t;
typedef struct EXT_STR_h101_ACTAF_t EXT_STR_h101_ACTAF;
typedef struct EXT_STR_h101_ACTAF_onion_t EXT_STR_h101_ACTAF_onion;
class ext_data_struct_info;
/**
 * A reader of ACTAF data with UCESB.
 * Receives mapped raw data and converts it to R3BRoot objects.
 * @author P. Rusell
 * @since Feb 17, 2025
 */
class R3BActafReader : public R3BReader
{
  public:
    // Standard constructor
    R3BActafReader(EXT_STR_h101_ACTAF_onion*, size_t);

    // Destructor
    virtual ~R3BActafReader() = default;

    // Setup structure information
    auto Init(ext_data_struct_info* /*unused*/) -> Bool_t override;

    // Read data from full event structure
    auto R3BRead() -> Bool_t override;

    // Reset
    void Reset() override;

    // Accessor to select online mode
    inline void SetOnline(Bool_t option) { fOnline = option; }

  private:
    // An event counter
    unsigned int fNEvent = 0;
    // Reader specific data structure from ucesb
    EXT_STR_h101_ACTAF_onion* fData;
    // Data offset
    size_t fOffset;
    // Don't store data for online
    Bool_t fOnline = kFALSE;
    // Output array
    std::unique_ptr<TClonesArray> fArray;

  public:
    ClassDefOverride(R3BActafReader, 0);
};
