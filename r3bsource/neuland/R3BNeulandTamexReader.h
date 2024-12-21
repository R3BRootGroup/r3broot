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
#include <Rtypes.h>

class TClonesArray;

struct EXT_STR_h101_raw_nnp_tamex_t;
typedef struct EXT_STR_h101_raw_nnp_tamex_t EXT_STR_h101_raw_nnp_tamex;
typedef struct EXT_STR_h101_raw_nnp_tamex_onion_t EXT_STR_h101_raw_nnp_tamex_onion;
class ext_data_struct_info;

class R3BNeulandTamexReader : public R3BReader
{
  public:
    // Standard constructor
    R3BNeulandTamexReader(EXT_STR_h101_raw_nnp_tamex_onion*, size_t);

    // Destructor
    virtual ~R3BNeulandTamexReader() override;

    // Setup structure information
    auto Init(ext_data_struct_info* /*unused*/) -> Bool_t override;

    // Read data from full event structure
    auto R3BRead() -> Bool_t override;

    // Reset
    void Reset() override;

    // Accessor to select online mode
    inline void SetOnline(Bool_t option) { fOnline = option; }

    // Set the maximum number of planes
    inline void SetMaxNbPlanes(UInt_t max) { fNofPlanes = max; }

    // Accessor to skip trigger times
    inline void SetSkipTriggerTimes() { fSkiptriggertimes = kTRUE; }

  private:
    EXT_STR_h101_raw_nnp_tamex_onion* fData; // Reader specific data structure from ucesb
    size_t fOffset;                          // Data offset
    TClonesArray* fArray = nullptr;          // Output array
    TClonesArray* fArrayTrigger = nullptr;   // Output array, cards' trigger
    UInt_t fNofPlanes;                       // Number of planes
    Bool_t fOnline;                          // Don't store data for online
    Bool_t fSkiptriggertimes;                // Skip trigger times

  public:
    ClassDefOverride(R3BNeulandTamexReader, 0);
};
