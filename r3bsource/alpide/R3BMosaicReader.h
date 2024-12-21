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

// ---------------------------------------------------------------------
// -----                       R3BMosaicReader                     -----
// -----             Created 04/02/24 by J.L. Rodriguez-Sanchez    -----
// ---------------------------------------------------------------------

#pragma once

#include "R3BReader.h"
#include <stdint.h>
#include <vector>

class TClonesArray;

struct EXT_STR_h101_MOSAIC_t;
typedef struct EXT_STR_h101_MOSAIC_onion_t EXT_STR_h101_MOSAIC_onion;
class ext_data_struct_info;

class R3BMosaicReader : public R3BReader
{
  public:
    // Standard constructor
    R3BMosaicReader(EXT_STR_h101_MOSAIC_onion*, size_t);

    // Destructor
    virtual ~R3BMosaicReader();

    // Setup structure information
    Bool_t Init(ext_data_struct_info*) override;

    // Read data from full event structure
    Bool_t R3BRead() override;

    // Reset
    void Reset() override;

    // Accessor to select online mode
    inline void SetOnline(bool option = true) { fOnline = option; }

  private:
    // An event counter
    unsigned int fNEvent = 1;
    // Reader specific data structure from ucesb
    EXT_STR_h101_MOSAIC_onion* fData;
    // Number of Mosaics
    int fNbMosaic = 1;
    // Data offset
    size_t fOffset = 0;
    // Don't store data for online
    bool fOnline = false;
    // Output array
    TClonesArray* fArray = nullptr;

    std::vector<int> map_mosaics = { 0, 0, 0, 0, 1, 2, 3, 0 };

    std::vector<int> map_sensors = { 1,  2,  3,  4,  5,  6,  7,  8,  9,  10, 11, 12,
                                     18, 17, 16, 15, 14, 13, 24, 23, 22, 21, 20, 19 };

  public:
    ClassDefOverride(R3BMosaicReader, 1);
};
