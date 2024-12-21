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
// -----                       R3BAlpideReader                     -----
// -----             Created 09/02/22 by J.L. Rodriguez-Sanchez    -----
// ---------------------------------------------------------------------

#ifndef R3BALPIDEREADER_H
#define R3BALPIDEREADER_H 1

#include "R3BReader.h"
#include <stdint.h>

class TClonesArray;

struct EXT_STR_h101_ALPIDE_t;
typedef struct EXT_STR_h101_ALPIDE_t EXT_STR_h101_ALPIDE;
typedef struct EXT_STR_h101_ALPIDE_onion_t EXT_STR_h101_ALPIDE_onion;
class ext_data_struct_info;

class R3BAlpideReader : public R3BReader
{
  public:
    // Standard constructor
    R3BAlpideReader(EXT_STR_h101_ALPIDE_onion*, size_t);

    // Destructor
    virtual ~R3BAlpideReader();

    // Setup structure information
    virtual Bool_t Init(ext_data_struct_info*) override;

    // Read data from full event structure
    virtual Bool_t R3BRead() override;

    // Reset
    virtual void Reset() override;

    // Accessor to select online mode
    void SetOnline(Bool_t option) { fOnline = option; }

    // Set number of sensors
    void SetNbDet(uint n) { fNbDet = n; }

  private:
    // An event counter
    unsigned int fNEvent;
    // Number of ALPIDE sensors
    uint fNbDet;
    // Reader specific data structure from ucesb
    EXT_STR_h101_ALPIDE_onion* fData;
    // Data offset
    size_t fOffset;
    // Don't store data for online
    Bool_t fOnline;
    // Output array
    TClonesArray* fArray;

  public:
    ClassDefOverride(R3BAlpideReader, 1);
};

#endif /* R3BALPIDEREADER_H */
