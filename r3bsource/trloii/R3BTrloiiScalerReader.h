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

// ----------------------------------------------------------------------
// -----                       R3BTrloiiScalerReader                -----
// -----             Created 12/05/18 by J.L. Rodriguez-Sanchez     -----
// ----------------------------------------------------------------------

#pragma once

#include "R3BReader.h"
#include <Rtypes.h>

class TClonesArray;

struct EXT_STR_h101_TRLO_t;
typedef struct EXT_STR_h101_TRLO_t EXT_STR_h101_TRLO;
typedef struct EXT_STR_h101_TRLO_onion_t EXT_STR_h101_TRLO_onion;
class ext_data_struct_info;

class R3BTrloiiScalerReader : public R3BReader
{
  public:
    // Standard constructor
    R3BTrloiiScalerReader(EXT_STR_h101_TRLO_onion*, size_t);

    // Destructor
    virtual ~R3BTrloiiScalerReader();

    // Setup structure information
    virtual Bool_t Init(ext_data_struct_info*) override;

    // Read data from full event structure
    virtual Bool_t R3BRead() override;

    // Reset
    virtual void Reset() override;

    // Accessor to select online mode
    inline void SetOnline(Bool_t option) { fOnline = option; }

  private:
    // An event counter
    unsigned int fNEvent;
    // Reader specific data structure from ucesb
    EXT_STR_h101_TRLO_onion* fData;
    // Data offset
    size_t fOffset;
    // Don't store data for online
    Bool_t fOnline;
    // Output array
    TClonesArray* fArray;

  public:
    ClassDefOverride(R3BTrloiiScalerReader, 0);
};
