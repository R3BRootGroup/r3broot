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
// -----                       R3BTrloiiSampReader                  -----
// ----------------------------------------------------------------------
// sampler reader modified to read the sampler storing master start timestamps
// along with sampler hits and both the high and low bits

#pragma once

#include "R3BReader.h"
#include <TClonesArray.h>

struct EXT_STR_h101_SAMP_t;
typedef struct EXT_STR_h101_SAMP_t EXT_STR_h101_SAMP;
typedef struct EXT_STR_h101_SAMP_onion_t EXT_STR_h101_SAMP_onion;

struct EXT_STR_h101_SAMPLOS_t;
typedef struct EXT_STR_h101_SAMPLOS_t EXT_STR_h101_SAMPLOS;
typedef struct EXT_STR_h101_SAMPLOS_onion_t EXT_STR_h101_SAMPLOS_onion;

struct EXT_STR_h101_SAMPMS_t;
typedef struct EXT_STR_h101_SAMPMS_t EXT_STR_h101_SAMPMS;
typedef struct EXT_STR_h101_SAMPMS_onion_t EXT_STR_h101_SAMPMS_onion;

struct EXT_STR_h101_SAMPLOSMS_t;
typedef struct EXT_STR_h101_SAMPLOSMS_t EXT_STR_h101_SAMPLOSMS;
typedef struct EXT_STR_h101_SAMPLOSMS_onion_t EXT_STR_h101_SAMPLOSMS_onion;

class R3BTrloiiSampReader : public R3BReader
{
  public:
    // Parametrized constructors
    R3BTrloiiSampReader(EXT_STR_h101_SAMP_onion*, size_t);
    R3BTrloiiSampReader(EXT_STR_h101_SAMPLOS_onion*, size_t);
    R3BTrloiiSampReader(EXT_STR_h101_SAMPMS_onion*, size_t);
    R3BTrloiiSampReader(EXT_STR_h101_SAMPLOSMS_onion*, size_t);

    // Destructor
    virtual ~R3BTrloiiSampReader();

    // Setup structure information
    virtual Bool_t Init(ext_data_struct_info*) override;

    // Read data from full event structure
    virtual Bool_t R3BRead() override;

    // Reset
    virtual void Reset() override;

  private:
    /* Reader specific data structure from ucesb */
    EXT_STR_h101_SAMP_onion* fSampData;
    EXT_STR_h101_SAMPLOS_onion* fSampLosData;
    EXT_STR_h101_SAMPMS_onion* fSampMSData;
    EXT_STR_h101_SAMPLOSMS_onion* fSampLosMSData;

    /* Offset of detector specific data in full data structure */
    size_t fOffset;

    /* number of samplers */
    UInt_t fSNum;

    /* number of sampler components */
    UInt_t fSCNum;

    /* Output arrays */
    TClonesArray* fArray;
    TClonesArray* fArrayLH;

    /* Should data for online be stored */
    Bool_t fOnline;

  public:
    ClassDefOverride(R3BTrloiiSampReader, 0);
};
