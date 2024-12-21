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

//
// SPMT trigger reader, since many fiber detectors share TAMEX cards and we
// don't want to recalibrate channels in many places.
//

#ifndef R3BBUNCHEDFIBERSPMTTRIGREADER_H
#define R3BBUNCHEDFIBERSPMTTRIGREADER_H

#include "R3BReader.h"
#include <Rtypes.h>

class TClonesArray;

struct EXT_STR_h101_FIB_t;
typedef struct EXT_STR_h101_FIB_t EXT_STR_h101_FIB;
typedef struct EXT_STR_h101_FIB_onion_t EXT_STR_h101_FIB_onion;
class ext_data_struct_info;

class R3BBunchedFiberSPMTTrigReader : public R3BReader
{
  public:
    // Standard constructor
    R3BBunchedFiberSPMTTrigReader(EXT_STR_h101_FIB*, size_t);

    // Destructor
    virtual ~R3BBunchedFiberSPMTTrigReader();

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
    EXT_STR_h101_FIB_onion* fData;
    // Data offset
    size_t fOffset;
    // Don't store data for online
    Bool_t fOnline;
    // Output array
    TClonesArray* fMappedArray;

  public:
    ClassDefOverride(R3BBunchedFiberSPMTTrigReader, 0);
};

#endif
