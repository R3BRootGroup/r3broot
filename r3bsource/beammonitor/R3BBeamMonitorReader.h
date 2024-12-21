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

#ifndef R3BBEAMMONITORREADER_H
#define R3BBEAMMONITORREADER_H 1

#include "FairTask.h"
#include "R3BReader.h"
#include <TClonesArray.h>

struct EXT_STR_h101_BMON_t;
typedef struct EXT_STR_h101_BMON_t EXT_STR_h101_BMON;

struct EXT_STR_h101_BMON_onion_t;
typedef struct EXT_STR_h101_BMON_onion_t EXT_STR_h101_BMON_onion;
class R3BBeamMonitorReader : public R3BReader
{
  public:
    // Standard constructor
    R3BBeamMonitorReader(EXT_STR_h101_BMON*, size_t);

    // Destructor
    virtual ~R3BBeamMonitorReader();

    // Setup structure information
    virtual Bool_t Init(ext_data_struct_info*) override;

    // Read data from full event structure
    virtual Bool_t R3BRead() override;

    // Reset
    virtual void Reset() override;

    // Accessor to select online mode
    void SetOnline(Bool_t option) { fOnline = option; }

  private:
    /* Reader specific data structure from ucesb */
    EXT_STR_h101_BMON_onion* fData;
    /* Data offset */
    size_t fOffset;
    // Don't store data for online
    Bool_t fOnline;
    TClonesArray* fArray; /**< Output array. */
    Int_t fNEvents = 0;

  public:
    ClassDefOverride(R3BBeamMonitorReader, 0);
};

#endif
