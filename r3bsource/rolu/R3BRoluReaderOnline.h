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

#ifndef R3BROLUREADERONLINE_H
#define R3BROLUREADERONLINE_H 1

#include "R3BReader.h"
#include <Rtypes.h>

class TClonesArray;

struct EXT_STR_h101_ROLU_NOTRIG_t;
typedef struct EXT_STR_h101_ROLU_NOTRIG_t EXT_STR_h101_ROLU_NOTRIG;
typedef struct EXT_STR_h101_ROLU_NOTRIG_onion_t EXT_STR_h101_ROLU_NOTRIG_onion;
class ext_data_struct_info;

class R3BRoluReaderOnline : public R3BReader
{
  public:
    // Standard constructor
    R3BRoluReaderOnline(EXT_STR_h101_ROLU_NOTRIG_onion*, size_t);

    // Destructor
    virtual ~R3BRoluReaderOnline();

    // Setup structure information
    virtual bool Init(ext_data_struct_info*) override;

    // Read data from full event structure
    virtual bool R3BRead() override;

    // Reset
    virtual void Reset() override;

    // Accessor to select online mode
    void SetOnline(bool option) { fOnline = option; }

    // Accessor to select the number of detectors
    void SetNbDet(unsigned int nbb) { fNbDet = nbb; }

  private:
    // Reader specific data structure from ucesb
    EXT_STR_h101_ROLU_NOTRIG_onion* fData;
    // Number of detectors
    unsigned int fNbDet;
    // Data offset
    size_t fOffset;
    // Output array
    TClonesArray* fArray;
    // Don't store data for online
    bool fOnline;

    int fNEvents = 0;

  public:
    ClassDefOverride(R3BRoluReaderOnline, 0);
};

#endif // R3BROLUREADERONLINE_H
