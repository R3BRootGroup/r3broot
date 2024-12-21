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

#ifndef R3BSCI8HITITEM_H
#define R3BSCI8HITITEM_H

#include "TObject.h"

// for the data analysis of the SCI8 detectors.

class R3BSci8HitData : public TObject
{
  public:
    // Default Constructor
    R3BSci8HitData();

    // Standard Constructor
    R3BSci8HitData(UChar_t detector, Double_t t, Double_t x, Double_t Z);

    // Destructor
    virtual ~R3BSci8HitData() {}

    // Getters
    inline const UChar_t& GetDetector() const { return fDetector; }

  private:
    UChar_t fDetector; // 1..n
    Double_t fTime_ns; // average time of the scintillators
    Double_t fX_cm;    // position in cm
    Double_t fZ;       // Charge

  public:
    ClassDef(R3BSci8HitData, 1)
};

#endif
