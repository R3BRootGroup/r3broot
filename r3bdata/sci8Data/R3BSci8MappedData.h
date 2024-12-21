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

#ifndef R3BSCI8MAPPEDITEM_H
#define R3BSCI8MAPPEDITEM_H

#include "TObject.h"

// for the data analysis of the SCI8 detectors.

class R3BSci8MappedData : public TObject
{
  public:
    // Default Constructor
    R3BSci8MappedData();

    // Standard Constructor
    R3BSci8MappedData(UInt_t detector, UInt_t channel, UInt_t type, UInt_t timeFine, UInt_t timeCoarse);

    // Destructor
    virtual ~R3BSci8MappedData() {}

    // Getters
    inline const UInt_t& GetDetector() const { return fDetector; }
    inline const UInt_t& GetChannel() const { return fChannel; }
    inline const UInt_t& GetType() const { return fType; }
    inline const UInt_t& GetTimeFine() const { return fTimeFine; }
    inline const UInt_t& GetTimeCoarse() const { return fTimeCoarse; }

    // private:
  public:
    UInt_t fDetector; // 1..n
    UInt_t fChannel;  // 1..n (n=4 or 8 or 18)
    UInt_t fType;     // 0 (VFTX), 1 (TAMEX leading), 2 (TAMEX trailing)
    UInt_t fTimeCoarse;
    UInt_t fTimeFine;
    Int_t fNEvents = 0;

  public:
    ClassDef(R3BSci8MappedData, 1)
};

#endif
