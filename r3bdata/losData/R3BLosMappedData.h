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

#ifndef R3BLOSMAPPEDITEM_H
#define R3BLOSMAPPEDITEM_H

#include "TObject.h"

// for the data analysis of the LOS detectors.
// Introduced by Ralf, Jan 2016

class R3BLosMappedData : public TObject
{
  public:
    // Default Constructor
    R3BLosMappedData();

    // Standard Constructor
    R3BLosMappedData(UInt_t, UInt_t, UInt_t, UInt_t, UInt_t);

    UInt_t GetDetector() const;
    UInt_t GetChannel() const;
    UInt_t GetType() const;
    UInt_t GetTimeFine() const;
    UInt_t GetTimeCoarse() const;

  public:
    UInt_t fDetector; // 1..n
    UInt_t fChannel;  // 1..n (n=4 or 8 or 18)
    UInt_t fType;     // 0 (VFTX), 1 (TAMEX leading), 2 (TAMEX trailing)
    UInt_t fTimeCoarse;
    UInt_t fTimeFine;

  public:
    ClassDef(R3BLosMappedData, 1) // 2
};

#endif
