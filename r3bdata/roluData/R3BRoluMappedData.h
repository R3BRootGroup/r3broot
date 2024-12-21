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

#ifndef R3BROLUMAPPEDITEM_H
#define R3BROLUMAPPEDITEM_H

#include "TObject.h"

// for the data analysis of ROLU following LOS

class R3BRoluMappedData : public TObject
{
  public:
    // Default Constructor
    R3BRoluMappedData();

    // Standard Constructor
    R3BRoluMappedData(UInt_t, UInt_t, Int_t, UInt_t, UInt_t);

    UInt_t GetDetector() const;
    UInt_t GetChannel() const;
    Int_t GetType() const;
    UInt_t GetTimeFine() const;
    UInt_t GetTimeCoarse() const;

  public:
    UInt_t fDetector; // 1..n
    UInt_t fChannel;  // 1..n (n=4 )
    Int_t fType;      // 0 (TAMEX leading), 1 (TAMEX trailing)
    UInt_t fTimeCoarse;
    UInt_t fTimeFine;

  public:
    ClassDef(R3BRoluMappedData, 1)
};

#endif // R3BROLUMAPPEDITEM_H
