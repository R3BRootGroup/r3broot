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

#ifndef R3BTOFICALITEM_H
#define R3BTOFICALITEM_H

#include "TObject.h"

class R3BTofiCalData : public TObject
{
  public:
    R3BTofiCalData();
    R3BTofiCalData(UInt_t, UInt_t, UInt_t, Double_t, Double_t);

    UInt_t GetDetectorId() const;
    UInt_t GetBarId() const;
    UInt_t GetSideId() const;
    Double_t GetTimeLeading_ns() const;
    Double_t GetTimeTrailing_ns() const;

  private:
    UInt_t fDetector; // 1..4
    UInt_t fBar;      // 1..44
    UInt_t fSide;     // 1..2
    Double_t fLeading_ns;
    Double_t fTrailing_ns;

    ClassDef(R3BTofiCalData, 2)
};

#endif
