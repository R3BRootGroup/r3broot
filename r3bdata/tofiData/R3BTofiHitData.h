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

#ifndef R3BTOFIHITDATA_H
#define R3BTOFIHITDATA_H

#include "R3BHit.h"

// for the data analysis of the Tofi detectors.
// Introduced by M.Heil, May 2016

class R3BTofiHitData : public R3BHit
{
  public:
    // Default Constructor
    R3BTofiHitData();

    // Standard Constructor
    R3BTofiHitData(Double_t t,
                   Double_t x,
                   Double_t y,
                   Double_t Z,
                   Double_t tdiff,
                   Double_t ELoss = 0,
                   Double_t ID = 0,
                   Int_t iBar = 0);

    // Destructor
    virtual ~R3BTofiHitData() {}

    Int_t GetBarId() const;

  private:
    Int_t fBarId;

    ClassDef(R3BTofiHitData, 2)
};

#endif
