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

#ifndef R3BROLUHITITEM_H
#define R3BROLUHITITEM_H

#include "TObject.h"

class R3BRoluHitData : public TObject
{
  public:
    // Default Constructor
    R3BRoluHitData();

    // Standard Constructor
    R3BRoluHitData(Int_t detector, Int_t channel, Double_t time, Double_t tot);

    // Destructor
    virtual ~R3BRoluHitData() {}

    // Getters
    inline const Int_t& GetDetector() const { return fDetector; }
    inline const Int_t& GetChannel() const { return fChannel; }
    inline const Double_t& GetTime() const { return ftime; }
    inline const Double_t& GetToT() const { return ftot; }

  private:
    Int_t fDetector; // 1..n
    Int_t fChannel;
    Double_t ftot;
    Double_t ftime;

    ClassDef(R3BRoluHitData, 1)
};

#endif
