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
// This class contains hit-level data for one fiber.
//

#ifndef R3BSFIBHITDATA_H
#define R3BSFIBHITDATA_H

#include "R3BHit.h"

class R3BSfibHitData : public R3BHit
{
  public:
    R3BSfibHitData();
    R3BSfibHitData(Int_t, Double_t, Double_t, Double_t, Double_t);

    virtual ~R3BSfibHitData();

    Int_t GetFiberId() const;
    Double_t GetTopTime_ns() const;
    Double_t GetBotTime_ns() const;
    Double_t GetTime_ns() const;
    Double_t GetTopToT_ns() const;
    Double_t GetBotToT_ns() const;

  private:
    Int_t fFiberId;
    Double_t fTopTime_ns;
    Double_t fBotTime_ns;
    Double_t fTime_ns;
    Double_t fTopToT_ns;
    Double_t fBotToT_ns;

    ClassDef(R3BSfibHitData, 1)
};

#endif
