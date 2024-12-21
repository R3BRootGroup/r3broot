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
// This class contains hit-level data for one fiber with MAPMT
// readout on both ends
//

#ifndef R3BFIBERMAPMTHITDATA_H
#define R3BFIBERMAPMTHITDATA_H 1

#include "R3BHit.h"

class R3BFiberMAPMTHitData : public R3BHit
{
  public:
    R3BFiberMAPMTHitData();
    R3BFiberMAPMTHitData(Int_t detId,
                         Double_t x,
                         Double_t y,
                         Double_t eloss,
                         Double_t time,
                         Int_t a_fiber_id,
                         Double_t a_bottom_time_ns,
                         Double_t a_top_time_ns,
                         Double_t a_bottom_tot_ns,
                         Double_t a_top_tot_ns);

    virtual ~R3BFiberMAPMTHitData();

    Int_t GetFiberId() const;
    Double_t GetBottomTime_ns() const;
    Double_t GetTopTime_ns() const;
    Double_t GetTime_ns() const;
    Double_t GetBottomToT_ns() const;
    Double_t GetTopToT_ns() const;

  private:
    Int_t fFiberId;
    Double_t fBottomTime_ns;
    Double_t fTopTime_ns;
    Double_t fTime_ns;
    Double_t fBottomToT_ns;
    Double_t fTopToT_ns;

    ClassDef(R3BFiberMAPMTHitData, 1)
};

#endif
