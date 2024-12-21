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

#ifndef R3BBUNCHEDFIBERHITDATA_H
#define R3BBUNCHEDFIBERHITDATA_H

#include "R3BHit.h"

class R3BBunchedFiberHitData : public R3BHit
{
  public:
    R3BBunchedFiberHitData();
    R3BBunchedFiberHitData(Int_t detId,
                           Double_t x,
                           Double_t y,
                           Double_t eloss,
                           Double_t time,
                           Int_t a_fiber_id,
                           Double_t a_mapmt_time_ns,
                           Double_t a_spmt_time_ns,
                           Double_t a_mapmt_tot_ns,
                           Double_t a_spmt_tot_ns);

    virtual ~R3BBunchedFiberHitData();

    Int_t GetFiberId() const;
    Double_t GetMAPMTTime_ns() const;
    Double_t GetSPMTTime_ns() const;
    Double_t GetTime_ns() const;
    Double_t GetMAPMTToT_ns() const;
    Double_t GetSPMTToT_ns() const;

  private:
    Int_t fFiberId;
    Double_t fMAPMTTime_ns;
    Double_t fSPMTTime_ns;
    Double_t fTime_ns;
    Double_t fMAPMTToT_ns;
    Double_t fSPMTToT_ns;

    ClassDef(R3BBunchedFiberHitData, 2)
};

#endif
