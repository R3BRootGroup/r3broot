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

// -----------------------------------------------------------------------------
// -----                              R3BHit                               -----
// -----                 Created on 01.03.2018 by D.Kresan                 -----
// -----------------------------------------------------------------------------

#ifndef R3BHIT_H
#define R3BHIT_H

#include "TObject.h"

class R3BHit : public TObject
{
  public:
    R3BHit();
    R3BHit(Int_t detId, Double_t x, Double_t y, Double_t eloss, Double_t time, Int_t hitId = -1);
    virtual ~R3BHit();

    inline const Int_t& GetDetId() const { return fDetId; }
    inline const Double_t& GetX() const { return fX; }
    inline const Double_t& GetY() const { return fY; }
    inline const Double_t& GetEloss() const { return fEloss; }
    inline const Double_t& GetTime() const { return fTime; }
    inline const Int_t& GetHitId() const { return fHitId; }

    inline void SetEloss(const Double_t& eloss) { fEloss = eloss; }
    inline void SetTime(const Double_t& time) { fTime = time; }
    inline void SetHitId(const Int_t& hitId) { fHitId = hitId; }

  protected:
    Int_t fDetId;
    Double_t fX;
    Double_t fY;
    Double_t fEloss;
    Double_t fTime;
    Int_t fHitId;

  public:
    ClassDef(R3BHit, 1)
};

#endif
