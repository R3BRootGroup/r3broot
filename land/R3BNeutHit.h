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

// ------------------------------------------------------------------------
// -----                        R3BNeutHit                            -----
// -----              Created 21.05.2012 by D.Kresan                  -----
// ------------------------------------------------------------------------

#ifndef R3BNEUTHIT_H
#define R3BNEUTHIT_H

#include "TMath.h"
#include "TObject.h"

class R3BNeutHit : public TObject
{

  public:
    R3BNeutHit();
    R3BNeutHit(Double_t x, Double_t y, Double_t z, Double_t t);
    virtual ~R3BNeutHit();

  private:
    Double32_t fX;
    Double32_t fY;
    Double32_t fZ;
    Double32_t fT;

  public:
    inline Double_t GetX() const { return fX; }
    inline Double_t GetY() const { return fY; }
    inline Double_t GetZ() const { return fZ; }
    inline Double_t GetT() const { return fT; }
    inline Double_t GetD() const { return TMath::Sqrt(fX * fX + fY * fY + fZ * fZ); }
    Double_t GetP() const;

    ClassDef(R3BNeutHit, 1)
};

#endif
