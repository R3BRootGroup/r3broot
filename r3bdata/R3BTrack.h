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
// -----                              R3BTrack                             -----
// -----                 Created on 09.03.2020 by M.Heil                   -----
// -----------------------------------------------------------------------------

#ifndef R3BTrack_H
#define R3BTrack_H

#include "TObject.h"

class R3BTrack : public TObject
{
  public:
    R3BTrack();
    R3BTrack(Double_t x,
             Double_t y,
             Double_t z,
             Double_t px,
             Double_t py,
             Double_t pz,
             Double_t q,
             Double_t AoZ,
             Double_t chix,
             Double_t chiy,
             Int_t quality);
    virtual ~R3BTrack();

    inline const Double_t& GetX() const { return fX; }
    inline const Double_t& GetY() const { return fY; }
    inline const Double_t& GetZ() const { return fZ; }
    inline const Double_t& GetPx() const { return fPx; }
    inline const Double_t& GetPy() const { return fPy; }
    inline const Double_t& GetPz() const { return fPz; }
    inline const Double_t& GetQ() const { return fQ; }
    inline const Double_t& GetAoZ() const { return fAoZ; }
    inline const Double_t& GetChix() const { return fChix; }
    inline const Double_t& GetChiy() const { return fChiy; }
    inline const Int_t& GetQuality() const { return fQuality; }

  protected:
    Double_t fX;
    Double_t fY;
    Double_t fZ;
    Double_t fPx;
    Double_t fPy;
    Double_t fPz;
    Double_t fQ;
    Double_t fAoZ;
    Double_t fChix;
    Double_t fChiy;
    Int_t fQuality;

  public:
    ClassDef(R3BTrack, 1)
};

#endif
