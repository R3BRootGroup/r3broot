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
// -----                        R3BPrimPart                                -----
// -----                 Created 12-07-2011 by D.Kresan                    -----
// -----------------------------------------------------------------------------

#ifndef R3BPRIMPART
#define R3BPRIMPART

#include "TObject.h"
#include "TVector3.h"

class R3BPrimPart : public TObject
{
  public:
    R3BPrimPart();
    R3BPrimPart(Int_t pdg,
                Double_t px,
                Double_t py,
                Double_t pz,
                Double_t x,
                Double_t y,
                Double_t z,
                Double_t t,
                Double_t A,
                Double_t M);
    virtual ~R3BPrimPart();

    Int_t GetPdgCode() const { return fPdg; }
    void Momentum(TVector3& mom) const { mom.SetXYZ(fPx, fPy, fPz); }
    void Position(TVector3& pos) const { pos.SetXYZ(fX, fY, fZ); }
    Double_t GetT() const { return fT; }
    Double_t GetA() const { return fA; }
    Double_t GetM() const { return fM; }
    Double_t GetE() const
    {
        TVector3 mom;
        Momentum(mom);
        return sqrt(mom.Mag2() + GetM2());
    }
    Double_t GetM2() const { return fM * fM; }
    Double_t GetBeta() const
    {
        TVector3 mom;
        Momentum(mom);
        return mom.Mag() / GetE();
    }
    Double_t GetGamma() const { return 1. / sqrt(1. - pow(GetBeta(), 2)); }

    void SetPz(const Double_t& pz) { fPz = pz; }
    void SetA(const Double_t& a) { fA = a; }
    void SetM(const Double_t& m) { fM = m; }

  private:
    Int_t fPdg;
    Double_t fPx;
    Double_t fPy;
    Double_t fPz;
    Double_t fX;
    Double_t fY;
    Double_t fZ;
    Double_t fT;
    Double_t fA;
    Double_t fM;

  public:
    ClassDef(R3BPrimPart, 1)
};

#endif
