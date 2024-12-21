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

#ifndef R3BNEULANDHIT_H
#define R3BNEULANDHIT_H

#include "TObject.h"
#include "TVector3.h"
#include <iostream>

class R3BNeulandHit : public TObject
{
  public:
    R3BNeulandHit() = default;
    R3BNeulandHit(Int_t paddle,
                  Double_t TdcL,
                  Double_t TdcR,
                  Double_t time,
                  Double_t QdcL,
                  Double_t QdcR,
                  Double_t energy,
                  const TVector3& pos,
                  const TVector3& pix);

    Int_t GetPaddle() const { return fPaddle; }
    Double_t GetTdcR() const { return fTdcR; }
    Double_t GetTdcL() const { return fTdcL; }
    Double_t GetT() const { return fT; }
    Double_t GetQdcR() const { return fQdcR; }
    Double_t GetQdcL() const { return fQdcL; }
    Double_t GetE() const { return fE; }
    TVector3 GetPosition() const { return fPosition; }
    TVector3 GetPixel() const { return fPixel; }

    Double_t GetBeta() const;
    Double_t GetEToF(Double_t mass = 939.565) const; // 939.565379

    void Print(const Option_t*) const override;

    bool operator==(const R3BNeulandHit& b) const
    {
        // TODO: Change this if multi-hit capability is introduced
        return this->GetPaddle() == b.GetPaddle();
    }

  private:
    Int_t fPaddle;
    Double_t fTdcL;
    Double_t fTdcR;
    Double_t fT;
    Double_t fQdcL;
    Double_t fQdcR;
    Double_t fE;
    TVector3 fPosition;
    TVector3 fPixel;

    ClassDefOverride(R3BNeulandHit, 1)
};

std::ostream& operator<<(std::ostream&, const R3BNeulandHit&); // Support easy printing

#endif // R3BNEULANDHIT_H
