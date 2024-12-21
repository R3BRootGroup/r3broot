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

#ifndef R3BNEULANDNEUTRON_H
#define R3BNEULANDNEUTRON_H

#include "R3BNeulandCluster.h"
#include "R3BNeulandHit.h"
#include "TObject.h"
#include "TVector3.h"
#include <iostream>

class R3BNeulandNeutron : public TObject
{
  public:
    R3BNeulandNeutron() = default;
    R3BNeulandNeutron(Int_t paddle, Double_t time, const TVector3& pos, const TVector3& pix);
    explicit R3BNeulandNeutron(const R3BNeulandHit&);
    explicit R3BNeulandNeutron(const R3BNeulandCluster&);

    Int_t GetPaddle() const { return fPaddle; }
    Double_t GetT() const { return fT; }
    TVector3 GetPosition() const { return fPosition; }
    TVector3 GetPixel() const { return fPixel; }

    Double_t GetEtot() const;
    Double_t GetEkin() const;
    Double_t GetGamma() const;
    TVector3 GetP() const;

    void Print(const Option_t*) const override;

  private:
    Int_t fPaddle;
    Double_t fT;
    TVector3 fPosition;
    TVector3 fPixel;

    ClassDefOverride(R3BNeulandNeutron, 1)
};

std::ostream& operator<<(std::ostream&, const R3BNeulandNeutron&); // Support easy printing

#endif // R3BNEULANDNEUTRON_H
