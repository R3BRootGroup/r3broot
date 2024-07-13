/******************************************************************************
 *   Copyright (C) 2019 GSI Helmholtzzentrum für Schwerionenforschung GmbH    *
 *   Copyright (C) 2019-2024 Members of R3B Collaboration                     *
 *                                                                            *
 *             This software is distributed under the terms of the            *
 *                 GNU General Public Licence (GPL) version 3,                *
 *                    copied verbatim in the file "LICENSE".                  *
 *                                                                            *
 * In applying this license GSI does not waive the privileges and immunities  *
 * granted to it by virtue of its status as an Intergovernmental Organization *
 * or submit itself to any jurisdiction.                                      *
 ******************************************************************************/

#ifndef R3BNEULANDCALDATA_H
#define R3BNEULANDCALDATA_H

#include "TObject.h"
#include <iostream>

class R3BNeulandCalData : public TObject
{
  public:
    R3BNeulandCalData() = default;
    R3BNeulandCalData(Int_t barId, Int_t side, Double_t time, Double_t triggertime, double qdc);

    Int_t GetBarId() const { return fBarId; }
    Int_t GetSide() const { return fSide; }
    Double_t GetTime() const { return fTime; }
    Double_t GetTriggerTime() const { return fTriggerTime; }
    [[nodiscard]] auto GetQdc() const -> double { return fQdc; }

    void Print(const Option_t*) const override;

  private:
    Int_t fBarId;
    Int_t fSide;
    Double_t fTime;
    Double_t fTriggerTime;
    double fQdc;

    ClassDefOverride(R3BNeulandCalData, 3)
};

std::ostream& operator<<(std::ostream&, const R3BNeulandCalData&); // Support easy printing

#endif // R3BNEULANDCALDATA_H
