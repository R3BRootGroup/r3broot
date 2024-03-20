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
#pragma once

#include "TObject.h"

class R3BFrsSciPosCalData : public TObject
{
  public:
    // Default Constructor
    R3BFrsSciPosCalData();

    // Standard Constructor
    R3BFrsSciPosCalData(UShort_t, Double_t, Double_t, Float_t, Float_t);

    // Destructor
    virtual ~R3BFrsSciPosCalData() {}

    // Getters
    inline const UShort_t& GetDetector() const { return fDetector; }
    inline const Double_t& GetRawTimeNs() const { return fRawTimeNs; }
    inline const Double_t& GetRawTimeNsWithTref() const { return fRawTimeNs_wTref; }
    inline const Float_t& GetRawPosNs() const { return fRawPosNs; }
    inline const Float_t& GetCalPosMm() const { return fCalPosMm; }

  private:
    UShort_t fDetector;
    Double_t fRawTimeNs;       // 0.5 * (Tright + Tleft) after selection of the multiplicity
    Double_t fRawTimeNs_wTref; // 0.5 * (Tright + Tleft) - Ttrig
    Float_t fRawPosNs;         // Tright - Tleft: x increasing for right to left
    Float_t fCalPosMm;         // calibrated position in Mm

  public:
    ClassDef(R3BFrsSciPosCalData, 3)
};
