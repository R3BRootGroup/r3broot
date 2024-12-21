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

#pragma once

#include <Rtypes.h>
#include <TObject.h>
#include <TString.h>

class R3BCoarseTimeStitch : public TObject
{
  public:
    // Default Constructor for sync two tamex channels
    R3BCoarseTimeStitch();

    // Destructor
    ~R3BCoarseTimeStitch() override = default;

    // Accessors with [[nodiscard]]
    [[nodiscard]] inline double GetRange1() const { return fRange1; }
    [[nodiscard]] inline double GetRange2() const { return fRange2; }
    [[nodiscard]] double GetTime(double, TString const& name1 = "tamex", TString const& name2 = "tamex");

    inline void SetRange1(Float_t range) { fRange1 = range; }
    inline void SetRange2(Float_t range) { fRange2 = range; }

    inline void SetClockTDC150() { fRangeClockTDC = fRangeClockTDC150; }

  private:
    double fRange1;
    double fRange2;
    double fRangeTamex;
    double fRangeVftx;
    double fRangeTrb;
    double fRangeClockTDC;
    double fRangeClockTDC150;

  public:
    ClassDefOverride(R3BCoarseTimeStitch, 1)
};
