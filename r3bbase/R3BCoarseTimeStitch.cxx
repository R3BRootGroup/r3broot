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

#include "R3BCoarseTimeStitch.h"
#include "R3BLogger.h"

#include <TMath.h>
#define IS_NAN(x) TMath::IsNaN(x)

// Standard constructur
R3BCoarseTimeStitch::R3BCoarseTimeStitch()
    : fRange1(2048 * 1000. / 200.)           // ns
    , fRange2(2048 * 1000. / 200.)           // ns
    , fRangeTamex(2048 * 1000. / 200.)       // ns
    , fRangeTrb(2048 * 1000. / 200.)         // ns
    , fRangeVftx(8192 * 1000. / 200.)        // ns
    , fRangeClockTDC(4096 * 1000. / 250.)    // ns
    , fRangeClockTDC150(4096 * 1000. / 150.) // ns
{
}

double R3BCoarseTimeStitch::GetTime(double time, TString const& name1, TString const& name2)
{
    if (name1 == "tamex")
    {
        fRange1 = fRangeTamex;
    }
    else if (name1 == "trb")
    {
        fRange1 = fRangeTrb;
    }
    else if (name1 == "vftx")
    {
        fRange1 = fRangeVftx;
    }
    else if (name1 == "clocktdc")
    {
        fRange1 = fRangeClockTDC;
    }
    else
    {
        R3BLOG(fatal, "Module " << name1 << " does not exist.");
    }

    if (name2 == "tamex")
    {
        fRange2 = fRangeTamex;
    }
    else if (name2 == "trb")
    {
        fRange2 = fRangeTrb;
    }
    else if (name2 == "vftx")
    {
        fRange2 = fRangeVftx;
    }
    else if (name2 == "clocktdc")
    {
        fRange2 = fRangeClockTDC;
    }
    else
    {
        R3BLOG(fatal, "Module " << name2 << " does not exist.");
    }

    const double const1 = TMath::Min(fRange1, fRange2);
    const double const2 = TMath::Max(fRange1, fRange2);

    R3BLOG(debug,
           "Time: " << time << " , range1(" << name1 << "): " << fRange1 << " , range2(" << name1 << "): " << fRange2);
    return fmod(time + const2 + const1 / 2., const1) - const1 / 2.;
}

ClassImp(R3BCoarseTimeStitch)
