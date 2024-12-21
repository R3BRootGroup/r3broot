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

#include "R3BRoluCalData.h"
#include "FairLogger.h"
#include "TMath.h"
#include <cmath>
#include <iostream>

#define ASSERT(a, op, b)                                                                                     \
    do                                                                                                       \
    {                                                                                                        \
        auto const a_ = a;                                                                                   \
        auto const b_ = b;                                                                                   \
        if (!(a_ op b_))                                                                                     \
        {                                                                                                    \
            std::cerr << "Assertion \"" #a << '=' << a_ << ' ' << #op << " " #b << '=' << b_ << "\" failed." \
                      << std::endl;                                                                          \
            abort();                                                                                         \
        }                                                                                                    \
    } while (0)
#define IS_NAN(x) TMath::IsNaN(x)
#define LENGTH(x) (sizeof x / sizeof x[0])

namespace
{
    double const c_range_TAMEX = 5 * 2048;
}
/*
 * for the data analysis of the Rolu detector following LOS
 */

R3BRoluCalData::R3BRoluCalData()
    : fDetector(0)

{
    Reset();
}
R3BRoluCalData::R3BRoluCalData(UInt_t detector)
    : fDetector(detector)
{
    Reset();
}

UInt_t R3BRoluCalData::GetDetector() const { return fDetector; }

void R3BRoluCalData::Reset()
{
    for (Int_t j = 0; j < 4; ++j)
    {
        fTimeL_ns[j] = 0.0 / 0.0;
        fTimeT_ns[j] = 0.0 / 0.0;
    }
}

Double_t R3BRoluCalData::GetTimeL_ns(UInt_t chann) const
{
    ASSERT(chann, <, LENGTH(fTimeL_ns));
    return fTimeL_ns[chann];
}

Double_t R3BRoluCalData::GetTimeT_ns(UInt_t chann) const
{
    ASSERT(chann, <, LENGTH(fTimeT_ns));
    return fTimeT_ns[chann];
    ;
}

Double_t R3BRoluCalData::GetMeanTime()
{
    Double_t mean = 0;
    Int_t num_times = 0;
    for (Int_t i = 0; i < 4; i++)
    {
        if (!std::isnan(fTimeL_ns[i]))
        {
            mean += fTimeL_ns[i];
            num_times++;
        }
        if (!std::isnan(fTimeT_ns[i]))
        {
            mean += fTimeT_ns[i];
            num_times++;
        }
    }

    return num_times ? mean / num_times : NAN;
}
Double_t R3BRoluCalData::GetMeanTimeTAMEXL()
{
    Double_t mean = 0;
    Int_t num_times = 0;
    for (Int_t i = 0; i < 4; i++)
    {
        if (!std::isnan(fTimeL_ns[i]))
        {
            mean += fTimeL_ns[i];
            num_times++;
        }
    }
    return num_times ? mean / num_times : NAN;
}

Double_t R3BRoluCalData::GetMeanTimeTAMEXT()
{
    Double_t mean = 0;
    Int_t num_times = 0;
    for (Int_t i = 0; i < 4; i++)
    {
        if (!std::isnan(fTimeT_ns[i]))
        {
            mean += fTimeT_ns[i];
            num_times++;
        }
    }

    return num_times ? mean / num_times : NAN;
}

UInt_t R3BRoluCalData::GetTAMEXLNcha()
{
    Int_t num_times = 0;
    for (Int_t i = 0; i < 4; i++)
    {
        if (!std::isnan(fTimeL_ns[i]))
        {
            num_times++;
        }
    }

    return num_times;
}

UInt_t R3BRoluCalData::GetTAMEXTNcha()
{
    Int_t num_times = 0;
    for (Int_t i = 0; i < 4; i++)
    {
        if (!std::isnan(fTimeT_ns[i]))
        {
            num_times++;
        }
    }

    return num_times;
}

ClassImp(R3BRoluCalData);
