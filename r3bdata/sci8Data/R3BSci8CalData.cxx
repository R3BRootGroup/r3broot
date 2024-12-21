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

#include "R3BSci8CalData.h"
#include "FairLogger.h"
#include <cmath>
#include <iostream>

/* Re isnan:
 * std::isnan is defined as an overloaded function in cmath
 * the macro isnan is defined in math.h
 * if you include math.h and do using namespace std;, my
 * compiler will complain about ambiguous definition of isnan
 * because cmath seems to get pulled in as a dependency.
 *
 * Solution: we explicitly include cmath and do not commit
 * what has been described as the
 * "environmental crime of namespace pollution",
 * in case anyone includes math.h for some reason.
 */

/* for the data analysis of Sci8 detector.
 */

R3BSci8CalData::R3BSci8CalData()
    : fDetector(0)
{
}

R3BSci8CalData::R3BSci8CalData(UInt_t detector)
    : fDetector(detector)
    , fTimeV_r_ns(0.0 / 0.0) // 0/0 produces a NAN. Macro NAN not available??
    , fTimeV_l_ns(0.0 / 0.0)
    , fTimeL_r_ns(0.0 / 0.0)
    , fTimeL_l_ns(0.0 / 0.0)
    , fTimeT_r_ns(0.0 / 0.0)
    , fTimeT_l_ns(0.0 / 0.0)
{
}

Double_t R3BSci8CalData::GetMeanTime()
{
    Double_t mean = 0;
    Int_t num_times = 0;
    if (!std::isnan(fTimeV_r_ns))
    {
        mean += fTimeV_r_ns;
        num_times++;
    }
    if (!std::isnan(fTimeV_l_ns))
    {
        mean += fTimeV_l_ns;
        num_times++;
    }
    if (!std::isnan(fTimeL_r_ns))
    {
        mean += fTimeL_r_ns;
        num_times++;
    }
    if (!std::isnan(fTimeL_l_ns))
    {
        mean += fTimeL_l_ns;
        num_times++;
    }
    if (!std::isnan(fTimeT_r_ns))
    {
        mean += fTimeT_r_ns;
        num_times++;
    }
    if (!std::isnan(fTimeT_l_ns))
    {
        mean += fTimeT_l_ns;
        num_times++;
    }

    return num_times ? mean / num_times : NAN;
}
Double_t R3BSci8CalData::GetMeanTimeTAMEXL()
{
    Double_t mean = 0;
    Int_t num_times = 0;
    if (!std::isnan(fTimeL_r_ns))
    {
        mean += fTimeL_r_ns;
        num_times++;
    }
    if (!std::isnan(fTimeL_l_ns))
    {
        mean += fTimeL_l_ns;
        num_times++;
    }

    return num_times ? mean / num_times : NAN;
}

Double_t R3BSci8CalData::GetMeanTimeTAMEXT()
{
    Double_t mean = 0;
    Int_t num_times = 0;
    if (!std::isnan(fTimeT_r_ns))
    {
        mean += fTimeT_r_ns;
        num_times++;
    }
    if (!std::isnan(fTimeT_l_ns))
    {
        mean += fTimeT_l_ns;
        num_times++;
    }

    return num_times ? mean / num_times : NAN;
}

UInt_t R3BSci8CalData::GetTAMEXLNcha()
{
    Int_t num_times = 0;
    if (!std::isnan(fTimeL_r_ns))
    {
        num_times++;
    }
    if (!std::isnan(fTimeL_l_ns))
    {
        num_times++;
    }

    return num_times;
}

UInt_t R3BSci8CalData::GetTAMEXTNcha()
{
    Int_t num_times = 0;
    if (!std::isnan(fTimeT_r_ns))
    {
        num_times++;
    }
    if (!std::isnan(fTimeT_l_ns))
    {
        num_times++;
    }

    return num_times;
}

Double_t R3BSci8CalData::GetMeanTimeVFTX()
{
    Double_t mean = 0;
    Int_t num_times = 0;
    if (!std::isnan(fTimeV_r_ns))
    {
        mean += fTimeV_r_ns;
        num_times++;
    }
    if (!std::isnan(fTimeV_l_ns))
    {
        mean += fTimeV_l_ns;
        num_times++;
    }

    return num_times ? mean / num_times : NAN;
}
UInt_t R3BSci8CalData::GetVFTXNcha()
{
    Int_t num_times = 0;
    if (!std::isnan(fTimeV_r_ns))
    {
        num_times++;
    }
    if (!std::isnan(fTimeV_l_ns))
    {
        num_times++;
    }

    return num_times;
}
ClassImp(R3BSci8CalData)
