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

#include "R3BLosCalData.h"
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
    double const c_range_VFTX = 5 * 8192;
} // namespace
/* for the data analysis of the Los detector.
 * Introduced by Ralf, Feb 2016
 *
 * March 2016
 * Rewrote the Cal structure to provide individual leafs for the
 * left, top, right and bottom signals. This allows to plot
 * the time differences via cbmsim->Draw(...) interactively (aka without
 * looping over all channels) which is crucial for a quick check of the
 * detector status during the experiment.
 *
 * September 2016
 * TAMEX readout of MCFD amplitudes added by Aleksandra.
 *
 * November 2016
 * Data structure changed and adapted to new LOS detector (with 8 PMs)
 *
 * Current version since Dez 2018 by Hans Toernqvist
 *
 */

R3BLosCalData::R3BLosCalData()
    : fDetector(-1)
{
    Reset();
}
R3BLosCalData::R3BLosCalData(UInt_t detector)
    : fDetector(detector)
{
    Reset();
}

UInt_t R3BLosCalData::GetDetector() const { return fDetector; }

void R3BLosCalData::Reset()
{
    for (Int_t j = 0; j < 8; ++j)
    {
        fTimeV_ns[j] = 0.0 / 0.0;
        fTimeL_ns[j] = 0.0 / 0.0;
        fTimeT_ns[j] = 0.0 / 0.0;
        fTimeM_ns[j] = 0.0 / 0.0;
    }
}

Double_t R3BLosCalData::GetTimeV_ns(UInt_t chann) const
{
    ASSERT(chann, <, LENGTH(fTimeV_ns));
    return fTimeV_ns[chann];
}

Double_t R3BLosCalData::GetTimeL_ns(UInt_t chann) const
{
    ASSERT(chann, <, LENGTH(fTimeL_ns));
    return fTimeL_ns[chann];
}

Double_t R3BLosCalData::GetTimeT_ns(UInt_t chann) const
{
    ASSERT(chann, <, LENGTH(fTimeT_ns));
    return fTimeT_ns[chann];
}

Double_t R3BLosCalData::GetTimeM_ns(UInt_t chann) const
{
    ASSERT(chann, <, LENGTH(fTimeM_ns));
    return fTimeM_ns[chann];
}

/*
Double_t R3BLosCalData::GetTime_ns(UInt_t chann, UInt_t i) const
{
  ASSERT(i, <, fTime_length);
  return fTime_ns[chann][i];
}

Double_t R3BLosCalData::GetTimeDiff(Double_t a_early, Double_t a_late, Double_t c_range) const
{
  return fmod(a_late - a_early + c_range, c_range);
}

void R3BLosCalData::Reset()
{
 for (Int_t j = 0; j < sizeof(fTime_ns)/sizeof(fTime_ns[0]); ++j){
  for (Int_t i = 0; i < sizeof(fTime_ns[0])/sizeof(fTime_ns[0][0]); ++i) {
    fTime_ns[j][i] = {0.0 / 0.0};
   }
 }
  fFirst = -1;
  fLast = -1;
}

Bool_t R3BLosCalData::SetTimeV_ns(UInt_t chann, Double_t time_ns, Double_t time_width)
{
  return SetTime_ns(0, chann, time_ns, time_width, c_range_VFTX);
}

Bool_t R3BLosCalData::SetTimeL_ns(UInt_t chann, Double_t time_ns, Double_t time_width)
{
  return SetTime_ns(1, chann, time_ns, time_width, c_range_TAMEX);
}

Bool_t R3BLosCalData::SetTimeT_ns(UInt_t chann, Double_t time_ns, Double_t time_width)
{
  return SetTime_ns(2, chann, time_ns, time_width, c_range_TAMEX);
}

Bool_t R3BLosCalData::SetTime_ns(UInt_t i, UInt_t icha, Double_t time_ns, Double_t
    time_width, Double_t c_range)
{

 //std::cout<<"SET TIME START: "<<i<<", "<<icha<<", "<<", "<<time_ns<<", "<<c_range<<", "<<fTime_ns[icha][i]<<std::endl;

  ASSERT(i*icha, <, fTime_length);//LENGTH(fTime_ns));
  if (time_ns >= c_range) {
     std::cout<<"Time for channel "<<icha<<", type "<<i<<" out of range "<<time_ns<<std::endl;
      return kFALSE;
  }
  ASSERT(time_ns, <=, c_range);

  if (!IS_NAN(fTime_ns[icha][i])) {

     std::cout<<"Time for channel "<<icha<<", type "<<i<<" already set: "<<time_ns<<std::endl;

    return kFALSE;
  }

  if (fFirst < 0) {
    // Nothing set, just grab this one.
    fLast = fFirst = time_ns;

    std::cout<<"Nothing yet set; start new item with time for channel "<<icha<<", type "<<i<<" time:
"<<time_ns<<std::endl; } else { auto fFirst_time_ns = GetTimeDiff(fFirst, time_ns, c_range); auto time_ns_fLast =
GetTimeDiff(time_ns, fLast, c_range);


    if (fFirst_time_ns < time_width) {
      if (time_ns_fLast > 0) {
        fLast = time_ns;

      }
    } else if (time_ns_fLast < time_width) {
      if (fFirst_time_ns > 0) {
        fFirst = time_ns;


      }
    } else {

    std::cout<<"Time for channel "<<icha<<", type "<<i<<" not compatible"<<std::endl;


      // New time is way off previous times, incompatible.
      return kFALSE;
    }
  }

//  fTime_ns[icha][i] = time_ns;

//   std::cout<<"SET TIME EXIT: "<<i<<", "<<time_ns<<", "<<fTime_ns[icha][i]<<", "<<fFirst<<", "<<fLast<<std::endl;

  return kTRUE;
}
*/

Double_t R3BLosCalData::GetMeanTime()
{
    Double_t mean = 0;
    Int_t num_times = 0;
    for (Int_t i = 0; i < 8; i++)
    {
        if (!std::isnan(fTimeV_ns[i]))
        {
            mean += fTimeV_ns[i];
            num_times++;
        }
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
Double_t R3BLosCalData::GetMeanTimeTAMEXL()
{
    Double_t mean = 0;
    Int_t num_times = 0;
    for (Int_t i = 0; i < 8; i++)
    {
        if (!std::isnan(fTimeL_ns[i]))
        {
            mean += fTimeL_ns[i];
            num_times++;
        }
    }
    return num_times ? mean / num_times : NAN;
}

Double_t R3BLosCalData::GetMeanTimeTAMEXT()
{
    Double_t mean = 0;
    Int_t num_times = 0;
    for (Int_t i = 0; i < 8; i++)
    {
        if (!std::isnan(fTimeT_ns[i]))
        {
            mean += fTimeT_ns[i];
            num_times++;
        }
    }

    return num_times ? mean / num_times : NAN;
}

UInt_t R3BLosCalData::GetTAMEXLNcha()
{
    Int_t num_times = 0;
    for (Int_t i = 0; i < 8; i++)
    {
        if (!std::isnan(fTimeL_ns[i]))
        {
            num_times++;
        }
    }

    return num_times;
}

UInt_t R3BLosCalData::GetTAMEXTNcha()
{
    Int_t num_times = 0;
    for (Int_t i = 0; i < 8; i++)
    {
        if (!std::isnan(fTimeT_ns[i]))
        {
            num_times++;
        }
    }

    return num_times;
}

Double_t R3BLosCalData::GetMeanTimeVFTX()
{
    Double_t mean = 0;
    Int_t num_times = 0;
    for (Int_t i = 0; i < 8; i++)
    {
        if (!std::isnan(fTimeV_ns[i]))
        {
            mean += fTimeV_ns[i];
            num_times++;
        }
    }

    return num_times ? mean / num_times : NAN;
}

Double_t R3BLosCalData::GetMeanTimeMTDC32()
{
    Double_t mean = 0;
    Int_t num_times = 0;
    for (Int_t i = 0; i < 8; i++)
    {
        if (!std::isnan(fTimeM_ns[i]))
        {
            mean += fTimeM_ns[i];
            num_times++;
        }
    }
    return num_times ? mean / num_times : NAN;
}

UInt_t R3BLosCalData::GetVFTXNcha()
{
    Int_t num_times = 0;
    for (Int_t i = 0; i < 8; i++)
    {
        if (!std::isnan(fTimeV_ns[i]))
        {
            num_times++;
        }
    }

    return num_times;
}
UInt_t R3BLosCalData::GetMTDC32Ncha()
{
    Int_t num_times = 0;
    for (Int_t i = 0; i < 8; i++)
    {
        if (!std::isnan(fTimeM_ns[i]))
        {
            num_times++;
        }
    }
    return num_times;
}

ClassImp(R3BLosCalData)
