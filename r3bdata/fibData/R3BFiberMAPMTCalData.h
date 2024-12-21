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

// Adapted from R3BunchedFiber class by V. Panin, May 2021
// This class contains calibrated data for single MAPMT channels of
// R3BFiberMAPMT detectors i.e. fib30,31,32,33,23A(B)

#ifndef R3BFIBERMAPMTCALDATA_H
#define R3BFIBERMAPMTCALDATA_H

#include "TObject.h"

class R3BFiberMAPMTCalData : public TObject
{
  public:
    R3BFiberMAPMTCalData();
    R3BFiberMAPMTCalData(Int_t, Int_t, Bool_t, Double_t);
    virtual ~R3BFiberMAPMTCalData();

    Int_t GetChannel() const;
    Double_t GetTime_ns() const;
    Int_t GetSide() const;
    Double_t GetTimeLeading_ns() const;
    Double_t GetTimeTrailing_ns() const;
    Bool_t IsBottom() const;
    Bool_t IsTop() const;
    Bool_t IsTrigger() const;
    Bool_t IsLeading() const;
    Bool_t IsTrailing() const;
    Bool_t IsSortable() const;

  private:
    Int_t fSide; // 0 = bottom, 1 = top, 2 = MAPMT trigger
    Int_t fChannel;
    Bool_t fIsLeading;
    Double_t fTime_ns;

    ClassDef(R3BFiberMAPMTCalData, 1)
};

#endif
