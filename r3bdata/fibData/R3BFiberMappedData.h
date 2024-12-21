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

// ----------------------------------------------------------------------
// -----              R3BFiberMappedData header file                -----
// -----         Created 13/04/22 by J.L. Rodriguez-Sanchez         -----
// ----------------------------------------------------------------------

#ifndef R3BFiberMappedData_H
#define R3BFiberMappedData_H 1

#include "TObject.h"

class R3BFiberMappedData : public TObject
{
  public:
    R3BFiberMappedData();
    R3BFiberMappedData(UInt_t, UInt_t, Bool_t, Int_t, Int_t);
    virtual ~R3BFiberMappedData(){};

    inline UInt_t GetSide() const { return fSide; }
    inline UInt_t GetChannel() const { return fChannel; }
    inline Int_t GetCoarse() const { return fCoarse; }
    inline Int_t GetFine() const { return fFine; }
    inline Bool_t IsTop() const { return 2 == fSide; }
    inline Bool_t IsBottom() const { return 1 == fSide; }
    inline Bool_t IsTrigger() const { return 3 == fSide; }
    inline Bool_t IsLeading() const { return fIsLeading; }
    inline Bool_t IsTrailing() const { return !fIsLeading; }
    inline Bool_t IsMAPMT() const { return 1 == fSide; }
    inline Bool_t IsSPMT() const { return 2 == fSide; }
    inline Bool_t IsMAPMTTrigger() const { return 3 == fSide; }
    inline Bool_t IsSPMTTrigger() const { return 4 == fSide; }

  protected:
    // fSide is defined as follows:
    // 1 = bottom, 2 = top, 3 = MAPMT trigger (Fib3X)
    // 1 = MAPMT, 2 = SPMT, 3 = MAPMT trigger, 4 = SPMT trigger (Fib10, 11,...)
    UInt_t fSide;
    UInt_t fChannel;
    Bool_t fIsLeading;
    Int_t fCoarse;
    Int_t fFine;

  public:
    ClassDef(R3BFiberMappedData, 1)
};

#endif /* R3BFiberMappedData_H */
