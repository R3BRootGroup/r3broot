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

// ----------------------------------------------------------------
// -----              R3BBunchedFiberMappedData                -----
// -----         Created May 2021 by V. Panin and H. Toernqvist         -----
// ----------------------------------------------------------------
//
// This class contains mapped data for one fiber electronics channel and
// logical signal edge, not for one fiber!
//

#ifndef R3BFIBERMAPMTMAPPEDDATA_H
#define R3BFIBERMAPMTMAPPEDDATA_H

#include "TObject.h"

class R3BFiberMAPMTMappedData : public TObject
{
  public:
    R3BFiberMAPMTMappedData();
    R3BFiberMAPMTMappedData(Int_t, Int_t, Bool_t, Int_t, Int_t);
    virtual ~R3BFiberMAPMTMappedData();

    Int_t GetChannel() const;
    Int_t GetCoarse() const;
    Int_t GetFine() const;
    Int_t GetSide() const;
    Bool_t IsTop() const;
    Bool_t IsBottom() const;
    Bool_t IsTrigger() const;
    Bool_t IsLeading() const;
    Bool_t IsTrailing() const;

  protected:
    Int_t fSide; // 0 = bottom, 1 = top, 2 = MAPMT trigger
    Int_t fChannel;
    Bool_t fIsLeading;
    Int_t fCoarse;
    Int_t fFine;

  public:
    ClassDef(R3BFiberMAPMTMappedData, 2)
};

#endif
