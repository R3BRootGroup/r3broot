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

#ifndef R3BTOFDMAPPEDITEM_H
#define R3BTOFDMAPPEDITEM_H

#include "TObject.h"

class R3BTofdMappedData : public TObject
{
  public:
    R3BTofdMappedData();
    R3BTofdMappedData(UInt_t, UInt_t, UInt_t, UInt_t, UInt_t, UInt_t);

    UInt_t GetDetectorId() const;
    UInt_t GetSideId() const;
    UInt_t GetBarId() const;
    UInt_t GetEdgeId() const;
    UInt_t GetTimeCoarse() const;
    UInt_t GetTimeFine() const;

  private:
    UInt_t fDetector; // 1..n
    UInt_t fSide;     // 1 = bottom, 2 = top
    UInt_t fBar;      // 1..n
    UInt_t fEdge;     // 1 = leading, 2 = trailing
    UInt_t fTimeCoarse;
    UInt_t fTimeFine;

  public:
    ClassDef(R3BTofdMappedData, 1)
};

#endif
