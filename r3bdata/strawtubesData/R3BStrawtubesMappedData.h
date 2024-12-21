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

#ifndef R3BSTRAWTUBESMAPPEDITEM_H
#define R3BSTRAWTUBESMAPPEDITEM_H

#include "TObject.h"

class R3BStrawtubesMappedData : public TObject
{
  public:
    R3BStrawtubesMappedData(UChar_t, UChar_t, UChar_t, UInt_t, UInt_t);
    virtual ~R3BStrawtubesMappedData();
    UChar_t GetPlane() const { return fPlane; }
    UChar_t GetSide() const { return fSide; }
    UChar_t GetStraw() const { return fStraw; }
    UInt_t GetTimeCoarse() const { return fTimeCoarse; }
    UInt_t GetTimeFine() const { return fTimeFine; }

  private:
    UChar_t fPlane;
    UChar_t fStraw;
    UChar_t fSide;
    UInt_t fTimeCoarse;
    UInt_t fTimeFine;

  public:
    ClassDef(R3BStrawtubesMappedData, 0)
};

#endif
