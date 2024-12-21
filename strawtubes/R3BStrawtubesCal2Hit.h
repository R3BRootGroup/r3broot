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

#ifndef R3BSTRAWTUBESCAL2HIT
#define R3BSTRAWTUBESCAL2HIT

#include "FairTask.h"

class TClonesArray;

class R3BStrawtubesCal2Hit : public FairTask
{
  public:
    R3BStrawtubesCal2Hit(const char* = "StrawtubesCal2Hit", Int_t = 1);
    virtual ~R3BStrawtubesCal2Hit();
    virtual void Exec(Option_t*);
    virtual void FinishEvent();
    virtual void FinishTask();
    virtual InitStatus Init();
    virtual InitStatus ReInit();

  private:
    TClonesArray* fCalItems;
    TClonesArray* fHitItems;
    UInt_t fNofHitItems;

  public:
    ClassDef(R3BStrawtubesCal2Hit, 0)
};

#endif
