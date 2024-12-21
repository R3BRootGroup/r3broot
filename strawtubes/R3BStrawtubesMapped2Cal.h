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

#ifndef R3BSTRAWTUBESMAPPED2CAL
#define R3BSTRAWTUBESMAPPED2CAL

#include "FairTask.h"
#include <TClonesArray.h>

class R3BTCalPar;

class R3BStrawtubesMapped2Cal : public FairTask
{
  public:
    R3BStrawtubesMapped2Cal(const char* = "StrawtubesTcal", Int_t = 1);
    virtual ~R3BStrawtubesMapped2Cal();
    virtual void Exec(Option_t*);
    virtual void FinishEvent();
    virtual void FinishTask();
    virtual InitStatus Init();
    virtual InitStatus ReInit();
    virtual void SetParContainers();

  private:
    TClonesArray* fMappedItems;
    TClonesArray* fCalItems;
    Int_t fNofCalItems;
    R3BTCalPar* fTcalPar;

  public:
    ClassDef(R3BStrawtubesMapped2Cal, 0)
};

#endif
