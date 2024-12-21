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

#ifndef R3BSTRAWTUBESMAPPED2CALPAR_H
#define R3BSTRAWTUBESMAPPED2CALPAR_H

#include "FairTask.h"

class R3BTCalEngine;
class R3BTCalPar;
class TClonesArray;

class R3BStrawtubesMapped2CalPar : public FairTask
{
  public:
    R3BStrawtubesMapped2CalPar(const char* = "R3BStrawtubesMapped2CalPar", Int_t = 1);
    virtual ~R3BStrawtubesMapped2CalPar();
    virtual void Exec(Option_t* option);
    virtual void FinishEvent();
    virtual void FinishTask();
    virtual InitStatus Init();
    void SetMinStats(Int_t);

  private:
    Int_t fMinStats;
    R3BTCalPar* fCalPar;
    R3BTCalEngine* fEngine;
    TClonesArray* fMapped;

  public:
    ClassDef(R3BStrawtubesMapped2CalPar, 0)
};

#endif
