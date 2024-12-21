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

#ifndef R3BNEULANDCALTEST_H
#define R3BNEULANDCALTEST_H

#include "FairTask.h"

class TClonesArray;
class TH1F;
class TH2F;

class R3BNeulandCalTest : public FairTask
{
  public:
    R3BNeulandCalTest();
    R3BNeulandCalTest(const char* name, Int_t iVerbose);
    virtual ~R3BNeulandCalTest();

    virtual InitStatus Init();

    virtual void Exec(Option_t* option);

    virtual void FinishTask();

  private:
    Int_t fnEvents;

    TClonesArray* fLandPmt;
    TClonesArray* fNeulandPmt;
    TClonesArray* fLosHit;

    TH2F* fh_los_corr;
    TH1F* fh_tof;
    TH2F* fh_qdctof;
    TH1F* fh_tdiff;

    TH2F* fh_qdctof_2;

    void CreateHistos();

    void WriteHistos();

  public:
    ClassDef(R3BNeulandCalTest, 1)
};

#endif
