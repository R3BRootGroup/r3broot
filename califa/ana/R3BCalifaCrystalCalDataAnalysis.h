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

#ifndef R3BCALIFACRYSTALCALDATAANALYSIS_H
#define R3BCALIFACRYSTALCALDATAANALYSIS_H

#include "FairTask.h"

class TClonesArray;
class TH1F;
class TH2F;

class R3BCalifaCrystalCalDataAnalysis : public FairTask
{
  public:
    R3BCalifaCrystalCalDataAnalysis();
    virtual ~R3BCalifaCrystalCalDataAnalysis();

    virtual InitStatus Init();

    virtual void Exec(Option_t* option);

    virtual void FinishTask();

  private:
    Int_t fnEvents;

    TClonesArray* fCrystalData;

    TH1F* thCrystalID;
    TH1F* thEnergy;
    TH1F* thNf;
    TH1F* thNs;
    TH1F* thTime;

    void CreateHistos();

    void WriteHistos();

  public:
    ClassDef(R3BCalifaCrystalCalDataAnalysis, 1)
};

#endif
