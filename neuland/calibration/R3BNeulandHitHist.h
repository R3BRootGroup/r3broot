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

// -----------------------------------------------------------------------------
// -----                            R3BNeulandHitHist                      -----
// -----                     Created 22-04-2014 by D.Kresan                -----
// -----------------------------------------------------------------------------

#ifndef R3BNEULANDHITHIST_H
#define R3BNEULANDHITHIST_H

#include "FairTask.h"

class TClonesArray;
class TH1F;
class TH2F;

class R3BNeulandHitHist : public FairTask
{
  public:
    R3BNeulandHitHist();
    R3BNeulandHitHist(const char* name, Int_t iVerbose);
    virtual ~R3BNeulandHitHist();

    virtual InitStatus Init();

    virtual void Exec(Option_t* option);

    virtual void FinishTask();

    inline void SetNofBars(Int_t nofBars) { fNofBars = nofBars; }

    inline void SetFirstPlaneHorisontal() { fFirstPlaneHorisontal = kTRUE; }

    inline void SetTimeOffset(Double_t t) { fTimeOffset = t; }

    inline void SetDistance(Double_t t) { fDistance = t; }

  private:
    Int_t fNofBars;
    Double_t fTimeOffset;
    Bool_t fFirstPlaneHorisontal;

    Int_t fnEvents;
    Double_t fDistance;

    TClonesArray* fLandDigi;
    TClonesArray* fLosHit;

    TH1F* fh_land_barid;
    TH2F* fh_land_qdcbarid;
    TH1F* fh_land_tof;
    TH2F* fh_land_qdctof;
    TH2F* fh_land_qdctof_norm;
    TH2F* fh_land_timebarid;
    TH2F* fh_land_tofbarid;
    TH2F* fh_land_betabarid;
    TH2F* fh_land_yx;
    TH2F* fh_land_yx1;
    TH2F* fh_land_yx2;
    TH2F* fh_land_lbarid;
    TH2F* fh_land_ltime;

    TH1F* fh_land_beta;
    TH1F* fh_land_qdc;
    TH1F* fh_land_qdc_cut;

    TH1F* fh_los_time;

    TH2F* fh_los_corr;

    TH1F* fh_land_norm_tof;

    void CreateHistos();

    void WriteHistos();

  public:
    ClassDef(R3BNeulandHitHist, 0)
};

#endif
