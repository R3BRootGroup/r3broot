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

#ifndef R3BNEULANDMAPPEDHIST_H
#define R3BNEULANDMAPPEDHIST_H

#include "FairTask.h"

class R3BEventHeader;
class TClonesArray;
class TH1F;
class TH2F;

class R3BNeulandMappedHist : public FairTask
{
  public:
    R3BNeulandMappedHist();
    R3BNeulandMappedHist(const char* name, Int_t iVerbose);
    virtual ~R3BNeulandMappedHist();

    virtual InitStatus Init();

    virtual void Exec(Option_t* option);

    virtual void FinishTask();

    inline Int_t GetNItemsTotal() { return fNItemsTotal; }

  private:
    Int_t fnEvents;
    Int_t fNItemsTotal;

    R3BEventHeader* fHeader;
    TClonesArray* fLandMappedData;
    TClonesArray* fNeulandTamexHitMapped;
    TClonesArray* fLosMappedData;

    TH1F* fh_trigger;

    TH1F* fh_land_mapped_barid;
    TH1F* fh_land_mapped_side;
    TH1F* fh_land_mapped_clock;
    TH1F* fh_land_mapped_tac;
    TH1F* fh_land_mapped_qdc;

    TH1F* fh_neuland_mapped_is17;
    TH1F* fh_neuland_mapped_planeid;
    TH1F* fh_neuland_mapped_barid;
    TH1F* fh_neuland_mapped_side;
    TH1F* fh_neuland_mapped_cle;
    TH1F* fh_neuland_mapped_cte;
    TH1F* fh_neuland_mapped_fle;
    TH1F* fh_neuland_mapped_fte;

    TH1F* fh_los_det;
    TH1F* fh_los_ch;
    TH1F* fh_los_tcoarse;
    TH1F* fh_los_tfine;

    void CreateHistos();

    void WriteHistos();

  public:
    ClassDef(R3BNeulandMappedHist, 0)
};

#endif
