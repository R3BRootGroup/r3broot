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
// -----                       R3BNeutronCalibr2D                          -----
// -----                 Created 27-02-2012 by D.Kresan                    -----
// -----------------------------------------------------------------------------

#ifndef R3BNEUTRONCALIBR2D_H
#define R3BNEUTRONCALIBR2D_H

#include "FairTask.h"

class TClonesArray;
class TH1F;
class TH2F;

class R3BNeutronCalibr2D : public FairTask
{
  public:
    /** Default constructor **/
    R3BNeutronCalibr2D();

    /** Destructor **/
    ~R3BNeutronCalibr2D();

    /** Virtual method Init **/
    virtual InitStatus Init();

    /** Virtual method Exec **/
    virtual void Exec(Option_t* opt);

    virtual void Finish();

    inline void UseBeam(Int_t beamEnergy) { fBeamEnergy = beamEnergy; }

  private:
    Int_t fEventNo;
    Int_t fBeamEnergy;
    TClonesArray* fArrayDigi;
    TClonesArray* fArrayCluster;
    TH1F* fh_etot;
    TH2F* fh_ncl_etot;

    void CreateHistograms();

  public:
    ClassDef(R3BNeutronCalibr2D, 1)
};

#endif
