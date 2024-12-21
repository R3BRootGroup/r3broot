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

#ifndef R3BTRA2PDIGITIZER_H
#define R3BTRA2PDIGITIZER_H 1

#include "FairTask.h"
#include "R3BTra2pDigi.h"
#include "R3BTra2pDigiPar.h"
#include <map>
#include <string>

class TClonesArray;
class TObjectArray;
class TH1F;
class TH2F;

class R3BTra2pDigitizer : public FairTask
{

  public:
    /** Default constructor **/
    R3BTra2pDigitizer();

    /** Destructor **/
    ~R3BTra2pDigitizer();

    /** Virtual method Init **/
    virtual InitStatus Init();

    /** Virtual method Exec **/
    virtual void Exec(Option_t* opt);

    virtual void Finish();
    virtual void Reset();

    R3BTra2pDigi* AddHit(Double_t ss03_se_p1,
                         Double_t ss03_spos_p1,
                         Double_t ss03_sbw_p1,
                         Double_t ss03_sarea_p1,
                         Double_t ss03_seta_p1,
                         Double_t ss03_ke_p1,
                         Double_t ss03_kpos_p1,
                         Double_t ss03_kbw_p1,
                         Double_t ss03_karea_p1,
                         Double_t ss03_keta_p1,
                         Double_t ss06_se_p1,
                         Double_t ss06_spos_p1,
                         Double_t ss06_sbw_p1,
                         Double_t ss06_sarea_p1,
                         Double_t ss06_seta_p1,
                         Double_t ss06_ke_p1,
                         Double_t ss06_kpos_p1,
                         Double_t ss06_kbw_p1,
                         Double_t ss06_karea_p1,
                         Double_t ss06_keta_p1,
                         Double_t ss03_se_p2,
                         Double_t ss03_spos_p2,
                         Double_t ss03_sbw_p2,
                         Double_t ss03_sarea_p2,
                         Double_t ss03_seta_p2,
                         Double_t ss03_ke_p2,
                         Double_t ss03_kpos_p2,
                         Double_t ss03_kbw_p2,
                         Double_t ss03_karea_p2,
                         Double_t ss03_keta_p2,
                         Double_t ss06_se_p2,
                         Double_t ss06_spos_p2,
                         Double_t ss06_sbw_p2,
                         Double_t ss06_sarea_p2,
                         Double_t ss06_seta_p2,
                         Double_t ss06_ke_p2,
                         Double_t ss06_kpos_p2,
                         Double_t ss06_kbw_p2,
                         Double_t ss06_karea_p2,
                         Double_t ss06_keta_p2);

  protected:
    TClonesArray* fTra2pPoints;
    TClonesArray* fTra2pMCTrack;
    TClonesArray* fTra2pDigi;
    TClonesArray* fMCTrack;

    // Parameter class
    R3BTra2pDigiPar* fTra2pDigiPar;

    //- Control Hitograms

    Int_t eventNoTra;

  private:
    virtual void SetParContainers();

    ClassDef(R3BTra2pDigitizer, 1);
};

#endif
