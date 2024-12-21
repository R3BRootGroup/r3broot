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

#ifndef R3BTARGETDIGITIZER_H
#define R3BTARGETDDIGITISER_H 1

#include "FairTask.h"
#include "R3BTargetDigi.h"
#include "R3BTargetDigiPar.h"
#include <map>
#include <string>

class TClonesArray;
class TObjectArray;
class TH1F;
class TH2F;

class R3BTargetDigitizer : public FairTask
{

  public:
    /** Default constructor **/
    R3BTargetDigitizer();

    /** Destructor **/
    ~R3BTargetDigitizer();

    /** Virtual method Init **/
    virtual InitStatus Init();

    /** Virtual method Exec **/
    virtual void Exec(Option_t* opt);

    virtual void Finish();
    virtual void Reset();

    R3BTargetDigi* AddHit(Double_t ss03_smul,
                          Double_t ss03_kmul,
                          Double_t ss06_smul,
                          Double_t ss06_kmul,
                          Double_t x0,
                          Double_t y0,
                          Double_t t0,
                          Double_t estar,
                          Double_t Pxf,
                          Double_t Pyf,
                          Double_t Pzf,
                          Double_t Pf_tot,
                          Double_t f_beta,
                          Double_t Pxp1,
                          Double_t Pyp1,
                          Double_t Pzp1,
                          Double_t Pp1_tot,
                          Double_t p1_beta);

  protected:
    TClonesArray* fTargetPoints;
    TClonesArray* fTargetMCTrack;
    TClonesArray* fTargetDigi;
    TClonesArray* fMCTrack;

    // Parameter class
    R3BTargetDigiPar* fTargetDigiPar;

    //- Control Hitograms

    TH1F* ExEnIn_his;

    Int_t eventNoTra;

  private:
    virtual void SetParContainers();

    ClassDef(R3BTargetDigitizer, 1);
};

#endif
