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

#ifndef R3BTRAFRADIGITIZER_H
#define R3BTRA2FRADIGITISER_H 1

#include "FairTask.h"
#include "R3BTraFraDigi.h"
#include "R3BTraFraDigiPar.h"
#include <map>
#include <string>

class TClonesArray;
class TObjectArray;
class TH1F;
class TH2F;

class R3BTraFraDigitizer : public FairTask
{

  public:
    /** Default constructor **/
    R3BTraFraDigitizer();

    /** Destructor **/
    ~R3BTraFraDigitizer();

    /** Virtual method Init **/
    virtual InitStatus Init();

    /** Virtual method Exec **/
    virtual void Exec(Option_t* opt);

    virtual void Finish();
    virtual void Reset();

    R3BTraFraDigi* AddHit(Double_t ss03_se_f,
                          Double_t ss03_spos_f,
                          Double_t ss03_sbw_f,
                          Double_t ss03_sarea_f,
                          Double_t ss03_seta_f,
                          Double_t ss03_ke_f,
                          Double_t ss03_kpos_f,
                          Double_t ss03_kbw_f,
                          Double_t ss03_karea_f,
                          Double_t ss03_keta_f,
                          Double_t ss06_se_f,
                          Double_t ss06_spos_f,
                          Double_t ss06_sbw_f,
                          Double_t ss06_sarea_f,
                          Double_t ss06_seta_f,
                          Double_t ss06_ke_f,
                          Double_t ss06_kpos_f,
                          Double_t ss06_kbw_f,
                          Double_t ss06_karea_f,
                          Double_t ss06_keta_f);

  protected:
    TClonesArray* fTraFraPoints;
    TClonesArray* fTraFraMCTrack;
    TClonesArray* fTraFraDigi;
    TClonesArray* fMCTrack;

    // Parameter class
    R3BTraFraDigiPar* fTraFraDigiPar;

    //- Control Hitograms

    Int_t eventNoTra;

  private:
    virtual void SetParContainers();

    ClassDef(R3BTraFraDigitizer, 1);
};

#endif
