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

#ifndef R3BTRADIGITIZER_H
#define R3BTRADDIGITISER_H 1

#include "FairTask.h"
#include "R3BTraDigi.h"
#include "R3BTraDigiPar.h"
#include <map>
#include <string>

class TClonesArray;
class TObjectArray;
class TH1F;
class TH2F;

class R3BTraDigitizer : public FairTask
{

  public:
    /** Default constructor **/
    R3BTraDigitizer();

    /** Destructor **/
    ~R3BTraDigitizer();

    /** Virtual method Init **/
    virtual InitStatus Init();

    /** Virtual method Exec **/
    virtual void Exec(Option_t* opt);

    virtual void Finish();
    virtual void Reset();

    R3BTraDigi* AddHit(Double_t ss03_se_f,
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
                       Double_t ss06_keta_f,
                       Double_t ss03_se_p1,
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
                       Double_t ss06_keta_p1);

  protected:
    TClonesArray* fTraPoints;
    TClonesArray* fTraMCTrack;
    TClonesArray* fTraDigi;
    TClonesArray* fMCTrack;

    // Parameter class
    R3BTraDigiPar* fTraDigiPar;

    //- Control Hitograms

    Int_t eventNoTra;

  private:
    virtual void SetParContainers();

    ClassDef(R3BTraDigitizer, 1);
};

#endif
