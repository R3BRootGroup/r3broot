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

// -------------------------------------------------------------------------
// -----                      R3BTargetDigi header file                  -----
// -----                  Created 30/06/11  by D.Bertini/Justyna               -----
// -------------------------------------------------------------------------

/**  R3BTargetDigi.h
 **/

#ifndef R3BTARGETDIGI_H
#define R3BTARGETDIGI_H

#include "TObject.h"
#include "TVector3.h"

class R3BTargetDigi : public TObject
{

  public:
    /** Default constructor **/
    R3BTargetDigi();
    R3BTargetDigi(Double_t ss03_smul,
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

    /** Copy constructor **/
    R3BTargetDigi(const R3BTargetDigi& point) { *this = point; };

    /** Destructor **/
    virtual ~R3BTargetDigi();

    /** Output to screen **/
    virtual void Print(const Option_t* opt) const;

    void Setss03_smul(Double_t mul) { Ss03_smul = mul; }
    Double_t Getss03_smul() { return Ss03_smul; }

    void Setss03_kmul(Double_t mul) { Ss03_kmul = mul; }
    Double_t Getss03_kmul() { return Ss03_kmul; }

    void Setss06_smul(Double_t mul) { Ss06_smul = mul; }
    Double_t Getss06_smul() { return Ss06_smul; }

    void Setss06_kmul(Double_t mul) { Ss06_kmul = mul; }
    Double_t Getss06_kmul() { return Ss06_kmul; }

    void Setx0(Double_t pos) { X0 = pos; }
    Double_t Getx0() { return X0; }

    void Sety0(Double_t pos) { Y0 = pos; }
    Double_t Gety0() { return Y0; }

    void Sett0(Double_t time) { T0 = time; }
    Double_t Gett0() { return T0; }

    void Setestar(Double_t energy) { Estar = energy; }
    Double_t Getestar() { return Estar; }

    void SetPxf(Double_t mom) { pxf = mom; }
    Double_t GetPxf() { return pxf; }
    void SetPyf(Double_t mom) { pyf = mom; }
    Double_t GetPyf() { return pyf; }
    void SetPzf(Double_t mom) { pzf = mom; }
    Double_t GetPzf() { return pzf; }
    void SetPf_tot(Double_t mom) { pf_tot = mom; }
    Double_t GetPf_tot() { return pf_tot; }
    void Setf_beta(Double_t vel) { f_Beta = vel; }
    Double_t Getf_beta() { return f_Beta; }

    void SetPxp1(Double_t mom) { pxp1 = mom; }
    Double_t GetPxp1() { return pxp1; }
    void SetPyp1(Double_t mom) { pyp1 = mom; }
    Double_t GetPyp1() { return pyp1; }
    void SetPzp1(Double_t mom) { pzp1 = mom; }
    Double_t GetPzp1() { return pzp1; }
    void SetPp1_tot(Double_t mom) { pp1_tot = mom; }
    Double_t GetPp1_tot() { return pp1_tot; }
    void Setp1_beta(Double_t vel) { p1_Beta = vel; }
    Double_t Getp1_beta() { return p1_Beta; }

  protected:
    Double32_t Ss03_smul;
    Double32_t Ss03_kmul;
    Double32_t Ss06_smul;
    Double32_t Ss06_kmul;

    Double32_t X0;
    Double32_t Y0;
    Double32_t T0;

    Double32_t Estar;

    Double32_t pxf;
    Double32_t pyf;
    Double32_t pzf;
    Double32_t pf_tot;
    Double32_t f_Beta;

    Double32_t pxp1;
    Double32_t pyp1;
    Double32_t pzp1;
    Double32_t pp1_tot;
    Double32_t p1_Beta;

    ClassDef(R3BTargetDigi, 1)
};

#endif
