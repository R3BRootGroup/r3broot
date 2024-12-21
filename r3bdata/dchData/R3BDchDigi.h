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
// -----                      R3BDchDigi header file                  -----
// -----                  Created 26/03/09  by D.Bertini               -----
// -------------------------------------------------------------------------

/**  R3BDchDigi.h
 **/

#ifndef R3BDCHDIGI_H
#define R3BDCHDIGI_H

#include "TObject.h"
#include "TVector3.h"

class R3BDchDigi : public TObject
{

  public:
    /** Default constructor **/
    R3BDchDigi();
    R3BDchDigi(Int_t pd1mul, Double_t Pdx1_p1, Double_t Pdy1_p1, Int_t pd2mul, Double_t Pdx2_p1, Double_t Pdy2_p1);

    /** Copy constructor **/
    R3BDchDigi(const R3BDchDigi& point);

    /** Destructor **/
    virtual ~R3BDchDigi();

    /** Output to screen **/
    virtual void Print(const Option_t* opt) const;

    void Setpd1mul(Int_t mul) { Pd1mul = mul; }
    Double_t Getpd1mul() { return Pd1mul; }

    void SetPdx1_p1(Double_t pos) { Pd1x_p1 = pos; }
    Double_t GetPdx1_p1() { return Pd1x_p1; }

    void SetPdy1_p1(Double_t pos) { Pd1y_p1 = pos; }
    Double_t GetPdy1_p1() { return Pd1y_p1; }

    void Setpd2mul(Int_t mul) { Pd2mul = mul; }
    Double_t Getpd2mul() { return Pd2mul; }

    void SetPdx2_p1(Double_t pos) { Pd2x_p1 = pos; }
    Double_t GetPdx2_p1() { return Pd2x_p1; }

    void SetPdy2_p1(Double_t pos) { Pd2y_p1 = pos; }
    Double_t GetPdy2_p1() { return Pd2y_p1; }

    R3BDchDigi& operator=(const R3BDchDigi&) { return *this; }

  protected:
    Int_t Pd1mul;
    Double32_t Pd1x_p1;
    Double32_t Pd1y_p1;

    Int_t Pd2mul;
    Double32_t Pd2x_p1;
    Double32_t Pd2y_p1;

    ClassDef(R3BDchDigi, 1)
};

#endif
