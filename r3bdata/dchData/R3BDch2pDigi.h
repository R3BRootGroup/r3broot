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
// -----                      R3BDch2pDigi header file                  -----
// -----                  Created 26/03/09  by D.Bertini               -----
// -------------------------------------------------------------------------

/**  R3BDch2pDigi.h
 **/

#ifndef R3BDCH2PDIGI_H
#define R3BDCH2PDIGI_H

#include "TObject.h"
#include "TVector3.h"

class R3BDch2pDigi : public TObject
{

  public:
    /** Default constructor **/
    R3BDch2pDigi();
    R3BDch2pDigi(Int_t pd1mul,
                 Double_t Pdx1_p1,
                 Double_t Pdy1_p1,
                 Double_t Pdx1_p2,
                 Double_t Pdy1_p2,
                 Int_t pd2mul,
                 Double_t Pdx2_p1,
                 Double_t Pdy2_p1,
                 Double_t Pdx2_p2,
                 Double_t Pdy2_p2);

    /** Copy constructor **/
    R3BDch2pDigi(const R3BDch2pDigi&);

    /** Destructor **/
    virtual ~R3BDch2pDigi();

    /** Output to screen **/
    virtual void Print(const Option_t* opt) const;

    void Setpd1mul(Int_t mul) { Pd1mul = mul; }
    Double_t Getpd1mul() { return Pd1mul; }

    void SetPdx1_p1(Double_t pos) { Pd1x_p1 = pos; }
    Double_t GetPdx1_p1() { return Pd1x_p1; }

    void SetPdy1_p1(Double_t pos) { Pd1y_p1 = pos; }
    Double_t GetPdy1_p1() { return Pd1y_p1; }

    void SetPdx1_p2(Double_t pos) { Pd1x_p2 = pos; }
    Double_t GetPdx1_p2() { return Pd1x_p2; }

    void SetPdy1_p2(Double_t pos) { Pd1y_p2 = pos; }
    Double_t GetPdy1_p2() { return Pd1y_p2; }

    void Setpd2mul(Int_t mul) { Pd2mul = mul; }
    Double_t Getpd2mul() { return Pd2mul; }

    void SetPdx2_p1(Double_t pos) { Pd2x_p1 = pos; }
    Double_t GetPdx2_p1() { return Pd2x_p1; }

    void SetPdy2_p1(Double_t pos) { Pd2y_p1 = pos; }
    Double_t GetPdy2_p1() { return Pd2y_p1; }

    void SetPdx2_p2(Double_t pos) { Pd2x_p2 = pos; }
    Double_t GetPdx2_p2() { return Pd2x_p2; }

    void SetPdy2_p2(Double_t pos) { Pd2y_p2 = pos; }
    Double_t GetPdy2_p2() { return Pd2y_p2; }

    R3BDch2pDigi& operator=(const R3BDch2pDigi&) { return *this; }

  protected:
    Int_t Pd1mul;
    Double32_t Pd1x_p1;
    Double32_t Pd1y_p1;
    Double32_t Pd1x_p2;
    Double32_t Pd1y_p2;

    Int_t Pd2mul;
    Double32_t Pd2x_p1;
    Double32_t Pd2y_p1;
    Double32_t Pd2x_p2;
    Double32_t Pd2y_p2;

    ClassDef(R3BDch2pDigi, 1)
};

#endif
