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
// -----                      R3BGfiDigi header file                  -----
// -----                  Created 30/06/11  by D.Bertini/Justyna               -----
// -------------------------------------------------------------------------

/**  R3BGfiDigi.h
 **/

#ifndef R3BGFIDIGI_H
#define R3BGFIDIGI_H

#include "TObject.h"
#include "TVector3.h"

class R3BGfiDigi : public TObject
{

  public:
    /** Default constructor **/
    R3BGfiDigi();
    R3BGfiDigi(Int_t gf1mul, Double_t gf1_01x, Int_t gf2mul, Double_t gf2_01x);

    /** Copy constructor **/
    R3BGfiDigi(const R3BGfiDigi&);

    /** Destructor **/
    virtual ~R3BGfiDigi();

    /** Output to screen **/
    virtual void Print(const Option_t* opt) const;

    void Setgf1mul(Int_t mul) { Gf1mul = mul; }
    Double_t Getgf1mul() { return Gf1mul; }

    void Setgf1_01x(Double_t pos) { Gf1_01x = pos; }
    Double_t Getgf1_01x() { return Gf1_01x; }

    void Setgf2mul(Int_t mul) { Gf2mul = mul; }
    Double_t Getgf2mul() { return Gf2mul; }

    void Setgf2_01x(Double_t pos) { Gf2_01x = pos; }
    Double_t Getgf2_01x() { return Gf2_01x; }

  protected:
    Int_t Gf1mul;
    Double32_t Gf1_01x;
    Int_t Gf2mul;
    Double32_t Gf2_01x;

    ClassDef(R3BGfiDigi, 1)
};

#endif
