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
// -----                      R3BTofDigi header file                  -----
// -----                  Created 29/06/11  by D.Bertini/Justyna               -----
// -------------------------------------------------------------------------

/**  R3BTofDigi.h
 **/

#ifndef R3BTOFDIGI_H
#define R3BTOFDIGI_H

#include "TObject.h"
#include "TVector3.h"

class R3BTofDigi : public TObject
{

  public:
    /** Default constructor **/
    R3BTofDigi();
    R3BTofDigi(Int_t tfmul,
               Double_t tfwx_p1,
               Double_t tfwy_p1,
               Double_t tfwt_p1,
               Double_t tfwpath_p1,
               Double_t tfwpx_p1,
               // Double_t tfwpy_p1,Double_t tfwpz_p1);
               Double_t tfwpy_p1,
               Double_t tfwpz_p1,
               Double_t tfwe_p1);

    /** Copy constructor **/
    R3BTofDigi(const R3BTofDigi& point) { *this = point; };

    /** Destructor **/
    virtual ~R3BTofDigi();

    /** Output to screen **/
    virtual void Print(const Option_t* opt) const;

    void Settfmul(Int_t mul) { Tfmul = mul; }
    Double_t Gettfmul() { return Tfmul; }

    void Settfwx_p1(Double_t pos) { Tfx_p1 = pos; }
    Double_t Gettfwx_p1() { return Tfx_p1; }

    void Settfwy_p1(Double_t pos) { Tfy_p1 = pos; }
    Double_t Gettfwy_p1() { return Tfy_p1; }

    void Settfwt_p1(Double_t time) { Tft_p1 = time; }
    Double_t Gettfwt_p1() { return Tft_p1; }

    void Settfwpath_p1(Double_t length) { Tfpath_p1 = length; }
    Double_t Gettfwpath_p1() { return Tfpath_p1; }

    void Settfwpx_p1(Int_t mom) { Tfwpx_p1 = mom; }
    Double_t Gettfwpx_p1() { return Tfwpx_p1; }

    void Settfwpy_p1(Int_t mom) { Tfwpy_p1 = mom; }
    Double_t Gettfwpy_p1() { return Tfwpy_p1; }

    void Settfwpz_p1(Int_t mom) { Tfwpz_p1 = mom; }
    Double_t Gettfwpz_p1() { return Tfwpz_p1; }

    void Settfwe_p1(Double_t en) { Tfe_p1 = en; }
    Double_t Gettfwe_p1() { return Tfe_p1; }

  protected:
    Int_t Tfmul;
    Double32_t Tfx_p1;
    Double32_t Tfy_p1;
    Double32_t Tft_p1;
    Double32_t Tfpath_p1;
    Double32_t Tfwpx_p1;
    Double32_t Tfwpy_p1;
    Double32_t Tfwpz_p1;
    Double32_t Tfe_p1;

    ClassDef(R3BTofDigi, 1)
};

#endif
