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
// -----                      R3BTofPoint source file                  -----
// -------------------------------------------------------------------------

#include "R3BTof2pDigi.h"

#include <iostream>

using std::cout;
using std::endl;
using std::flush;

// -----   Default constructor   -------------------------------------------
R3BTof2pDigi::R3BTof2pDigi()
{

    Tfmul = 0;
    Tfx_p1 = 0.;
    Tfy_p1 = 0.;
    Tft_p1 = 0.;
    Tfpath_p1 = 0.;
    Tfwpx_p1 = 0.;
    Tfwpy_p1 = 0.;
    Tfwpz_p1 = 0.;
    Tfe_p1 = 0.;

    Tfx_p2 = 0.;
    Tfy_p2 = 0.;
    Tft_p2 = 0.;
    Tfpath_p2 = 0.;
    Tfwpx_p2 = 0.;
    Tfwpy_p2 = 0.;
    Tfwpz_p2 = 0.;
    Tfe_p2 = 0.;
}

R3BTof2pDigi::R3BTof2pDigi(
    Int_t tfmul,
    Double_t tfwx_p1,
    Double_t tfwy_p1,
    Double_t tfwt_p1,
    Double_t tfwx_p2,
    Double_t tfwy_p2,
    Double_t tfwt_p2,
    Double_t tfwpath_p1,
    Double_t tfwpath_p2,
    Double_t tfwpx_p1,
    // Double_t tfwpy_p1,Double_t tfwpz_p1,Double_t tfwpx_p2,Double_t tfwpy_p2,Double_t tfwpz_p2) {
    Double_t tfwpy_p1,
    Double_t tfwpz_p1,
    Double_t tfwpx_p2,
    Double_t tfwpy_p2,
    Double_t tfwpz_p2,
    Double_t tfwe_p1,
    Double_t tfwe_p2)
{

    Tfmul = tfmul;
    Tfx_p1 = tfwx_p1;
    Tfy_p1 = tfwy_p1;
    Tft_p1 = tfwt_p1;
    Tfpath_p1 = tfwpath_p1;
    Tfwpx_p1 = tfwpx_p1;
    Tfwpy_p1 = tfwpy_p1;
    Tfwpz_p1 = tfwpz_p1;
    Tfe_p1 = tfwe_p1;

    Tfx_p2 = tfwx_p2;
    Tfy_p2 = tfwy_p2;
    Tft_p2 = tfwt_p2;
    Tfpath_p2 = tfwpath_p2;
    Tfwpx_p2 = tfwpx_p2;
    Tfwpy_p2 = tfwpy_p2;
    Tfwpz_p2 = tfwpz_p2;
    Tfe_p2 = tfwe_p2;
}

// -----   Destructor   ----------------------------------------------------
R3BTof2pDigi::~R3BTof2pDigi() {}

// -----   Public method Print   -------------------------------------------
void R3BTof2pDigi::Print(const Option_t* opt) const {}
// -------------------------------------------------------------------------

ClassImp(R3BTof2pDigi)
