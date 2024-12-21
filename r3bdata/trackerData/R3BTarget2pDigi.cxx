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
// -----                      R3BTarget2pPoint source file                  -----
// -------------------------------------------------------------------------

#include "R3BTarget2pDigi.h"

#include <iostream>

using std::cout;
using std::endl;
using std::flush;

// -----   Default constructor   -------------------------------------------
R3BTarget2pDigi::R3BTarget2pDigi()
{

    Ss03_smul = 0;
    Ss03_kmul = 0;
    Ss06_smul = 0;
    Ss06_kmul = 0;

    X0 = 0.;
    Y0 = 0.;
    T0 = 0.;

    Estar = 0;

    pxf = 0;
    pyf = 0;
    pzf = 0;
    pf_tot = 0;
    f_Beta = 0;

    pxp1 = 0;
    pyp1 = 0;
    pzp1 = 0;
    pp1_tot = 0;
    p1_Beta = 0;

    pxp2 = 0;
    pyp2 = 0;
    pzp2 = 0;
    pp2_tot = 0;
    p2_Beta = 0;

    Alpha = 0;

    Ppmul = 0;
}

R3BTarget2pDigi::R3BTarget2pDigi(Double_t ss03_smul,
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
                                 Double_t p1_beta,
                                 Double_t Pxp2,
                                 Double_t Pyp2,
                                 Double_t Pzp2,
                                 // Double_t Pp2_tot,Double_t p2_beta) {
                                 Double_t Pp2_tot,
                                 Double_t p2_beta,
                                 Double_t alpha,
                                 Int_t ppmul)
{

    Ss03_smul = ss03_smul;
    Ss03_kmul = ss03_kmul;
    Ss06_smul = ss06_smul;
    Ss06_kmul = ss06_kmul;

    X0 = x0;
    Y0 = y0;
    T0 = t0;

    Estar = estar;

    pxf = Pxf;
    pyf = Pyf;
    pzf = Pzf;
    pf_tot = Pf_tot;
    f_Beta = f_beta;

    pxp1 = Pxp1;
    pyp1 = Pyp1;
    pzp1 = Pzp1;
    pp1_tot = Pp1_tot;
    p1_Beta = p1_beta;

    pxp2 = Pxp2;
    pyp2 = Pyp2;
    pzp2 = Pzp2;
    pp2_tot = Pp2_tot;
    p2_Beta = p2_beta;

    Alpha = alpha;

    Ppmul = ppmul;
}

// -----   Destructor   ----------------------------------------------------
R3BTarget2pDigi::~R3BTarget2pDigi() {}

// -----   Public method Print   -------------------------------------------
void R3BTarget2pDigi::Print(const Option_t* opt) const {}
// -------------------------------------------------------------------------

ClassImp(R3BTarget2pDigi)
