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
// -----                      R3BTra2pPoint source file                  -----
// -------------------------------------------------------------------------

#include "R3BTra2pDigi.h"

#include <iostream>

using std::cout;
using std::endl;
using std::flush;

// -----   Default constructor   -------------------------------------------
R3BTra2pDigi::R3BTra2pDigi()
{

    Ss03_se_p1 = 0.;
    Ss03_spos_p1 = 0.;
    Ss03_sbw_p1 = 0.;
    Ss03_sarea_p1 = 0.;
    Ss03_seta_p1 = 0.;

    Ss03_ke_p1 = 0.;
    Ss03_kpos_p1 = 0.;
    Ss03_kbw_p1 = 0.;
    Ss03_karea_p1 = 0.;
    Ss03_keta_p1 = 0.;

    Ss06_se_p1 = 0.;
    Ss06_spos_p1 = 0.;
    Ss06_sbw_p1 = 0.;
    Ss06_sarea_p1 = 0.;
    Ss06_seta_p1 = 0.;

    Ss06_ke_p1 = 0.;
    Ss06_kpos_p1 = 0.;
    Ss06_kbw_p1 = 0.;
    Ss06_karea_p1 = 0.;
    Ss06_keta_p1 = 0.;

    Ss03_se_p2 = 0.;
    Ss03_spos_p2 = 0.;
    Ss03_sbw_p2 = 0.;
    Ss03_sarea_p2 = 0.;
    Ss03_seta_p2 = 0.;

    Ss03_ke_p2 = 0.;
    Ss03_kpos_p2 = 0.;
    Ss03_kbw_p2 = 0.;
    Ss03_karea_p2 = 0.;
    Ss03_keta_p2 = 0.;

    Ss06_se_p2 = 0.;
    Ss06_spos_p2 = 0.;
    Ss06_sbw_p2 = 0.;
    Ss06_sarea_p2 = 0.;
    Ss06_seta_p2 = 0.;

    Ss06_ke_p2 = 0.;
    Ss06_kpos_p2 = 0.;
    Ss06_kbw_p2 = 0.;
    Ss06_karea_p2 = 0.;
    Ss06_keta_p2 = 0.;
}

R3BTra2pDigi::R3BTra2pDigi(Double_t ss03_se_p1,
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
                           Double_t ss06_keta_p1,
                           Double_t ss03_se_p2,
                           Double_t ss03_spos_p2,
                           Double_t ss03_sbw_p2,
                           Double_t ss03_sarea_p2,
                           Double_t ss03_seta_p2,
                           Double_t ss03_ke_p2,
                           Double_t ss03_kpos_p2,
                           Double_t ss03_kbw_p2,
                           Double_t ss03_karea_p2,
                           Double_t ss03_keta_p2,
                           Double_t ss06_se_p2,
                           Double_t ss06_spos_p2,
                           Double_t ss06_sbw_p2,
                           Double_t ss06_sarea_p2,
                           Double_t ss06_seta_p2,
                           Double_t ss06_ke_p2,
                           Double_t ss06_kpos_p2,
                           Double_t ss06_kbw_p2,
                           Double_t ss06_karea_p2,
                           Double_t ss06_keta_p2)
{

    Ss03_se_p1 = ss03_se_p1;
    Ss03_spos_p1 = ss03_spos_p1;
    Ss03_sbw_p1 = ss03_sbw_p1;
    Ss03_sarea_p1 = ss03_sarea_p1;
    Ss03_seta_p1 = ss03_seta_p1;

    Ss03_ke_p1 = ss03_ke_p1;
    Ss03_kpos_p1 = ss03_kpos_p1;
    Ss03_kbw_p1 = ss03_kbw_p1;
    Ss03_karea_p1 = ss03_karea_p1;
    Ss03_keta_p1 = ss03_keta_p1;

    Ss06_se_p1 = ss06_se_p1;
    Ss06_spos_p1 = ss06_spos_p1;
    Ss06_sbw_p1 = ss06_sbw_p1;
    Ss06_sarea_p1 = ss06_sarea_p1;
    Ss06_seta_p1 = ss06_seta_p1;

    Ss06_ke_p1 = ss06_ke_p1;
    Ss06_kpos_p1 = ss06_kpos_p1;
    Ss06_kbw_p1 = ss06_kbw_p1;
    Ss06_karea_p1 = ss06_karea_p1;
    Ss06_keta_p1 = ss06_keta_p1;

    Ss03_se_p2 = ss03_se_p2;
    Ss03_spos_p2 = ss03_spos_p2;
    Ss03_sbw_p2 = ss03_sbw_p2;
    Ss03_sarea_p2 = ss03_sarea_p2;
    Ss03_seta_p2 = ss03_seta_p2;

    Ss03_ke_p2 = ss03_ke_p2;
    Ss03_kpos_p2 = ss03_kpos_p2;
    Ss03_kbw_p2 = ss03_kbw_p2;
    Ss03_karea_p2 = ss03_karea_p2;
    Ss03_keta_p2 = ss03_keta_p2;

    Ss06_se_p2 = ss06_se_p2;
    Ss06_spos_p2 = ss06_spos_p2;
    Ss06_sbw_p2 = ss06_sbw_p2;
    Ss06_sarea_p2 = ss06_sarea_p2;
    Ss06_seta_p2 = ss06_seta_p2;

    Ss06_ke_p2 = ss06_ke_p2;
    Ss06_kpos_p2 = ss06_kpos_p2;
    Ss06_kbw_p2 = ss06_kbw_p2;
    Ss06_karea_p2 = ss06_karea_p2;
    Ss06_keta_p2 = ss06_keta_p2;
}

// -----   Destructor   ----------------------------------------------------
R3BTra2pDigi::~R3BTra2pDigi() {}

// -----   Public method Print   -------------------------------------------
void R3BTra2pDigi::Print(const Option_t* opt) const {}
// -------------------------------------------------------------------------

ClassImp(R3BTra2pDigi)
