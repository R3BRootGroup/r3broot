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
// -----                      R3BTraPoint source file                  -----
// -------------------------------------------------------------------------

#include "R3BTraDigi.h"

#include <iostream>

using std::cout;
using std::endl;
using std::flush;

// -----   Default constructor   -------------------------------------------
R3BTraDigi::R3BTraDigi()
{

    Ss03_se_f = 0.;
    Ss03_spos_f = 0.;
    Ss03_sbw_f = 0.;
    Ss03_sarea_f = 0.;
    Ss03_seta_f = 0.;

    Ss03_ke_f = 0.;
    Ss03_kpos_f = 0.;
    Ss03_kbw_f = 0.;
    Ss03_karea_f = 0.;
    Ss03_keta_f = 0.;

    Ss06_se_f = 0.;
    Ss06_spos_f = 0.;
    Ss06_sbw_f = 0.;
    Ss06_sarea_f = 0.;
    Ss06_seta_f = 0.;

    Ss06_ke_f = 0.;
    Ss06_kpos_f = 0.;
    Ss06_kbw_f = 0.;
    Ss06_karea_f = 0.;
    Ss06_keta_f = 0.;

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
}

R3BTraDigi::R3BTraDigi(Double_t ss03_se_f,
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
                       Double_t ss06_keta_p1)
{

    Ss03_se_f = ss03_se_f;
    Ss03_spos_f = ss03_spos_f;
    Ss03_sbw_f = ss03_sbw_f;
    Ss03_sarea_f = ss03_sarea_f;
    Ss03_seta_f = ss03_seta_f;

    Ss03_ke_f = ss03_ke_f;
    Ss03_kpos_f = ss03_kpos_f;
    Ss03_kbw_f = ss03_kbw_f;
    Ss03_karea_f = ss03_karea_f;
    Ss03_keta_f = ss03_keta_f;

    Ss06_se_f = ss06_se_f;
    Ss06_spos_f = ss06_spos_f;
    Ss06_sbw_f = ss06_sbw_f;
    Ss06_sarea_f = ss06_sarea_f;
    Ss06_seta_f = ss06_seta_f;

    Ss06_ke_f = ss06_ke_f;
    Ss06_kpos_f = ss06_kpos_f;
    Ss06_kbw_f = ss06_kbw_f;
    Ss06_karea_f = ss06_karea_f;
    Ss06_keta_f = ss06_keta_f;

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
}

// -----   Destructor   ----------------------------------------------------
R3BTraDigi::~R3BTraDigi() {}

// -----   Public method Print   -------------------------------------------
void R3BTraDigi::Print(const Option_t* opt) const {}
// -------------------------------------------------------------------------

ClassImp(R3BTraDigi)
