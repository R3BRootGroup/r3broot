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
// -----                      R3BTofdPoint source file                  -----
// -------------------------------------------------------------------------

#include "R3BTofdDigi.h"

#include <iostream>

using std::cout;
using std::endl;
using std::flush;

// -----   Default constructor   -------------------------------------------
R3BTofdDigi::R3BTofdDigi()
{

    Ntmul = 0;
    Ntx = 0.;
    Nty = 0.;
    Ntt = 0.;
    Ntpath = 0.;
    Ntfpx = 0.;
    Ntfpy = 0.;
    Ntfpz = 0.;
    Nte = 0.;
}

R3BTofdDigi::R3BTofdDigi(Int_t ntmul,
                         Double_t ntfx,
                         Double_t ntfy,
                         Double_t ntft,
                         Double_t ntfpath,
                         Double_t ntfpx,
                         // Double_t ntfpy,Double_t ntfpz) {
                         Double_t ntfpy,
                         Double_t ntfpz,
                         Double_t ntfe)
{

    Ntmul = ntmul;
    Ntx = ntfx;
    Nty = ntfy;
    Ntt = ntft;
    Ntpath = ntfpath;
    Ntfpx = ntfpx;
    Ntfpy = ntfpy;
    Ntfpz = ntfpz;
    Nte = ntfe;
}

// -----   Destructor   ----------------------------------------------------
R3BTofdDigi::~R3BTofdDigi() {}

// -----   Public method Print   -------------------------------------------
void R3BTofdDigi::Print(const Option_t* opt) const {}
// -------------------------------------------------------------------------

ClassImp(R3BTofdDigi)
