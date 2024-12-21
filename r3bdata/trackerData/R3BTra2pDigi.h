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
// -----                      R3BTra2pDigi header file                  -----
// -----                  Created 13/07/11  by D.Bertini/Justyna               -----
// -------------------------------------------------------------------------

/**  R3BTra2pDigi.h
 **/

#ifndef R3BTRA2PDIGI_H
#define R3BTRA2PDIGI_H

#include "TObject.h"
#include "TVector3.h"

class R3BTra2pDigi : public TObject
{

  public:
    /** Default constructor **/
    R3BTra2pDigi();
    R3BTra2pDigi(Double_t ss03_se_p1,
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
                 Double_t ss06_keta_p2);

    /** Copy constructor **/
    R3BTra2pDigi(const R3BTra2pDigi& point) { *this = point; };

    /** Destructor **/
    virtual ~R3BTra2pDigi();

    /** Output to screen **/
    virtual void Print(const Option_t* opt) const;

    void Setss03_se_p1(Double_t en) { Ss03_se_p1 = en; }
    Double_t Getss03_se_p1() { return Ss03_se_p1; }

    void Setss03_spos_p1(Double_t pos) { Ss03_spos_p1 = pos; }
    Double_t Getss03_spos_p1() { return Ss03_spos_p1; }

    void Setss03_sbw_p1(Double_t bw) { Ss03_sbw_p1 = bw; }
    Double_t Getss03_sbw_p1() { return Ss03_sbw_p1; }

    void Setss03_sarea_p1(Double_t area) { Ss03_sarea_p1 = area; }
    Double_t Getss03_sarea_p1() { return Ss03_sarea_p1; }

    void Setss03_seta_p1(Double_t eta) { Ss03_seta_p1 = eta; }
    Double_t Getss03_seta_p1() { return Ss03_seta_p1; }

    void Setss03_ke_p1(Double_t en) { Ss03_ke_p1 = en; }
    Double_t Getss03_ke_p1() { return Ss03_ke_p1; }

    void Setss03_kpos_p1(Double_t pos) { Ss03_kpos_p1 = pos; }
    Double_t Getss03_kpos_p1() { return Ss03_kpos_p1; }

    void Setss03_kbw_p1(Double_t bw) { Ss03_kbw_p1 = bw; }
    Double_t Getss03_kbw_p1() { return Ss03_kbw_p1; }

    void Setss03_karea_p1(Double_t area) { Ss03_karea_p1 = area; }
    Double_t Getss03_karea_p1() { return Ss03_karea_p1; }

    void Setss03_keta_p1(Double_t eta) { Ss03_keta_p1 = eta; }
    Double_t Getss03_keta_p1() { return Ss03_keta_p1; }

    void Setss06_se_p1(Double_t en) { Ss06_se_p1 = en; }
    Double_t Getss06_se_p1() { return Ss06_se_p1; }

    void Setss06_spos_p1(Double_t pos) { Ss06_spos_p1 = pos; }
    Double_t Getss06_spos_p1() { return Ss06_spos_p1; }

    void Setss06_sbw_p1(Double_t bw) { Ss06_sbw_p1 = bw; }
    Double_t Getss06_sbw_p1() { return Ss06_sbw_p1; }

    void Setss06_sarea_p1(Double_t area) { Ss06_sarea_p1 = area; }
    Double_t Getss06_sarea_p1() { return Ss06_sarea_p1; }

    void Setss06_seta_p1(Double_t eta) { Ss06_seta_p1 = eta; }
    Double_t Getss06_seta_p1() { return Ss06_seta_p1; }

    void Setss06_ke_p1(Double_t en) { Ss06_ke_p1 = en; }
    Double_t Getss06_ke_p1() { return Ss06_ke_p1; }

    void Setss06_kpos_p1(Double_t pos) { Ss06_kpos_p1 = pos; }
    Double_t Getss06_kpos_p1() { return Ss06_kpos_p1; }

    void Setss06_kbw_p1(Double_t bw) { Ss06_kbw_p1 = bw; }
    Double_t Getss06_kbw_p1() { return Ss06_kbw_p1; }

    void Setss06_karea_p1(Double_t area) { Ss06_karea_p1 = area; }
    Double_t Getss06_karea_p1() { return Ss06_karea_p1; }

    void Setss06_keta_p1(Double_t eta) { Ss06_keta_p1 = eta; }
    Double_t Getss06_keta_p1() { return Ss06_keta_p1; }

    void Setss03_se_p2(Double_t en) { Ss03_se_p2 = en; }
    Double_t Getss03_se_p2() { return Ss03_se_p2; }

    void Setss03_spos_p2(Double_t pos) { Ss03_spos_p2 = pos; }
    Double_t Getss03_spos_p2() { return Ss03_spos_p2; }

    void Setss03_sbw_p2(Double_t bw) { Ss03_sbw_p2 = bw; }
    Double_t Getss03_sbw_p2() { return Ss03_sbw_p2; }

    void Setss03_sarea_p2(Double_t area) { Ss03_sarea_p2 = area; }
    Double_t Getss03_sarea_p2() { return Ss03_sarea_p2; }

    void Setss03_seta_p2(Double_t eta) { Ss03_seta_p2 = eta; }
    Double_t Getss03_seta_p2() { return Ss03_seta_p2; }

    void Setss03_ke_p2(Double_t en) { Ss03_ke_p2 = en; }
    Double_t Getss03_ke_p2() { return Ss03_ke_p2; }

    void Setss03_kpos_p2(Double_t pos) { Ss03_kpos_p2 = pos; }
    Double_t Getss03_kpos_p2() { return Ss03_kpos_p2; }

    void Setss03_kbw_p2(Double_t bw) { Ss03_kbw_p2 = bw; }
    Double_t Getss03_kbw_p2() { return Ss03_kbw_p2; }

    void Setss03_karea_p2(Double_t area) { Ss03_karea_p2 = area; }
    Double_t Getss03_karea_p2() { return Ss03_karea_p2; }

    void Setss03_keta_p2(Double_t eta) { Ss03_keta_p2 = eta; }
    Double_t Getss03_keta_p2() { return Ss03_keta_p2; }

    void Setss06_se_p2(Double_t en) { Ss06_se_p2 = en; }
    Double_t Getss06_se_p2() { return Ss06_se_p2; }

    void Setss06_spos_p2(Double_t pos) { Ss06_spos_p2 = pos; }
    Double_t Getss06_spos_p2() { return Ss06_spos_p2; }

    void Setss06_sbw_p2(Double_t bw) { Ss06_sbw_p2 = bw; }
    Double_t Getss06_sbw_p2() { return Ss06_sbw_p2; }

    void Setss06_sarea_p2(Double_t area) { Ss06_sarea_p2 = area; }
    Double_t Getss06_sarea_p2() { return Ss06_sarea_p2; }

    void Setss06_seta_p2(Double_t eta) { Ss06_seta_p2 = eta; }
    Double_t Getss06_seta_p2() { return Ss06_seta_p2; }

    void Setss06_ke_p2(Double_t en) { Ss06_ke_p2 = en; }
    Double_t Getss06_ke_p2() { return Ss06_ke_p2; }

    void Setss06_kpos_p2(Double_t pos) { Ss06_kpos_p2 = pos; }
    Double_t Getss06_kpos_p2() { return Ss06_kpos_p2; }

    void Setss06_kbw_p2(Double_t bw) { Ss06_kbw_p2 = bw; }
    Double_t Getss06_kbw_p2() { return Ss06_kbw_p2; }

    void Setss06_karea_p2(Double_t area) { Ss06_karea_p2 = area; }
    Double_t Getss06_karea_p2() { return Ss06_karea_p2; }

    void Setss06_keta_p2(Double_t eta) { Ss06_keta_p2 = eta; }
    Double_t Getss06_keta_p2() { return Ss06_keta_p2; }

  protected:
    Double32_t Ss03_se_p1;
    Double32_t Ss03_spos_p1;
    Double32_t Ss03_sbw_p1;
    Double32_t Ss03_sarea_p1;
    Double32_t Ss03_seta_p1;

    Double32_t Ss03_ke_p1;
    Double32_t Ss03_kpos_p1;
    Double32_t Ss03_kbw_p1;
    Double32_t Ss03_karea_p1;
    Double32_t Ss03_keta_p1;

    Double32_t Ss06_se_p1;
    Double32_t Ss06_spos_p1;
    Double32_t Ss06_sbw_p1;
    Double32_t Ss06_sarea_p1;
    Double32_t Ss06_seta_p1;

    Double32_t Ss06_ke_p1;
    Double32_t Ss06_kpos_p1;
    Double32_t Ss06_kbw_p1;
    Double32_t Ss06_karea_p1;
    Double32_t Ss06_keta_p1;

    Double32_t Ss03_se_p2;
    Double32_t Ss03_spos_p2;
    Double32_t Ss03_sbw_p2;
    Double32_t Ss03_sarea_p2;
    Double32_t Ss03_seta_p2;

    Double32_t Ss03_ke_p2;
    Double32_t Ss03_kpos_p2;
    Double32_t Ss03_kbw_p2;
    Double32_t Ss03_karea_p2;
    Double32_t Ss03_keta_p2;

    Double32_t Ss06_se_p2;
    Double32_t Ss06_spos_p2;
    Double32_t Ss06_sbw_p2;
    Double32_t Ss06_sarea_p2;
    Double32_t Ss06_seta_p2;

    Double32_t Ss06_ke_p2;
    Double32_t Ss06_kpos_p2;
    Double32_t Ss06_kbw_p2;
    Double32_t Ss06_karea_p2;
    Double32_t Ss06_keta_p2;

    ClassDef(R3BTra2pDigi, 1)
};

#endif
