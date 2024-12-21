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
// -----                  R3BLandFirstHits header file                 -----
// -----                    Created 20/02/1  by M.Heil                 -----
// -------------------------------------------------------------------------

/**  R3BLandFirstHits.h
 **/

#ifndef R3BLANDFIRSTHITS_H
#define R3BLANDFIRSTHITS_H

#include "TObject.h"
#include "TVector3.h"

class R3BLandFirstHits : public TObject
{

  public:
    /** Default constructor **/
    R3BLandFirstHits();
    R3BLandFirstHits(Double_t x0,
                     Double_t y0,
                     Double_t z0,
                     Double_t t0,
                     Double_t x1,
                     Double_t y1,
                     Double_t z1,
                     Double_t t1,
                     Double_t x2,
                     Double_t y2,
                     Double_t z2,
                     Double_t t2,
                     Double_t x3,
                     Double_t y3,
                     Double_t z3,
                     Double_t t3,
                     Double_t x4,
                     Double_t y4,
                     Double_t z4,
                     Double_t t4,
                     Double_t x5,
                     Double_t y5,
                     Double_t z5,
                     Double_t t5);

    /** Copy constructor **/
    R3BLandFirstHits(const R3BLandFirstHits&);

    /** Destructor **/
    virtual ~R3BLandFirstHits();

    R3BLandFirstHits& operator=(const R3BLandFirstHits&) { return *this; }

    /** Output to screen **/
    virtual void Print(const Option_t* opt) const;

    void SetX0(Double_t pos) { fx0 = pos; }
    Double_t GetX0() { return fx0; }

    void SetY0(Double_t pos) { fy0 = pos; }
    Double_t GetY0() { return fy0; }

    void SetZ0(Double_t pos) { fz0 = pos; }
    Double_t GetZ0() { return fz0; }

    void SetT0(Double_t time) { fT0 = time; }
    Double_t GetT0() { return fT0; }

    void SetX1(Double_t pos) { fx1 = pos; }
    Double_t GetX1() { return fx1; }

    void SetY1(Double_t pos) { fy1 = pos; }
    Double_t GetY1() { return fy1; }

    void SetZ1(Double_t pos) { fz1 = pos; }
    Double_t GetZ1() { return fz1; }

    void SetT1(Double_t time) { fT1 = time; }
    Double_t GetT1() { return fT1; }

    void SetX2(Double_t pos) { fx2 = pos; }
    Double_t GetX2() { return fx2; }

    void SetY2(Double_t pos) { fy2 = pos; }
    Double_t GetY2() { return fy2; }

    void SetZ2(Double_t pos) { fz2 = pos; }
    Double_t GetZ2() { return fz2; }

    void SetT2(Double_t time) { fT2 = time; }
    Double_t GetT2() { return fT2; }

    void SetX3(Double_t pos) { fx3 = pos; }
    Double_t GetX3() { return fx3; }

    void SetY3(Double_t pos) { fy3 = pos; }
    Double_t GetY3() { return fy3; }

    void SetZ3(Double_t pos) { fz3 = pos; }
    Double_t GetZ3() { return fz3; }

    void SetT3(Double_t time) { fT3 = time; }
    Double_t GetT3() { return fT3; }

    void SetX4(Double_t pos) { fx4 = pos; }
    Double_t GetX4() { return fx4; }

    void SetY4(Double_t pos) { fy4 = pos; }
    Double_t GetY4() { return fy4; }

    void SetZ4(Double_t pos) { fz4 = pos; }
    Double_t GetZ4() { return fz4; }

    void SetT4(Double_t time) { fT4 = time; }
    Double_t GetT4() { return fT4; }

    void SetX5(Double_t pos) { fx5 = pos; }
    Double_t GetX5() { return fx5; }

    void SetY5(Double_t pos) { fy5 = pos; }
    Double_t GetY5() { return fy5; }

    void SetZ5(Double_t pos) { fz5 = pos; }
    Double_t GetZ5() { return fz5; }

    void SetT5(Double_t time) { fT5 = time; }
    Double_t GetT5() { return fT5; }

  protected:
    Double32_t fx0;
    Double32_t fy0;
    Double32_t fz0;
    Double32_t fT0;
    Double32_t fx1;
    Double32_t fy1;
    Double32_t fz1;
    Double32_t fT1;
    Double32_t fx2;
    Double32_t fy2;
    Double32_t fz2;
    Double32_t fT2;
    Double32_t fx3;
    Double32_t fy3;
    Double32_t fz3;
    Double32_t fT3;
    Double32_t fx4;
    Double32_t fy4;
    Double32_t fz4;
    Double32_t fT4;
    Double32_t fx5;
    Double32_t fy5;
    Double32_t fz5;
    Double32_t fT5;

    ClassDef(R3BLandFirstHits, 1)
};

#endif
