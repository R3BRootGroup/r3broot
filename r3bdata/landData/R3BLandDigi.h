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
// -----                      R3BLandDigi header file                  -----
// -----                  Created 26/03/09  by D.Bertini               -----
// -------------------------------------------------------------------------

/**  R3BLandDigi.h
 **/

#ifndef R3BLANDDIGI_H
#define R3BLANDDIGI_H

#include "TObject.h"
#include "TVector3.h"

class R3BLandDigi : public TObject
{

  public:
    /** Default constructor **/
    R3BLandDigi();
    R3BLandDigi(Int_t paddle,
                Double_t tdcL,
                Double_t tdcR,
                Double_t tdc,
                Double_t qdcL,
                Double_t qdcR,
                Double_t qdc,
                Double_t xx,
                Double_t yy,
                Double_t zz);

    R3BLandDigi(Int_t paddle, Double_t tdcR, Double_t tdcL, Double_t qdcR, Double_t qdcL);

    /** Copy constructor **/
    R3BLandDigi(const R3BLandDigi&);

    /** Destructor **/
    virtual ~R3BLandDigi();

    R3BLandDigi& operator=(const R3BLandDigi&) { return *this; }

    /** Output to screen **/
    virtual void Print(const Option_t* opt) const;

    void SetTdcR(Double_t time) { fTdcR = time; }
    Double_t GetTdcR() { return fTdcR; }

    void SetTdcL(Double_t time) { fTdcL = time; }
    Double_t GetTdcL() { return fTdcL; }

    void SetTdc(Double_t time) { fTdc = time; }
    Double_t GetTdc() const { return fTdc; }

    void SetQdcR(Double_t charge) { fQdcR = charge; }
    Double_t GetQdcR() { return fQdcR; }

    void SetQdcL(Double_t charge) { fQdcL = charge; }
    Double_t GetQdcL() { return fQdcL; }

    void SetQdc(Double_t charge) { fQdc = charge; }
    Double_t GetQdc() { return fQdc; }

    void SetXX(Double_t pos) { fxx = pos; }
    Double_t GetXX() { return fxx; }

    void SetYY(Double_t pos) { fyy = pos; }
    Double_t GetYY() { return fyy; }

    void SetZZ(Double_t pos) { fzz = pos; }
    Double_t GetZZ() { return fzz; }

    void SetPaddleNr(Int_t paddleNb) { fPaddleNb = paddleNb; }
    Int_t GetPaddleNr() { return fPaddleNb; }

  protected:
    Int_t fPaddleNb;
    Double32_t fTdcL;
    Double32_t fTdcR;
    Double32_t fTdc;
    Double32_t fQdcL;
    Double32_t fQdcR;
    Double32_t fQdc;
    Double32_t fxx;
    Double32_t fyy;
    Double32_t fzz;

    ClassDef(R3BLandDigi, 1)
};

#endif
