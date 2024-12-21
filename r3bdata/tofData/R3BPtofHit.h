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

#ifndef R3BPTOFHIT_H
#define R3BPTOFHIT_H

#include "TObject.h"

class R3BPtofHit : public TObject
{

  public:
    /** Default constructor **/
    R3BPtofHit();
    R3BPtofHit(Int_t paddle, Double_t charge, Double_t tdc, Double_t xx, Double_t yy, Double_t zz);

    /** Copy constructor **/
    R3BPtofHit(const R3BPtofHit&);

    /** Destructor **/
    virtual ~R3BPtofHit();

    R3BPtofHit& operator=(const R3BPtofHit&) { return *this; }

    /** Output to screen **/
    virtual void Print(const Option_t* opt) const;

    void SetCharge(Double_t charge) { fCharge = charge; }
    Double_t GetCharge() const { return fCharge; }

    void SetTdc(Double_t time) { fTdc = time; }
    Double_t GetTdc() const { return fTdc; }

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
    Double32_t fCharge;
    Double32_t fTdc;
    Double32_t fxx;
    Double32_t fyy;
    Double32_t fzz;

    ClassDef(R3BPtofHit, 1)
};

#endif
