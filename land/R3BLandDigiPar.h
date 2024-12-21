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

//
// C++ Interface: PndTutPar
//
#ifndef R3BLANDDIGIPAR_H
#define R3BLANDDIGIPAR_H

#include <TObjString.h>
#include <TVector3.h>

#include "FairParGenericSet.h"
#include "FairParamList.h"

class R3BLandDigiPar : public FairParGenericSet
{
  public:
    R3BLandDigiPar(const char* name = "R3BLandDigiParTest",
                   const char* title = "Tutorial  parameter",
                   const char* context = "TestDefaultContext");
    ~R3BLandDigiPar(void){};
    void clear(void){};
    void putParams(FairParamList* list);
    Bool_t getParams(FairParamList* list);

    void Print(Option_t* option = "") const;
    /** Accessor functions **/
    const Int_t GetMaxPaddle() { return nMaxPaddle; }
    const Int_t GetMaxPlane() { return nMaxPlane; }
    const Double_t GetPaddleLength() { return paddleLength; }
    const Double_t GetPaddleHeight() { return paddleHeight; }
    const Double_t GetPaddleDepth() { return paddleDepth; }
    const Double_t GetPaddleSpacing() { return paddleSpacing; }
    const Double_t GetPaddleWrapping() { return paddleWrapping; }
    const TString GetGeometryFileName() { return fGeometryFileName; }
    void SetMaxPaddle(Int_t paddle) { nMaxPaddle = paddle; }
    void SetMaxPlane(Int_t pla) { nMaxPlane = pla; }
    void SetPaddleLength(Double_t length) { paddleLength = length; }
    void SetPaddleHeight(Double_t height) { paddleHeight = height; }
    void SetPaddleDepth(Double_t depth) { paddleDepth = depth; }
    void SetPaddleSpacing(Double_t space) { paddleSpacing = space; }
    void SetPaddleWrapping(Double_t wrap) { paddleWrapping = wrap; }
    void SetGeometryFileName(TString name) { fGeometryFileName = name; }

  private:
    // Digi. Parameters
    Int_t nMaxPaddle;
    Int_t nMaxPlane;
    Double_t paddleLength;
    Double_t paddleHeight;
    Double_t paddleDepth;
    Double_t paddleSpacing;
    Double_t paddleWrapping;
    TString fGeometryFileName;

    ClassDef(R3BLandDigiPar, 2); //
};

#endif /* !R3BLANDDIGIPAR_H*/
