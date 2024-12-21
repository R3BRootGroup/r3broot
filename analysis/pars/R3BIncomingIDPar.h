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

#ifndef R3BIncomingIDPar_H
#define R3BIncomingIDPar_H

#include "R3BTcutPar.h"

#include "FairParGenericSet.h"
#include "TCutG.h"
#include "TObject.h"
#include <Rtypes.h>

#include "TArrayF.h"
#include "TObjArray.h"
#include <TObjString.h>

class FairParamList;

class R3BIncomingIDPar : public FairParGenericSet
{
  public:
    /** Standard constructor **/
    R3BIncomingIDPar(const TString& name = "IncomingIDPar",
                     const TString& title = "Incoming particle ID Parameters",
                     const TString& context = "IncomingIDParContext");

    /** Destructor **/
    virtual ~R3BIncomingIDPar();

    /** Reset all parameters **/
    virtual void clear();

    /** Store all parameters using FairRuntimeDB **/
    virtual void putParams(FairParamList* list);

    /** Retrieve all parameters using FairRuntimeDB**/
    Bool_t getParams(FairParamList* list);

    /** Print values of parameters to the standard output **/
    virtual void print();
    void printParams();

    /** Accessor functions **/
    const Float_t GetToFoffset(Int_t det) { return fToFoffset->GetAt(det - 1); }
    const Float_t GetPosS2Left(Int_t det) { return fPosS2Left->GetAt(det - 1); }
    const Float_t GetPosS2Right(Int_t det) { return fPosS2Right->GetAt(det - 1); }
    const Float_t GetBrho0_S2toCC(Int_t det) { return fBrho0_S2toCC->GetAt(det - 1); }
    const Float_t GetDispersionS2(Int_t det) { return fDispersionS2->GetAt(det - 1); }
    const Float_t Getx0_point() { return fx0_point; }
    const Float_t Gety0_point() { return fy0_point; }
    const Float_t Getrot_ang() { return frot_ang; }
    const Float_t Getx0_Aq() { return fx0_Aq; }
    const Float_t Gety0_Aq() { return fy0_Aq; }
    const Float_t Getang_Aq() { return fang_Aq; }
    const Float_t GetTof2InvV_p0(Int_t det) { return fTof2InvV_p0->GetAt(det - 1); }
    const Float_t GetTof2InvV_p1(Int_t det) { return fTof2InvV_p1->GetAt(det - 1); }
    const Float_t GetBeta_min() { return fBeta_min; }
    const Float_t GetBeta_max() { return fBeta_max; }
    TCutG* GetCutS2() { return fTcutparS2->GetCut(); }
    TCutG* GetCutCave() { return fTcutparCave->GetCut(); }

    void SetToFoffset(Float_t value, Int_t det) { fToFoffset->AddAt(value, det - 1); }
    void SetPosS2Left(Float_t value, Int_t det) { fPosS2Left->AddAt(value, det - 1); }
    void SetPosS2Right(Float_t value, Int_t det) { fPosS2Right->AddAt(value, det - 1); }
    void SetBrho0_S2toCC(Float_t value, Int_t det) { fBrho0_S2toCC->AddAt(value, det - 1); }
    void SetDispersionS2(Float_t value, Int_t det) { fDispersionS2->AddAt(value, det - 1); }
    void Setx0_point(Float_t value) { fx0_point = value; }
    void Sety0_point(Float_t value) { fy0_point = value; }
    void Setrot_ang(Float_t value) { frot_ang = value; }
    void SetTof2InvV_p0(Float_t value, Int_t det) { fTof2InvV_p0->AddAt(value, det - 1); }
    void SetTof2InvV_p1(Float_t value, Int_t det) { fTof2InvV_p1->AddAt(value, det - 1); }
    void SetCutS2(TCutG* cut) { fTcutparS2->SetCut(cut); }
    void SetCutCave(TCutG* cut) { fTcutparCave->SetCut(cut); }

  private:
    Int_t fNumDet;
    TArrayF* fToFoffset;
    TArrayF *fPosS2Left, *fPosS2Right;
    TArrayF *fTof2InvV_p0, *fTof2InvV_p1;
    TArrayF* fBrho0_S2toCC;
    TArrayF* fDispersionS2;
    Float_t fx0_point, fy0_point, frot_ang;
    Float_t fx0_Aq, fy0_Aq, fang_Aq;
    Float_t fBeta_max, fBeta_min;

    R3BTcutPar* fTcutparS2;
    R3BTcutPar* fTcutparCave;

    const R3BIncomingIDPar& operator=(const R3BIncomingIDPar&);
    R3BIncomingIDPar(const R3BIncomingIDPar&);

    ClassDef(R3BIncomingIDPar, 1);
};

#endif // R3BIncomingIDPar_H
