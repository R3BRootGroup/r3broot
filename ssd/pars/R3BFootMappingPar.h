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

// ----------------------------------------------------------------
// -----        R3BFootMappingPar source file                 -----
// -----    Created 05/11/21 by J.L. Rodriguez-Sanchez        -----
// ----------------------------------------------------------------

#pragma once

#include <FairParGenericSet.h>

#include <TArrayF.h>
#include <TArrayI.h>
#include <TObject.h>

class FairParamList;

class R3BFootMappingPar : public FairParGenericSet
{
  public:
    /** Standard constructor **/
    R3BFootMappingPar(const TString& name = "footMappingPar",
                      const TString& title = "FOOT Mapping Parameters",
                      const TString& context = "footMappingContext");

    /** Destructor **/
    virtual ~R3BFootMappingPar();

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
    const Int_t GetGeometry() { return fGeo; }
    const Int_t GetNumDets() { return fNumDet; }
    const Float_t GetDist2target(Int_t det) { return fDistance2target->GetAt(det - 1); }
    const Float_t GetAngleTheta(Int_t det) { return fAngleTheta->GetAt(det - 1); }
    const Float_t GetAnglePhi(Int_t det) { return fAnglePhi->GetAt(det - 1); }
    const Float_t GetOffsetX(Int_t det) { return fOffsetX->GetAt(det - 1); }
    const Float_t GetOffsetY(Int_t det) { return fOffsetY->GetAt(det - 1); }
    const Int_t GetPolPar() { return fPolPar; }
    const Float_t GetEnevsPosCorr(Int_t order) { return fEnevsPosCorr->GetAt(order - 1); }

    void SetGeometry(Int_t geo) { fGeo = geo; }
    void SetNumDets(Int_t numberDet) { fNumDet = numberDet; }
    void SetDist2target(Float_t value, Int_t det) { fDistance2target->AddAt(value, det - 1); }
    void SetAngleTheta(Float_t value, Int_t det) { fAngleTheta->AddAt(value, det - 1); }
    void SetAnglePhi(Float_t value, Int_t det) { fAnglePhi->AddAt(value, det - 1); }
    void SetOffsetX(Float_t value, Int_t det) { fOffsetX->AddAt(value, det - 1); }
    void SetOffsetY(Float_t value, Int_t det) { fOffsetY->AddAt(value, det - 1); }
    void SetPolPar(Int_t parPol) { fPolPar = parPol; }
    void SetEnevsPosCorr(Float_t value, Int_t order) { fEnevsPosCorr->AddAt(value, order - 1); }

  private:
    Int_t fGeo;                // Geometry of the foot detectors
    Int_t fNumDet;             // Number of foot detectors (from 1 to 10 for 2021 experiments)
    TArrayF* fDistance2target; // Distance to target
    TArrayF* fAngleTheta;      // Angle with respect to beam direction
    TArrayF* fAnglePhi;        // Rotation Angle around the Z direction
    TArrayF* fOffsetX;         // Offset for x
    TArrayF* fOffsetY;         // Offset for y
    Int_t fPolPar;             // Number of parameters for polynomial for Energy Vs Position Correction
    TArrayF* fEnevsPosCorr;    // Polynomial parameters for Energy Vs Position Correction

    const R3BFootMappingPar& operator=(const R3BFootMappingPar&);
    R3BFootMappingPar(const R3BFootMappingPar&);

  public:
    ClassDef(R3BFootMappingPar, 1);
};
