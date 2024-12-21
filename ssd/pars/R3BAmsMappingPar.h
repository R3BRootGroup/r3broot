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
// -----        R3BAmsMappingPar source file                  -----
// -----    Created 20/01/20  by J.L. Rodriguez-Sanchez       -----
// ----------------------------------------------------------------

#ifndef R3BAmsMappingPar_H
#define R3BAmsMappingPar_H

#include "FairParGenericSet.h"
#include "TObject.h"

#include "TArrayF.h"
#include "TArrayI.h"
#include "TObjArray.h"
#include <TObjString.h>

class FairParamList;

class R3BAmsMappingPar : public FairParGenericSet
{
  public:
    /** Standard constructor **/
    R3BAmsMappingPar(const TString& name = "amsMappingPar",
                     const TString& title = "AMS Mapping Parameters",
                     const TString& context = "AmsMappingContext");

    /** Destructor **/
    virtual ~R3BAmsMappingPar();

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
    const Int_t GetSam(Int_t det) { return fSam->GetAt(det - 1); }
    const Int_t GetGtb(Int_t det) { return fGtb->GetAt(det - 1); }
    const Int_t GetSiderem(Int_t det) { return fSiderem->GetAt(det - 1); }
    const Float_t GetDist2target(Int_t det) { return fDistance2target->GetAt(det - 1); }
    const Float_t GetAngleTheta(Int_t det) { return fAngleTheta->GetAt(det - 1); }
    const Float_t GetOffsetY(Int_t det) { return fOffsetY->GetAt(det - 1); }

    void SetGeometry(Int_t geo) { fGeo = geo; }
    void SetNumDets(Int_t numberDet) { fNumDet = numberDet; }
    void SetSam(Int_t value, Int_t det) { fSam->AddAt(value, det - 1); }
    void SetGtb(Int_t value, Int_t det) { fGtb->AddAt(value, det - 1); }
    void SetSiderem(Int_t value, Int_t det) { fSiderem->AddAt(value, det - 1); }
    void SetDist2target(Float_t value, Int_t det) { fDistance2target->AddAt(value, det - 1); }
    void SetAngleTheta(Float_t value, Int_t det) { fAngleTheta->AddAt(value, det - 1); }
    void SetOffsetY(Float_t value, Int_t det) { fOffsetY->AddAt(value, det - 1); }

  private:
    Int_t fGeo;                // geometry of the ams detectors per year: 2019 or 2020
    Int_t fNumDet;             // number of ams detectors (from 1 to 6 for 2020 experiments)
    TArrayI* fSam;             // sam number (4 or 5)
    TArrayI* fGtb;             // gtb number (0 or 1)
    TArrayI* fSiderem;         // siderem number (1 or 2)
    TArrayF* fDistance2target; // Distance to target
    TArrayF* fAngleTheta;      // Angle with respect to beam direction (-45 and 45 degrees)
    TArrayF* fOffsetY;         // offset for y (-20, 0 or 20 mm)

    const R3BAmsMappingPar& operator=(const R3BAmsMappingPar&);
    R3BAmsMappingPar(const R3BAmsMappingPar&);

  public:
    ClassDef(R3BAmsMappingPar, 1);
};

#endif /* R3BAmsMappingPar_H */
