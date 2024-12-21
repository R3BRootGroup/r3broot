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

// -------------------------------------------------------------
// -----        R3BTofDMappingPar source file              -----
// -----    Created 18/03/22 by J.L. Rodriguez-Sanchez     -----
// -------------------------------------------------------------

#pragma once

#include "FairParGenericSet.h"

#include "TArrayI.h"
#include "TObjArray.h"
#include <Rtypes.h>
#include <stdint.h>
#include <vector>

class FairParamList;

class R3BTofDMappingPar : public FairParGenericSet
{
  public:
    /** Standard constructor **/
    R3BTofDMappingPar(const TString& name = "tofdMappingPar",
                      const TString& title = "TofD Mapping parameters",
                      const TString& context = "tofdMappingContext");

    /** Destructor **/
    virtual ~R3BTofDMappingPar();

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
    const Int_t GetNbPlanes() { return fNumPlanes; }
    const Int_t GetNbPaddles() { return fNumPaddles; }
    // GetTrigMap in 1-base for plane(1-4), paddle(1-44) and pmt(1-2)
    const Int_t GetTrigMap(UInt_t plane, UInt_t paddle, UInt_t pmt)
    {
        return fTrigmap[pmt - 1][plane - 1]->GetAt(paddle - 1);
    }

    void SetNbPlanes(Int_t p) { fNumPlanes = p; }
    void SetNbPaddles(Int_t p) { fNumPaddles = p; }
    // SetTrigMap in 1-base for plane(1-4), paddle(1-44) and pmt(1-2)
    void SetTrigMap(Int_t value, UInt_t plane, UInt_t paddle, UInt_t pmt)
    {
        fTrigmap[pmt - 1][plane - 1]->AddAt(value, paddle - 1);
    }

  private:
    Int_t fNumPlanes;
    Int_t fNumPaddles;
    Int_t fNumPmts;
    std::vector<TArrayI*> fTrigmap[2]; // Two PMTs per paddle

    const R3BTofDMappingPar& operator=(const R3BTofDMappingPar&);
    R3BTofDMappingPar(const R3BTofDMappingPar&);

  public:
    ClassDef(R3BTofDMappingPar, 1);
};
