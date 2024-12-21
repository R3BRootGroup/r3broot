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
// -----        R3BFiberMappingPar source file             -----
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

class R3BFiberMappingPar : public FairParGenericSet
{
  public:
    /** Standard constructor **/
    R3BFiberMappingPar(const TString& name = "fiberMappingPar",
                       const TString& title = "Fiber Mapping parameters",
                       const TString& context = "fiberMappingContext");

    /** Destructor **/
    virtual ~R3BFiberMappingPar();

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
    const Int_t GetNbChannels() { return fNbChannels; }
    const Int_t GetNbSides() { return fNbSides; }
    // GetTrigMap in 1-base for side(1-2) and channel(1-X)
    const Int_t GetTrigMap(UInt_t side, UInt_t ch) { return fTrigmap[side - 1]->GetAt(ch - 1); }

    void SetNbChannels(Int_t p) { fNbChannels = p; }
    void SetNbSides(Int_t p) { fNbSides = p; }
    // SetTrigMap in 1-base for side(1-2) and channel(1-X)
    void SetTrigMap(Int_t value, UInt_t side, UInt_t ch) { fTrigmap[side - 1]->AddAt(value, ch - 1); }

  private:
    Int_t fNbChannels;
    Int_t fNbSides;
    TArrayI* fTrigmap[2]; // Two sides per fiber

    const R3BFiberMappingPar& operator=(const R3BFiberMappingPar&);
    R3BFiberMappingPar(const R3BFiberMappingPar&);

  public:
    ClassDef(R3BFiberMappingPar, 1);
};
