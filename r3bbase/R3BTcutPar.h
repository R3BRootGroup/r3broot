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

#ifndef R3BTCUTPAR_H
#define R3BTCUTPAR_H 1

#include "TCutG.h"
#include "TObject.h"
#include <Rtypes.h>

class FairParamList;

class R3BTcutPar : public TNamed
{
  public:
    /** Default constructor **/
    R3BTcutPar(const TString& cutname = "Cut");

    /** Destructor **/
    virtual ~R3BTcutPar() {}

    /** Method to print values of parameters **/
    void print();

    /** Method to store all parameters using FairRuntimeDB **/
    virtual void putParams(FairParamList* list);

    /** Method to retrieve all parameters using FairRuntimeDB**/
    TCutG* getParams(FairParamList* list);

    /** Accessor functions **/
    TCutG* GetCut() { return fCut; }

    void SetMaxPoints(UInt_t p) { fMaxPoints = p; }
    void SetCut(TCutG* c) { fCut = c; }
    TString GetNameObj() { return GetName(); }

  private:
    UInt_t fMaxPoints;
    TCutG* fCut;

  public:
    ClassDef(R3BTcutPar, 0);
};

#endif // R3BTCUTPAR_H
