/******************************************************************************
 *   Copyright (C) 2019 GSI Helmholtzzentrum für Schwerionenforschung GmbH    *
 *   Copyright (C) 2019 Members of R3B Collaboration                          *
 *                                                                            *
 *             This software is distributed under the terms of the            *
 *                 GNU General Public Licence (GPL) version 3,                *
 *                    copied verbatim in the file "LICENSE".                  *
 *                                                                            *
 * In applying this license GSI does not waive the privileges and immunities  *
 * granted to it by virtue of its status as an Intergovernmental Organization *
 * or submit itself to any jurisdiction.                                      *
 ******************************************************************************/

#ifndef R3B_OPTIMIZEMAGFIELDS494_H
#define R3B_OPTIMIZEMAGFIELDS494_H

#include "FairTask.h"

#include <string>
#include <vector>
#include "Math/Factory.h"
#include "Math/Functor.h"
#include "Math/Minimizer.h"
#include "Minuit2/Minuit2Minimizer.h"

//#define numPoints 1580      // 2023 data
//#define numPoints 1596      // 2022 data
#define numPoints 3176      // 2022+2023 data

//#define numPoints 3432

class TClonesArray;
class R3BFieldPar;

class TH1F;
class TH2F;

class R3BOptimizeMagfieldS494 : public FairTask
{
  public:
    R3BOptimizeMagfieldS494(const char* name, Bool_t vis = kFALSE, Int_t verbose = 1);
    virtual ~R3BOptimizeMagfieldS494();

    virtual InitStatus Init();
    virtual void Exec(const Option_t* = "");

    double Chi2();

    virtual void Finish();

  private:
    Bool_t InitPropagator();

    R3BFieldPar* fFieldPar;
   	Double_t totalChi2Mass = 0;
   	Double_t totalChi2P = 0;
   	Double_t d0[numPoints][6] = {0};
   	Double_t d[numPoints][6] = {0};


    ClassDef(R3BOptimizeMagfieldS494, 1)
};

#endif
