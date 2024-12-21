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

#ifndef R3BTSPLINEPAR_H
#define R3BTSPLINEPAR_H 1

#include "TObject.h"
#include "TSpline.h"
#include <Rtypes.h>

class FairParamList;

class R3BTsplinePar : public TNamed
{
  public:
    /** Default constructor **/
    R3BTsplinePar(const TString& name = "Spline");

    /** Destructor **/
    virtual ~R3BTsplinePar() {}

    /** Method to print values of parameters **/
    void print();

    /** Method to store all parameters using FairRuntimeDB **/
    virtual void putParams(FairParamList* list);

    /** Method to retrieve all parameters using FairRuntimeDB**/
    TSpline3* getParams(FairParamList* list);

    /** Accessor functions **/
    TString GetNameObj() { return GetName(); }
    TSpline3* GetSpline() { return fSpline; }

    void SetMaxPoints(UInt_t p) { fMaxPoints = p; }
    void SetSpline(TSpline3* c) { fSpline = c; }

  private:
    UInt_t fMaxPoints;
    TSpline3* fSpline;
    Double_t fxmin;
    Double_t fxmax;

  public:
    ClassDef(R3BTsplinePar, 0);
};

#endif // R3BTSPLINEPAR_H
