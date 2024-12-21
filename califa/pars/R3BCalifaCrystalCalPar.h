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

#pragma once

#include <FairParGenericSet.h>

#include <TArrayF.h>
#include <TObjArray.h>
#include <TObjString.h>
#include <TObject.h>
#include <cassert>

class FairParamList;

class R3BCalifaCrystalCalPar : public FairParGenericSet
{
  public:
    /** Standard constructor **/
    R3BCalifaCrystalCalPar(const char* name = "califaCrystalCalPar",
                           const char* title = "Califa CrystalCal Parameters",
                           const char* context = "CalifaCalParContext");

    /** Destructor **/
    virtual ~R3BCalifaCrystalCalPar();

    /** Method to reset all parameters **/
    virtual void clear();

    /** Method to store all parameters using FairRuntimeDB **/
    virtual void putParams(FairParamList* list);

    /** Method to retrieve all parameters using FairRuntimeDB**/
    Bool_t getParams(FairParamList* list);

    /** Method to print values of parameters to the standard output **/
    virtual void print();
    void printParams();
    void printCalCrystalInfo(const int cryID);

    /** Accessor functions **/
    [[nodiscard]] inline const int GetNumCrystals() const { return fNumCrystals; }
    [[nodiscard]] inline const int GetNumParametersFit() const { return fNumParamsFit; }
    [[nodiscard]] inline TArrayF* GetCryCalParams() { return fCryCalParams; }

    inline void SetNumCrystals(int numberCry) { fNumCrystals = numberCry; }
    inline void SetNumParametersFit(int numberParams) { fNumParamsFit = numberParams; }
    inline void SetCryCalParams(Float_t cc, int cry)
    {
        assert(std::clamp(cry, 1, fNumCrystals * fNumParamsFit) == cry && "Number of crystal out of range");
        fCryCalParams->AddAt(cc, cry);
    }

    /** Create more Methods if you need them! **/

  private:
    TArrayF* fCryCalParams;  /*< Calibration Parameters of Crystals>*/
    int fNumCrystals = 5088; /*< number of crystals>*/
    int fNumParamsFit = 2;   /*< number of cal parameters in the fit
                 pol1: A_fit & B_fit
                 pol2: A_fit, B_fit & C_fit>*/

    const R3BCalifaCrystalCalPar& operator=(const R3BCalifaCrystalCalPar&);
    R3BCalifaCrystalCalPar(const R3BCalifaCrystalCalPar&);

  public:
    ClassDef(R3BCalifaCrystalCalPar, 1);
};
