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
#include <TArrayD.h>
#include <memory>

class FairParamList;

class R3BFrsSciTcalPar : public FairParGenericSet
{
  public:
    /** Standard constructor **/
    R3BFrsSciTcalPar(const char* name = "FrsSciTcalPar",
                     const char* title = "FrsSciTcal Parameters",
                     const char* context = "FrsSciTcalParContext");

    /** Destructor **/
    virtual ~R3BFrsSciTcalPar();

    /** Method to reset all parameters **/
    virtual void clear();

    /** Method to store all parameters using FairRuntimeDB **/
    virtual void putParams(FairParamList* list);

    /** Method to retrieve all parameters using FairRuntimeDB**/
    Bool_t getParams(FairParamList* list);

    /** Method to print values of parameters to the standard output **/
    void printParams();

    /** Accessor functions **/

    inline void SetNumDets(int nDets = 3) { fNumDets = nDets; }
    inline void SetNumPmts(int nPmts = 3) { fNumPmts = nPmts; }
    inline void SetNumPars(int num = 1000) { fNumPars = num; }
    inline void SetMinStat(int stat = 0) { fMinStat = stat; }
    inline void SetOneTcalParam(Double_t ft_ns, UInt_t rank) { fAllTcalParams->AddAt(ft_ns, rank); }

    [[nodiscard]] inline const int GetNumDets() const { return fNumDets; }
    [[nodiscard]] inline const int GetNumPmts() const { return fNumPmts; }
    [[nodiscard]] inline const int GetNumPars() const { return fNumPars; }
    [[nodiscard]] inline const int GetMinStat() const { return fMinStat; }

    const TArrayD& GetAllTcalParams() const { return *fAllTcalParams; }
    Double_t GetOneTcalParam(UInt_t rank) { return (Double_t)fAllTcalParams->GetAt(rank); }

  private:
    int fNumDets = 3;    // number of FrsSci detectors
    int fNumPmts = 3;    // number of Pmts  (=3)
    int fNumPars = 1000; // 1000 parameters per signal for VFTX calibration
    int fMinStat = 0;
    std::unique_ptr<TArrayD> fAllTcalParams; // Calibration Parameters for all signals of all detectors

    const R3BFrsSciTcalPar& operator=(const R3BFrsSciTcalPar&);
    R3BFrsSciTcalPar(const R3BFrsSciTcalPar&);

    ClassDef(R3BFrsSciTcalPar, 1);
};
