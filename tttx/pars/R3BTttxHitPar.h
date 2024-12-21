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

// ------------------------------------------------------------------
// -----                  R3BTttxHitPar source file             -----
// -----         Created 18/01/24 by J.L. Rodriguez-Sanchez     -----
// ------------------------------------------------------------------

#pragma once

#include <FairParGenericSet.h>

#include <Rtypes.h>
#include <TArrayF.h>
#include <TObject.h>
#include <cassert>

class FairParamList;

class R3BTttxHitPar : public FairParGenericSet
{
  public:
    /** Standard constructor **/
    R3BTttxHitPar(const char* name = "tttxHitPar",
                  const char* title = "TTTX Hit Parameters",
                  const char* context = "TttxHitParContext");

    /** Destructor **/
    virtual ~R3BTttxHitPar();

    /** Method to reset all parameters **/
    void clear() override;

    /** Method to store all parameters using FairRuntimeDB **/
    void putParams(FairParamList* list) override;

    /** Method to retrieve all parameters using FairRuntimeDB**/
    bool getParams(FairParamList* list) override;

    /** Method to print values of parameters to the standard output **/
    void print() override;

    /** Accessor functions **/
    [[nodiscard]] inline const int GetNumDets() const { return fNumDets; }
    [[nodiscard]] inline const int GetNumParsZfit() const { return fNumParsZfit; }
    [[nodiscard]] inline TArrayF* GetZfitParams() { return fZfitpar; }

    inline void SetNumDets(int num) { fNumDets = num; }
    inline void SetNumParsZfit(int num) { fNumParsZfit = num; }
    inline void SetZfitParams(double value, int par, int det)
    {
        assert(std::clamp(det, 1, fNumDets) == det && "Number of tttx detector out of range");
        assert(std::clamp(par, 1, fNumParsZfit) == par && "Index of Z-fit parameters out of range");
        int index = (det - 1) * fNumParsZfit + par - 1;
        fZfitpar->AddAt(value, index);
    }

    /** Create more Methods if you need them! **/

  private:
    int fNumDets = 2;     // Number of detectors
    int fNumParsZfit = 2; // Number of parameters used to calibrate the charge Z
    TArrayF* fZfitpar;    // Parameters used to calibrate the charge Z

    const R3BTttxHitPar& operator=(const R3BTttxHitPar&);
    R3BTttxHitPar(const R3BTttxHitPar&);

  public:
    ClassDefOverride(R3BTttxHitPar, 1);
};
