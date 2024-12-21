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
// -----                R3BFootCalPar source file               -----
// -----         Created 07/09/21 by J.L. Rodriguez-Sanchez     -----
// ------------------------------------------------------------------

#pragma once

#include <FairParGenericSet.h>

#include <Rtypes.h>
#include <TArrayF.h>
#include <TObject.h>

class FairParamList;

class R3BFootCalPar : public FairParGenericSet
{
  public:
    /** Standard constructor **/
    R3BFootCalPar(const char* name = "footStripCalPar",
                  const char* title = "FOOT Strip Parameters",
                  const char* context = "FootCalParContext");

    /** Destructor **/
    virtual ~R3BFootCalPar();

    /** Method to reset all parameters **/
    void clear() override;

    /** Method to store all parameters using FairRuntimeDB **/
    void putParams(FairParamList* list) override;

    /** Method to retrieve all parameters using FairRuntimeDB**/
    Bool_t getParams(FairParamList* list) override;

    /** Method to print values of parameters to the standard output **/
    void print() override;

    /** Accessor functions **/
    [[nodiscard]] inline const int GetNumDets() const { return fNumDets; }
    [[nodiscard]] inline const int GetNumStrips() const { return fNumStrips; }
    [[nodiscard]] inline const int GetNumParsFit() const { return fNumParsFit; }
    TArrayF* GetStripCalParams() { return fStripCalParams; }

    inline void SetNumDets(int ndet) { fNumDets = ndet; }
    inline void SetNumStrips(int nstrip) { fNumStrips = nstrip; }
    inline void SetNumParsFit(int npar) { fNumParsFit = npar; }
    inline void SetStripCalParams(float c, int i) { fStripCalParams->AddAt(c, i); }

    /** Create more Methods if you need them! **/

  private:
    int fNumDets = 16;        // Number of detectors
    int fNumStrips = 640;     // Number of strips
    int fNumParsFit = 2;      // Number of parameters for the fit
    TArrayF* fStripCalParams; // Calibration Parameters of strip

    const R3BFootCalPar& operator=(const R3BFootCalPar&);
    R3BFootCalPar(const R3BFootCalPar&);

  public:
    ClassDefOverride(R3BFootCalPar, 1);
};
