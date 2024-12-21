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
// -----                R3BTttxStripCalPar source file          -----
// -----         Created 18/01/24 by J.L. Rodriguez-Sanchez     -----
// ------------------------------------------------------------------

#pragma once

#include <FairParGenericSet.h>

#include <Rtypes.h>
#include <TArrayF.h>
#include <TArrayI.h>
#include <TObject.h>
#include <cassert>

class FairParamList;

class R3BTttxStripCalPar : public FairParGenericSet
{
  public:
    /** Standard constructor **/
    R3BTttxStripCalPar(const char* name = "tttxStripCalPar",
                       const char* title = "TTTX Strip Parameters",
                       const char* context = "TttxCalParContext");

    /** Destructor **/
    virtual ~R3BTttxStripCalPar();

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
    [[nodiscard]] inline const int GetNumStrips() const { return fNumStrips; }
    [[nodiscard]] inline const int GetNumParsFit() const { return fNumParsFit; }
    [[nodiscard]] inline TArrayF* GetStripCalParams() { return fStripCalParams; }
    [[nodiscard]] inline int GetInUse(int det, int strip)
    {
        assert(std::clamp(det, 1, fNumDets) == det && "Number of detectors out of range");
        assert(std::clamp(strip, 1, fNumStrips) == strip && "Number of strips out of range");
        return fIn_use->GetAt(fNumStrips * (det - 1) + strip - 1);
    }

    inline void SetNumDets(int n) { fNumDets = n; }
    inline void SetNumStrips(int n) { fNumStrips = n; }
    inline void SetNumParsFit(int n) { fNumParsFit = n; }
    inline void SetStripCalParams(double value, int det, int strip)
    {
        assert(std::clamp(det, 1, fNumDets) == det && "Number of detectors out of range");
        assert(std::clamp(strip, 1, fNumStrips) == strip && "Number of strip out of range");
        fStripCalParams->AddAt(value, fNumStrips * (det - 1) + strip - 1);
    }
    inline void SetInUse(int value, int det, int strip)
    {
        assert(std::clamp(det, 1, fNumDets) == det && "Number of detectors out of range");
        assert(std::clamp(strip, 1, fNumStrips) == strip && "Number of strip out of range");
        fIn_use->AddAt(value, fNumStrips * (det - 1) + strip - 1);
    }

  private:
    TArrayF* fStripCalParams; // Calibration Parameters of strip
    TArrayI* fIn_use;         // 1: in use
    int fNumDets = 2;         // Number of detectors
    int fNumStrips = 32;      // Number of strips
    int fNumParsFit = 2;      // Number of parameters for the fit

    const R3BTttxStripCalPar& operator=(const R3BTttxStripCalPar&);
    R3BTttxStripCalPar(const R3BTttxStripCalPar&);

  public:
    ClassDefOverride(R3BTttxStripCalPar, 1);
};
