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

#define NCHMAX 5000

class FairParamList;

/**
 * Container for time calibration of a single module. This class is used for
 * storage of FiberMAPMT hit parameter
 */

class R3BFiberMAPMTHitModulePar : public FairParGenericSet
{
  public:
    /**
     * Standard constructor.
     * @param name a name of container.
     * @param title a title of container.
     * @param context context/purpose for parameters and conditions.
     * @param own class ownership, if flag is kTRUE FairDB has the par.
     */
    R3BFiberMAPMTHitModulePar(const char* name = "FiberMAPMTHitModulePar",
                              const char* title = "Fiber hit calibration of a module",
                              const char* context = "TestDefaultContext",
                              Bool_t own = kTRUE);

    /**
     * Destructor.
     * Frees the memory allocated by the object.
     */
    virtual ~R3BFiberMAPMTHitModulePar(void);

    /**
     * A method to reset the parameter values. Sets all parameters to 0.
     */
    void clear(void);

    /**
     * A method to write parameters using RuntimeDB.
     * @param list a list of parameters.
     */
    void putParams(FairParamList* list);

    /**
     * A method to read parameters using RuntimeDB.
     * @param list a list of parameters.
     * @return kTRUE if successfull, else kFALSE.
     */
    Bool_t getParams(FairParamList* list);

    /**
     * A method to print value of parameters to the standard
     * output using FairLogger.
     */
    void printParams();

    /**
     * A method to draw vaues of parameters on the current Canvas.
     */
    void DrawParams();

    /** Accessor functions **/
    Int_t GetFiber() const { return fFiber; }
    Double_t GetSync() const { return fSync; }
    Double_t GetOffsetUp() const { return fOffsetUp; }
    Double_t GetOffsetDown() const { return fOffsetDown; }
    Double_t GetVeff() const { return fVeff; }
    Double_t GetGainUp() const { return fGainUp; }
    Double_t GetGainDown() const { return fGainDown; }

    void SetFiber(Int_t i) { fFiber = i; }
    void SetSync(Double_t t) { fSync = t; }
    void SetOffsetUp(Double_t t) { fOffsetUp = t; }
    void SetOffsetDown(Double_t t) { fOffsetDown = t; }
    void SetGainUp(Double_t gainUp) { fGainUp = gainUp; }
    void SetGainDown(Double_t gainDown) { fGainDown = gainDown; }
    void SetVeff(Double_t v) { fVeff = v; }

  private:
    Int_t fFiber;         /**< Index of a fiber. */
    Double_t fOffsetUp;   /**< time offset of PM1 */
    Double_t fOffsetDown; /**< time offset of PM2 */
    Double_t fVeff;       /**< effective velocity of light in paddle for position determination */
    Double_t fSync;       /**< time offset between paddles */
    Double_t fGainUp;
    Double_t fGainDown;

  public:
    ClassDef(R3BFiberMAPMTHitModulePar, 1);
};
