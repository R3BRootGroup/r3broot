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
 * storage of BunchedFiber hit parameter
 */

class R3BBunchedFiberHitModulePar : public FairParGenericSet
{
  public:
    /**
     * Standard constructor.
     * @param name a name of container.
     * @param title a title of container.
     * @param context context/purpose for parameters and conditions.
     * @param own class ownership, if flag is kTRUE FairDB has the par.
     */
    R3BBunchedFiberHitModulePar(const char* name = "BunchedFiberHitModulePar",
                                const char* title = "Fiber hit calibration of a module",
                                const char* context = "TestDefaultContext",
                                Bool_t own = kTRUE);

    /**
     * Destructor.
     * Frees the memory allocated by the object.
     */
    virtual ~R3BBunchedFiberHitModulePar(void);

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
    Double_t GetOffset1() const { return fOffset1; }
    Double_t GetOffset2() const { return fOffset2; }
    Double_t GetVeff() const { return fVeff; }
    Double_t GetGainMA() const { return fGainMA; }
    Double_t GetGainS() const { return fGainS; }

    void SetFiber(Int_t i) { fFiber = i; }
    void SetSync(Double_t t) { fSync = t; }
    void SetOffset1(Double_t t) { fOffset1 = t; }
    void SetOffset2(Double_t t) { fOffset2 = t; }
    void SetGainMA(Double_t gainMA) { fGainMA = gainMA; }
    void SetGainS(Double_t gainS) { fGainS = gainS; }
    void SetVeff(Double_t v) { fVeff = v; }

  private:
    Int_t fFiber;      /**< Index of a fiber. */
    Double_t fOffset1; /**< time offset of PM1 */
    Double_t fOffset2; /**< time offset of PM2 */
    Double_t fVeff;    /**< effective velocity of light in paddle for position determination */
    Double_t fSync;    /**< time offset between paddles */
    Double_t fGainMA;
    Double_t fGainS;

  public:
    ClassDef(R3BBunchedFiberHitModulePar, 1);
};
