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

#ifndef R3BTOFDHITMODULEPAR_H
#define R3BTOFDHITMODULEPAR_H

#include "FairParGenericSet.h"

#define NCHMAX 5000

class FairParamList;

/**
 * Container for time calibration of a single module. This class is used for
 * storage of Tofd hit parameter
 */

class R3BTofdHitModulePar : public FairParGenericSet
{
  public:
    /**
     * Standard constructor.
     * @param name a name of container.
     * @param title a title of container.
     * @param context context/purpose for parameters and conditions.
     * @param own class ownership, if flag is kTRUE FairDB has the par.
     */
    R3BTofdHitModulePar(const char* name = "TofdHitModulePar",
                        const char* title = "TOFd hit calibration of a module",
                        const char* context = "TestDefaultContext",
                        Bool_t own = kTRUE);

    /**
     * Destructor.
     * Frees the memory allocated by the object.
     */
    virtual ~R3BTofdHitModulePar(void);

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
    Int_t GetPlane() const { return fPlane; }
    Int_t GetPaddle() const { return fPaddle; }
    Double_t GetSync() const { return fSync; }
    Double_t GetOffset1() const { return fOffset1; }
    Double_t GetOffset2() const { return fOffset2; }
    Double_t GetToTOffset1() const { return fToTOffset1; }
    Double_t GetToTOffset2() const { return fToTOffset2; }
    Double_t GetVeff() const { return fVeff; }
    Double_t GetLambda() const { return fLambda; }
    Double_t GetPar1a() const { return fPar1a; }
    Double_t GetPar1b() const { return fPar1b; }
    Double_t GetPar1c() const { return fPar1c; }
    Double_t GetPar1d() const { return fPar1d; }
    Double_t GetPar2a() const { return fPar2a; }
    Double_t GetPar2b() const { return fPar2b; }
    Double_t GetPar2c() const { return fPar2c; }
    Double_t GetPar2d() const { return fPar2d; }
    Double_t GetPola() const { return fPola; }
    Double_t GetPolb() const { return fPolb; }
    Double_t GetPolc() const { return fPolc; }
    Double_t GetPold() const { return fPold; }
    Double_t GetPar1za() const { return fPar1za; }
    Double_t GetPar1zb() const { return fPar1zb; }
    Double_t GetPar1zc() const { return fPar1zc; }
    Double_t GetPar1Walk() const { return fPar1walk; }
    Double_t GetPar2Walk() const { return fPar2walk; }
    Double_t GetPar3Walk() const { return fPar3walk; }
    Double_t GetPar4Walk() const { return fPar4walk; }
    Double_t GetPar5Walk() const { return fPar5walk; }
    void SetPlane(Int_t i) { fPlane = i; }
    void SetPaddle(Int_t i) { fPaddle = i; }
    void SetSync(Double_t t) { fSync = t; }
    void SetOffset1(Double_t t) { fOffset1 = t; }
    void SetOffset2(Double_t t) { fOffset2 = t; }
    void SetToTOffset1(Double_t t) { fToTOffset1 = t; }
    void SetToTOffset2(Double_t t) { fToTOffset2 = t; }
    void SetVeff(Double_t v) { fVeff = v; }
    void SetLambda(Double_t v) { fLambda = v; }
    void SetPar1a(Double_t par1a) { fPar1a = par1a; }
    void SetPar1b(Double_t par1b) { fPar1b = par1b; }
    void SetPar1c(Double_t par1c) { fPar1c = par1c; }
    void SetPar1d(Double_t par1d) { fPar1d = par1d; }
    void SetPar2a(Double_t par2a) { fPar2a = par2a; }
    void SetPar2b(Double_t par2b) { fPar2b = par2b; }
    void SetPar2c(Double_t par2c) { fPar2c = par2c; }
    void SetPar2d(Double_t par2d) { fPar2d = par2d; }
    void SetPola(Double_t pola) { fPola = pola; }
    void SetPolb(Double_t polb) { fPolb = polb; }
    void SetPolc(Double_t polc) { fPolc = polc; }
    void SetPold(Double_t pold) { fPold = pold; }
    void SetPar1za(Double_t par1za) { fPar1za = par1za; }
    void SetPar1zb(Double_t par1zb) { fPar1zb = par1zb; }
    void SetPar1zc(Double_t par1zc) { fPar1zc = par1zc; }
    void SetPar1Walk(Double_t par1Walk) { fPar1walk = par1Walk; }
    void SetPar2Walk(Double_t par2Walk) { fPar2walk = par2Walk; }
    void SetPar3Walk(Double_t par3Walk) { fPar3walk = par3Walk; }
    void SetPar4Walk(Double_t par4Walk) { fPar4walk = par4Walk; }
    void SetPar5Walk(Double_t par5Walk) { fPar5walk = par5Walk; }

  private:
    Int_t fPlane;                            /**< Index of a plane. */
    Int_t fPaddle;                           /**< Index of a paddle. */
    Double_t fOffset1;                       /**< time offset of PM1 */
    Double_t fOffset2;                       /**< time offset of PM2 */
    Double_t fToTOffset1;                    /**< ToT offset of PM1 */
    Double_t fToTOffset2;                    /**< ToT offset of PM2 */
    Double_t fVeff;                          /**< effective velocity of light in paddle for position determination */
    Double_t fLambda;                        /**< light attenuation in paddle for position determination */
    Double_t fSync;                          /**< time offset between paddles */
    Double_t fPar1a, fPar1b, fPar1c, fPar1d; /**< double exponential parameters for position dependent charge */
    Double_t fPar2a, fPar2b, fPar2c, fPar2d; /**< double exponential parameters for position dependent charge */
    Double_t fPola, fPolb, fPolc, fPold;     /**< polynomial parameters for position dependent charge */
    Double_t fPar1za, fPar1zb, fPar1zc, fPar1zd;                    /**< quench correction parameters */
    Double_t fPar1walk, fPar2walk, fPar3walk, fPar4walk, fPar5walk; /**< walk correction parameters */

    ClassDef(R3BTofdHitModulePar, 1);
};

#endif
