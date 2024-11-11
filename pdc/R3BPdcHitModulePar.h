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

#ifndef R3BPDCHITMODULEPAR_H
#define R3BPDCHITMODULEPAR_H

#include "FairParGenericSet.h"
#include <iostream>
#include <iomanip>

using namespace std;

#define NCHMAX 5000
#define XTC_MAXPOINTS 300

class FairParamList;

/**
 * Container for time calibration of a single module. This class is used for
 * storage of Pdc hit parameter
 */

class R3BPdcHitModulePar : public FairParGenericSet
{
  public:
    /**
     * Standard constructor.
     * @param name a name of container.
     * @param title a title of container.
     * @param context context/purpose for parameters and conditions.
     * @param own class ownership, if flag is kTRUE FairDB has the par.
     */
    R3BPdcHitModulePar(const char* name = "PdcHitModulePar",
                                const char* title = "Pdc hit calibration of a module",
                                const char* context = "TestDefaultContext",
                                Bool_t own = kTRUE);

    /**
     * Destructor.
     * Frees the memory allocated by the object.
     */
    virtual ~R3BPdcHitModulePar(void);

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
    Int_t GetWire() const { return fWire; }
    Double_t GetdTmin() const { return fdTmin; }
    Double_t GetdTmax() const { return fdTmax; }
    Double_t GetTmin() const { return fTmin; }
    Double_t GetTmax() const { return fTmax; }
    Double_t GetSync() const { return fSync; }
    Double_t GetNPoints() const {return fnpoints;}
    
    const Double_t* GetXT_xArray() const ;//{return fXT_xArray;}
    const Double_t* GetXT_tArray() const ;//{return fXT_tArray;}
    

    void SetWire(Int_t i) { fWire = i; }
    void SetdTmin(Double_t dtmin) { fdTmin = dtmin; }
    void SetdTmax(Double_t dtmax) { fdTmax = dtmax; }
    void SetTmin(Double_t tmin) { fTmin = tmin; }
    void SetTmax(Double_t tmax) { fTmax = tmax; }
    void SetSync(Double_t tsync) { fSync = tsync; }
    void SetNPoints(UInt_t npoints) {fnpoints = npoints;}
    void SetXT_xArray(Double_t* fXT_xArray, Int_t isize);
  /*  {
		for(Int_t i=0;i<isize;i++){
			fXT_xArray[i]=xt_xarray[i];
		}
	}*/
	void SetXT_tArray(Double_t* fXT_tArray, Int_t isize);
  /*  {
		for(Int_t i=0;i<isize;i++){
			fXT_tArray[i]=xt_tarray[i];
		}
	}*/
  private:
    Int_t fWire;      /**< Index of a fiber. */
    Double_t fdTmin;    /**< minimum drift time */
    Double_t fdTmax;    /**< maximum drift time */
    Double_t fTmin;    /**< minimum drift time */
    Double_t fTmax;    /**< maximum drift time */
	Double_t fSync;   /** TSync */
	UInt_t fnpoints;
    Double_t fXT_tArray[XTC_MAXPOINTS]={0};
    Double_t fXT_xArray[XTC_MAXPOINTS]={0};

    ClassDef(R3BPdcHitModulePar, 1);
};

#endif
