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

#ifndef R3BPDCHITPAR_H
#define R3BPDCHITPAR_H


#define N_WIRE_MAX 144
#define N_PLANE_MAX_PDC 4

#define NMODULEMAX 6000

#include "FairParGenericSet.h" // for FairParGenericSet
#include "R3BPdcHitModulePar.h"
#include "TObjArray.h"
#include <map>

using namespace std;

class FairParamList;

/**
 * Parameter container for storing time calibration of a detector.
 * This class contains array of time calibration containers for each
 * module (of type R3BTCalModulePar). Instance of this class has to be
 * created using FairRuntimeDB::getContainer("name") method. Supported
 * names: LandTCalPar, LosTCalPar.
 * @author D. Kresan
 * @since September 3, 2015
 */
class R3BPdcHitPar : public FairParGenericSet
{

  public:
    /**
     * Standard constructor. Creates instance of this class.
     * @param name a name of container.
     * @param title a title of container.
     * @param context context/purpose for parameters and conditions.
     * @param own class ownership, if flag is kTRUE FairDB has the par.
     */
    R3BPdcHitPar(const char* name = "PdcHitPar",
                          const char* title = "Pdc Hit Parameters",
                          const char* context = "TestDefaultContext",
                          Bool_t own = kTRUE);
    /** 
     * Destructor. Cleares the memory used by the object.
     */
    virtual ~R3BPdcHitPar(void);

    /**
     * Method to reset the values of parameters.
     */
    void clear(void);

    /**
     * Method to store parameters using FairRuntimeDB.
     * @param list a list of parameters.
     */
    void putParams(FairParamList* list);

    /**
     * Method to retrieve parameters using FairRuntimeDB.
     * @param list a list of parameters.
     * @return kTRUE if successful, else kFALSE.
     */
    Bool_t getParams(FairParamList* list);

    /**
     * Method to print value of parameters to the standard output.
     * Calls printParams() for every module container.
     * Accessible from Context menu in TBrowser.
     */
    void printParams(); // *MENU*

    /**
     * Method to print value of parameters for a specific module.
     * @param plane an index of detector plane
     * @param paddle a paddle index within the plane
     * @param side a side of a paddle
     * Accessible from Context menu in TBrowser.
     */
    void PrintModuleParams(Int_t wire); // *MENU*

    /**
     * Method to draw value of parameters for a specific module on a current Canvas.
     * @param plane an index of detector plane
     * @param paddle a paddle index within the plane
     * @param side a side of a paddle
     * Accessible from Context menu in TBrowser.
     */
    void DrawModuleParams(Int_t wire); // *MENU*

    /**
     * Method to add parameter container for a module.
     * Extends the array.
     * @param tch a parameter container for a detector module.
     */
    void AddModulePar(R3BPdcHitModulePar* tch);

    /**
     * Method to retrieve the arrray with module containers.
     * @return an array with parameter containers of type R3BTCalModulePar.
     */
    TObjArray* GetListOfModulePar() { return fHitParams; }

    /**
     * Method to get number of modules storred in array.
     * @return size of array.
     */
    Int_t GetNumModulePar() { return fHitParams->GetEntries(); }

    /**
     * Method to get single parameter container for a specific module.
     * @param idx an index of a module.
     * @return parameter container of this module.
     */
    R3BPdcHitModulePar* GetModuleParAt(Int_t wire);

  private:
    const R3BPdcHitPar& operator=(const R3BPdcHitPar&); /**< an assignment operator */
    R3BPdcHitPar(const R3BPdcHitPar&);                  /**< a copy constructor */

    TObjArray* fHitParams; /**< an array with parameter containers of all modules */

    Bool_t fMapInit;             /**< a boolean flag for indication whether the indexing map is initialized */
    map<Int_t, Int_t> fIndexMap; /**< a map between index of a container in array and plane,paddle,side */

    ClassDef(R3BPdcHitPar, 1);
};

#endif
