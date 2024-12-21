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

#ifndef R3BNEULANDHITPAR_H
#define R3BNEULANDHITPAR_H

#include "FairParGenericSet.h" // for FairParGenericSet
#include "R3BNeulandHitModulePar.h"
#include "TObjArray.h"

#include <array>

class FairParamList;

class R3BNeulandHitPar : public FairParGenericSet
{

  public:
    /**
     * Standard constructor. Creates instance of this class.
     * @param name a name of container.
     * @param title a title of container.
     * @param context context/purpose for parameters and conditions.
     * @param own class ownership, if flag is kTRUE FairDB has the par.
     */
    R3BNeulandHitPar(const char* name = "NeulandHitPar",
                     const char* title = "Calibration",
                     const char* context = "TestDefaultContext",
                     Bool_t own = kTRUE);
    /**
     * Destructor. Cleares the memory used by the object.
     */
    virtual ~R3BNeulandHitPar(void);

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
     */
    void printParams();

    /**
     * Method to add parameter container for a module.
     * Extends the array.
     * @param tch a parameter container for a detector module.
     */
    void AddModulePar(R3BNeulandHitModulePar* tch) { fParams->Add(tch); }

    /**
     * Method to retrieve the arrray with module containers.
     * @return an array with parameter containers of type R3BNeulandHitModulePar.
     */
    TObjArray* GetListOfModulePar() const { return fParams; }

    /**
     * Method to get number of modules storred in array.
     * @return size of array.
     */
    Int_t GetNumModulePar() const { return fParams->GetEntriesFast(); }

    /**
     * Method to get single parameter container for a specific module.
     * @param idx an index of a module.
     * @return parameter container of this module.
     */
    R3BNeulandHitModulePar* GetModuleParAt(Int_t idx) const
    {
        return dynamic_cast<R3BNeulandHitModulePar*>(fParams->At(idx));
    }

    // Global time offset in ns
    inline Double_t GetGlobalTimeOffset() const { return fGlobalTimeOffset; }
    inline void SetGlobalTimeOffset(Double_t t0_ns) { fGlobalTimeOffset = t0_ns; }

    // Distance to Target in cm
    inline Double_t GetDistanceToTarget(Int_t plane = 0) const
    {
        return fDistanceToTarget + fDistancesToFirstPlane.at(plane);
    }
    inline void SetDistanceToTarget(Double_t distance_cm) { fDistanceToTarget = distance_cm; }

    // Energy Cutoff for all Hits in MeV
    inline Double_t GetEnergyCutoff() const { return fEnergyCut; }
    inline void SetEnergyCutoff(Double_t energy_MeV) { fEnergyCut = energy_MeV; }

    // Number of Planes
    inline Int_t GetNumberOfPlanes() const { return fDistancesToFirstPlane.size(); }
    inline void SetNumberOfPlanes(const Int_t nPlanes) { fDistancesToFirstPlane.resize(nPlanes); }

    // Distance to first Plane in cm
    inline Double_t GetDistanceToFirstPlane(Int_t plane) const { return fDistancesToFirstPlane.at(plane); }
    inline void SetDistanceToFirstPlane(Int_t plane, Double_t distance_cm)
    {
        fDistancesToFirstPlane.at(plane) = distance_cm;
    }

    inline std::vector<Double_t> GetDistancesToFirstPlane() const { return fDistancesToFirstPlane; }

  private:
    TObjArray* fParams; /**< an array with parameter containers of all modules */

    Double_t fGlobalTimeOffset;
    Double_t fDistanceToTarget;
    Double_t fEnergyCut;
    std::vector<Double_t> fDistancesToFirstPlane;

    ClassDef(R3BNeulandHitPar, 3);
};

#endif /* !R3BNEULANDHITPAR_H*/
