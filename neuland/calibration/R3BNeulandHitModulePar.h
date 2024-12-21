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

#ifndef R3BNEULANDHITMODULEPAR_H
#define R3BNEULANDHITMODULEPAR_H

#include "FairParGenericSet.h"

#include "R3BNeulandCommon.h"

#include "assert.h"

// FIXME Do we need this???
// #define NCHMAX 5000

class FairParamList;

class R3BNeulandHitModulePar : public FairParGenericSet
{
  public:
    /**
     * Standard constructor.
     * @param name a name of container.
     * @param title a title of container.
     * @param context context/purpose for parameters and conditions.
     * @param own class ownership, if flag is kTRUE FairDB has the par.
     */
    R3BNeulandHitModulePar(const char* name = "NeulandHitModulePar",
                           const char* title = "Calibration of a NeuLAND module",
                           const char* context = "TestDefaultContext",
                           Bool_t own = kTRUE);

    /**
     * Destructor.
     * Frees the memory allocated by the object.
     */
    virtual ~R3BNeulandHitModulePar(void);

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

    /** Accessor functions **/
    Int_t GetModuleId() const { return fModuleId; }
    Double_t GetTDiff() const { return fTDiff; }
    Double_t GetTSync() const { return fTSync; }
    Double_t GetTimeOffset(const Int_t side) const
    {
        assert(side == 1 || side == 2);
        return -fTSync + (1.5 - side) * fTDiff;
    }
    Double_t GetEffectiveSpeed() const { return fEffectiveSpeed; }
    Int_t GetPedestal(const Int_t side) const
    {
        assert(side == 1 || side == 2);
        return fPedestal[side - 1];
    }
    Double_t GetEnergyGain(const Int_t side) const
    {
        assert(side == 1 || side == 2);
        return fEnergyGain[side - 1];
    }
    Double_t GetLightAttenuationLength() const { return fLightAttenuationLength; }
    Double_t GetPMTSaturation(const Int_t side) const
    {
        assert(side == 1 || side == 2);
        return fPMTSaturation[side - 1];
    }
    Double_t GetPMTThreshold(const Int_t side) const
    {
        assert(side == 1 || side == 2);
        return fPMTThreshold[side - 1];
    }

    void SetModuleId(const Int_t i) { fModuleId = i; }
    void SetTDiff(const Double_t val) { fTDiff = val; }
    void SetTSync(const Double_t val) { fTSync = val; }
    void SetEffectiveSpeed(const Double_t val) { fEffectiveSpeed = val; }
    void SetPedestal(const Int_t val, const Int_t side)
    {
        assert(side == 1 || side == 2);
        fPedestal[side - 1] = val;
    }
    void SetEnergyGain(const Double_t val, const Int_t side)
    {
        assert(side == 1 || side == 2);
        fEnergyGain[side - 1] = val;
    }
    void SetLightAttenuationLength(const Double_t val) { fLightAttenuationLength = val; }
    void SetPMTSaturation(const Double_t val, const Int_t side)
    {
        assert(side == 1 || side == 2);
        fPMTSaturation[side - 1] = val;
    }
    void SetPMTThreshold(const Double_t val, const Int_t side)
    {
        assert(side == 1 || side == 2);
        fPMTThreshold[side - 1] = val;
    }

  private:
    Int_t fModuleId; /**< Index of a detector module. */
    Double_t fTDiff;
    Double_t fTSync;
    Double_t fEffectiveSpeed;
    Int_t fPedestal[2];      // R/B, L/T
    Double_t fEnergyGain[2]; // R/B, L/T
    Double_t fLightAttenuationLength;
    Double_t fPMTSaturation[2]; // R/B, L/T
    Double_t fPMTThreshold[2];  // R/B, L/T
    ClassDef(R3BNeulandHitModulePar, 3);
};

#endif /* !R3BNEULANDHITMODULEPAR_H*/
