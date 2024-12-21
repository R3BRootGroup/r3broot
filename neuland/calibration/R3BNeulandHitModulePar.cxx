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

#include "R3BNeulandHitModulePar.h"

#include "FairLogger.h"
#include "FairParamList.h" // for FairParamList

#include "TString.h"

R3BNeulandHitModulePar::R3BNeulandHitModulePar(const char* name, const char* title, const char* context, Bool_t own)
    : FairParGenericSet(name, title, context, own)
    , fModuleId(0)
{
    // Reset all parameters
    clear();
}

R3BNeulandHitModulePar::~R3BNeulandHitModulePar() {}

void R3BNeulandHitModulePar::putParams(FairParamList* list)
{
    LOG(info) << "R3BNeulandHitModulePar::putParams() called";
    if (!list)
    {
        return;
    }
    list->add("module_id", fModuleId);
    list->add("tdiff", fTDiff);
    list->add("tsync", fTSync);
    list->add("pedestal_1", fPedestal[0]);
    list->add("pedestal_2", fPedestal[1]);
    list->add("energygain_1", fEnergyGain[0]);
    list->add("energygain_2", fEnergyGain[1]);
    list->add("effectivespeed", fEffectiveSpeed);
    list->add("lightattenuationlength", fLightAttenuationLength);
    list->add("pmtsaturation_1", fPMTSaturation[0]);
    list->add("pmtsaturation_2", fPMTSaturation[1]);
    list->add("pmtthreshold_1", fPMTThreshold[0]);
    list->add("pmtthreshold_2", fPMTThreshold[1]);
}

Bool_t R3BNeulandHitModulePar::getParams(FairParamList* list)
{
    if (!list)
    {
        return kFALSE;
    }
    if (!list->fill("module_id", &fModuleId))
    {
        return kFALSE;
    }
    if (!list->fill("tdiff", &fTDiff))
    {
        return kFALSE;
    }
    if (!list->fill("tsync", &fTSync))
    {
        return kFALSE;
    }
    if (!list->fill("pedestal_1", &fPedestal[0]))
    {
        return kFALSE;
    }
    if (!list->fill("pedestal_2", &fPedestal[1]))
    {
        return kFALSE;
    }
    if (!list->fill("energygain_1", &fEnergyGain[0]))
    {
        return kFALSE;
    }
    if (!list->fill("energygain_2", &fEnergyGain[1]))
    {
        return kFALSE;
    }
    if (!list->fill("effectivespeed", &fEffectiveSpeed))
    {
        return kFALSE;
    }
    if (!list->fill("lightattenuationlength", &fLightAttenuationLength))
    {
        return kFALSE;
    }
    if (!list->fill("pmtsaturation_1", &fPMTSaturation[0]))
    {
        return kFALSE;
    }
    if (!list->fill("pmtsaturation_2", &fPMTSaturation[1]))
    {
        return kFALSE;
    }
    if (!list->fill("pmtthreshold_1", &fPMTThreshold[0]))
    {
        return kFALSE;
    }
    if (!list->fill("pmtthreshold_2", &fPMTThreshold[1]))
    {
        return kFALSE;
    }

    return kTRUE;
}

void R3BNeulandHitModulePar::clear() {}

void R3BNeulandHitModulePar::printParams()
{
    LOG(info) << "R3BNeulandHitModulePar: Calibration Parameters:";
    LOG(info) << "   Bar ID: " << fModuleId;
    LOG(info) << TString::Format("   TDiff:     %4.2f", fTDiff);
    LOG(info) << TString::Format("   TSync:     %4.2f", fTSync);
    LOG(info) << TString::Format("   Effective Speed: %4.2f", fEffectiveSpeed);
    LOG(info) << TString::Format("   Light Attenuation Length: %4.1f", fLightAttenuationLength);
    LOG(info) << TString::Format("   Energy Gain R/B:    %4.2f", fEnergyGain[0]);
    LOG(info) << TString::Format("   Energy Gain L/T:    %4.2f\n", fEnergyGain[1]);
    LOG(info) << TString::Format("   Pedestal R/B:     %4d\n", fPedestal[0]);
    LOG(info) << TString::Format("   Pedestal L/T:     %4d\n", fPedestal[1]);
    LOG(info) << TString::Format("   PMT Saturation R/B:  %4.2f\n", fPMTSaturation[0]);
    LOG(info) << TString::Format("   PMT Saturation L/T:  %4.2f\n", fPMTSaturation[1]);
    LOG(info) << TString::Format("   PMT Threshold R/B:   %4.2f\n", fPMTThreshold[0]);
    LOG(info) << TString::Format("   PMT Threshold L/T:   %4.2f\n", fPMTThreshold[1]);
}

ClassImp(R3BNeulandHitModulePar);
