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

#include "R3BCalifaMappingPar.h"

#include "FairLogger.h"
#include "FairParamList.h"

#include "TMath.h"
#include "TString.h"

#include <iostream>

// ---- Standard Constructor ---------------------------------------------------
R3BCalifaMappingPar::R3BCalifaMappingPar(const char* name, const char* title, const char* context)
    : FairParGenericSet(name, title, context)
    , fNumCrystals(5088)
{
    fHalf = new TArrayI(fNumCrystals);
    fRing = new TArrayI(fNumCrystals);
    fPreamp = new TArrayI(fNumCrystals);
    fChannel = new TArrayI(fNumCrystals);
    fCrystal_type = new TArrayI(fNumCrystals);
    fApd_number = new TArrayI(fNumCrystals);
    fVoltage = new TArrayF(fNumCrystals);
    fFebex_pc = new TArrayI(fNumCrystals);
    fFebex_slot = new TArrayI(fNumCrystals);
    fFebex_mod = new TArrayI(fNumCrystals);
    fFebex_channel = new TArrayI(fNumCrystals);
    fLab = new TArrayI(fNumCrystals);
    fMrcc_module = new TArrayI(fNumCrystals);
    fMrcc_bus = new TArrayI(fNumCrystals);
    fMrcc_preamp = new TArrayI(fNumCrystals);
    fIn_use = new TArrayI(fNumCrystals);
}

// ----  Destructor ------------------------------------------------------------
R3BCalifaMappingPar::~R3BCalifaMappingPar()
{
    if (fHalf)
        delete fHalf;
    if (fRing)
        delete fRing;
    if (fPreamp)
        delete fPreamp;
    if (fChannel)
        delete fChannel;
    if (fCrystal_type)
        delete fCrystal_type;
    if (fApd_number)
        delete fApd_number;
    if (fVoltage)
        delete fVoltage;
    if (fFebex_pc)
        delete fFebex_pc;
    if (fFebex_slot)
        delete fFebex_slot;
    if (fFebex_mod)
        delete fFebex_mod;
    if (fFebex_channel)
        delete fFebex_channel;
    if (fLab)
        delete fLab;
    if (fMrcc_module)
        delete fMrcc_module;
    if (fMrcc_bus)
        delete fMrcc_bus;
    if (fMrcc_preamp)
        delete fMrcc_preamp;
    if (fIn_use)
        delete fIn_use;
}

// ----  Method clear ----------------------------------------------------------
void R3BCalifaMappingPar::clear()
{
    status = kFALSE;
    resetInputVersions();
}

// ----  Method putParams ------------------------------------------------------
void R3BCalifaMappingPar::putParams(FairParamList* list)
{
    LOG(info) << "R3BCalifaMappingPar::putParams() called";
    if (!list)
    {
        return;
    }
    fHalf->Set(fNumCrystals);
    list->add("califaHalfPar", *fHalf);
    fRing->Set(fNumCrystals);
    list->add("califaRingPar", *fRing);
    fPreamp->Set(fNumCrystals);
    list->add("califaPreampPar", *fPreamp);
    fChannel->Set(fNumCrystals);
    list->add("califaChannelPar", *fChannel);
    fCrystal_type->Set(fNumCrystals);
    list->add("califaCrystalTypePar", *fCrystal_type);
    fApd_number->Set(fNumCrystals);
    list->add("califaApdNumberPar", *fApd_number);
    fVoltage->Set(fNumCrystals);
    list->add("califaVoltagePar", *fVoltage);
    fFebex_pc->Set(fNumCrystals);
    list->add("califaFebexPCPar", *fFebex_pc);
    fFebex_slot->Set(fNumCrystals);
    list->add("califaFebexSlotPar", *fFebex_slot);
    fFebex_mod->Set(fNumCrystals);
    list->add("califaFebexModPar", *fFebex_mod);
    fFebex_channel->Set(fNumCrystals);
    list->add("califaFebexChannelPar", *fFebex_channel);
    fLab->Set(fNumCrystals);
    list->add("califaLabPar", *fLab);
    fMrcc_module->Set(fNumCrystals);
    list->add("califaMrccModulesPar", *fMrcc_module);
    fMrcc_bus->Set(fNumCrystals);
    list->add("califaMrccBusPar", *fMrcc_bus);
    fMrcc_preamp->Set(fNumCrystals);
    list->add("califaMrccPreampPar", *fMrcc_preamp);
    fIn_use->Set(fNumCrystals);
    list->add("califaInUsePar", *fIn_use);
    list->add("califaCrystalNumberPar", fNumCrystals);
}

// ----  Method getParams ------------------------------------------------------
Bool_t R3BCalifaMappingPar::getParams(FairParamList* list)
{
    LOG(info) << "R3BCalifaMappingPar::getParams() called";
    if (!list)
    {
        return kFALSE;
    }
    if (!list->fill("califaCrystalNumberPar", &fNumCrystals))
    {
        return kFALSE;
    }
    fHalf->Set(fNumCrystals);
    if (!(list->fill("califaHalfPar", fHalf)))
    {
        LOG(info) << "---Could not initialize califaHalfPar";
        return kFALSE;
    }
    fRing->Set(fNumCrystals);
    if (!(list->fill("califaRingPar", fRing)))
    {
        LOG(info) << "---Could not initialize califaRingPar";
        return kFALSE;
    }
    fPreamp->Set(fNumCrystals);
    if (!(list->fill("califaPreampPar", fPreamp)))
    {
        LOG(info) << "---Could not initialize califaPreampPar";
        return kFALSE;
    }
    fChannel->Set(fNumCrystals);
    if (!(list->fill("califaChannelPar", fChannel)))
    {
        LOG(info) << "---Could not initialize califaChannelPar";
        return kFALSE;
    }
    fCrystal_type->Set(fNumCrystals);
    if (!(list->fill("califaCrystalTypePar", fCrystal_type)))
    {
        LOG(info) << "---Could not initialize califaCrystalTypePar";
        return kFALSE;
    }
    fApd_number->Set(fNumCrystals);
    if (!(list->fill("califaApdNumberPar", fApd_number)))
    {
        LOG(info) << "---Could not initialize califaApdNumberPar";
        return kFALSE;
    }
    fVoltage->Set(fNumCrystals);
    if (!(list->fill("califaVoltagePar", fVoltage)))
    {
        LOG(info) << "---Could not initialize califaVoltagePar";
        return kFALSE;
    }
    fFebex_pc->Set(fNumCrystals);
    if (!(list->fill("califaFebexPCPar", fFebex_pc)))
    {
        LOG(info) << "---Could not initialize califaFebexPCPar";
        return kFALSE;
    }
    fFebex_slot->Set(fNumCrystals);
    if (!(list->fill("califaFebexSlotPar", fFebex_slot)))
    {
        LOG(info) << "---Could not initialize califaFebexSlotPar";
        return kFALSE;
    }
    fFebex_mod->Set(fNumCrystals);
    if (!(list->fill("califaFebexModPar", fFebex_mod)))
    {
        LOG(info) << "---Could not initialize califaFebexModPar";
        return kFALSE;
    }
    fFebex_channel->Set(fNumCrystals);
    if (!(list->fill("califaFebexChannelPar", fFebex_channel)))
    {
        LOG(info) << "---Could not initialize califaFebexChannelPar";
        return kFALSE;
    }
    fLab->Set(fNumCrystals);
    if (!(list->fill("califaLabPar", fLab)))
    {
        LOG(info) << "---Could not initialize califaLabPar";
        return kFALSE;
    }
    fMrcc_module->Set(fNumCrystals);
    if (!(list->fill("califaMrccModulesPar", fMrcc_module)))
    {
        LOG(info) << "---Could not initialize califaMrccModulesPar";
        return kFALSE;
    }
    fMrcc_bus->Set(fNumCrystals);
    if (!(list->fill("califaMrccBusPar", fMrcc_bus)))
    {
        LOG(info) << "---Could not initialize califaMrccBusPar";
        return kFALSE;
    }
    fMrcc_preamp->Set(fNumCrystals);
    if (!(list->fill("califaMrccPreampPar", fMrcc_preamp)))
    {
        LOG(info) << "---Could not initialize califaMrccPreampPar";
        return kFALSE;
    }
    fIn_use->Set(fNumCrystals);
    if (!(list->fill("califaInUsePar", fIn_use)))
    {
        LOG(info) << "---Could not initialize califaHacalifaInUseParlfPar";
        return kFALSE;
    }

    return kTRUE;
}

// ----  Method print ----------------------------------------------------------
void R3BCalifaMappingPar::print() { printParams(); }

// ----  Method printCrystalInfo -----------------------------------------------
void R3BCalifaMappingPar::printMapCrystalInfo(const UInt_t cryID)
{
    if (cryID < 1)
    {
        LOG(error) << "crystal_id must be given in 1-base";
        return;
    }
    else if (cryID > fNumCrystals)
    {
        LOG(error) << "crystal_id does not exist, crystal_id<=" << fNumCrystals;
        return;
    }
    auto index = cryID - 1;
    LOG(info) << "crystal_id " << cryID << " , "
              << "half " << fHalf->GetAt(index) << " , "
              << "ring " << fRing->GetAt(index) << " , "
              << "preamp " << fPreamp->GetAt(index) << " , "
              << "channel " << fChannel->GetAt(index) << " , "
              << "crystal_type " << fCrystal_type->GetAt(index) << " , "
              << "apd_number " << fApd_number->GetAt(index) << " , "
              << "voltage " << fVoltage->GetAt(index) << " , "
              << "febex_pc " << fFebex_pc->GetAt(index) << " , "
              << "febex_slot " << fFebex_slot->GetAt(index) << " , "
              << "febex_mod " << fFebex_mod->GetAt(index) << " , "
              << "febex_channel " << fFebex_channel->GetAt(index) << " , "
              << "lab " << fLab->GetAt(index) << " "
              << "mrcc_module " << fMrcc_module->GetAt(index) << " , "
              << "mrcc_bus " << fMrcc_bus->GetAt(index) << " , "
              << "mrcc_preamp " << fMrcc_preamp->GetAt(index) << " , "
              << "in_use " << fIn_use->GetAt(index);
}

// ----  Method printParams ----------------------------------------------------
void R3BCalifaMappingPar::printParams()
{
    LOG(info) << "R3BCalifaMappingPar: Califa Mapping Parameters:";

    LOG(info) << "crystal_id"
              << " "
              << "half"
              << " "
              << "ring"
              << " "
              << "preamp"
              << " "
              << "channel"
              << " "
              << "crystal_type"
              << " "
              << "apd_number"
              << " "
              << "voltage"
              << " "
              << "febex_pc"
              << " "
              << "febex_slot"
              << " "
              << "febex_mod"
              << " "
              << "febex_channel"
              << " "
              << "lab"
              << " "
              << "mrcc_module"
              << " "
              << "mrcc_bus"
              << " "
              << "mrcc_preamp"
              << " "
              << "in_use";
    for (Int_t i = 0; i < fNumCrystals; i++)
    {
        LOG(info) << i + 1 << " " << fHalf->GetAt(i) << " " << fRing->GetAt(i) << " " << fPreamp->GetAt(i) << " "
                  << fChannel->GetAt(i) << " " << fCrystal_type->GetAt(i) << " " << fApd_number->GetAt(i) << " "
                  << fVoltage->GetAt(i) << " " << fFebex_pc->GetAt(i) << " " << fFebex_slot->GetAt(i) << " "
                  << fFebex_mod->GetAt(i) << " " << fFebex_channel->GetAt(i) << " " << fLab->GetAt(i) << " "
                  << fMrcc_module->GetAt(i) << " " << fMrcc_bus->GetAt(i) << " " << fMrcc_preamp->GetAt(i) << " "
                  << fIn_use->GetAt(i);
    }
}

ClassImp(R3BCalifaMappingPar);
