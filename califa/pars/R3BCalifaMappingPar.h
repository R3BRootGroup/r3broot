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

#include <TArrayF.h>
#include <TArrayI.h>
#include <TObjArray.h>
#include <TObjString.h>
#include <TObject.h>

class FairParamList;

class R3BCalifaMappingPar : public FairParGenericSet
{
  public:
    /** Standard constructor **/
    R3BCalifaMappingPar(const char* name = "califaMappingPar",
                        const char* title = "Califa Mapping Parameters",
                        const char* context = "CalifaMappingContext");

    /** Destructor **/
    virtual ~R3BCalifaMappingPar();

    /** Reset all parameters **/
    virtual void clear();

    /** Store all parameters using FairRuntimeDB **/
    virtual void putParams(FairParamList* list);

    /** Retrieve all parameters using FairRuntimeDB**/
    Bool_t getParams(FairParamList* list);

    /** Print values of parameters to the standard output **/
    virtual void print();
    void printParams();
    void printMapCrystalInfo(const UInt_t cryID);

    /** Accessor functions **/
    [[nodiscard]] inline const int GetNumCrystals() const { return fNumCrystals; }
    const Int_t GetHalf(Int_t crystal) { return fHalf->GetAt(crystal - 1); }
    const Int_t GetRing(Int_t crystal) { return fRing->GetAt(crystal - 1); }
    const Int_t GetPreamp(Int_t crystal) { return fPreamp->GetAt(crystal - 1); }
    const Int_t GetChannel(Int_t crystal) { return fChannel->GetAt(crystal - 1); }
    const Int_t GetCrystalType(Int_t crystal) { return fCrystal_type->GetAt(crystal - 1); }
    const Int_t GetApdNumber(Int_t crystal) { return fApd_number->GetAt(crystal - 1); }
    const Float_t GetVoltage(Int_t crystal) { return fVoltage->GetAt(crystal - 1); }
    const Int_t GetFebexPC(Int_t crystal) { return fFebex_pc->GetAt(crystal - 1); }
    const Int_t GetFebexSlot(Int_t crystal) { return fFebex_slot->GetAt(crystal - 1); }
    const Int_t GetFebexMod(Int_t crystal) { return fFebex_mod->GetAt(crystal - 1); }
    const Int_t GetFebexChannel(Int_t crystal) { return fFebex_channel->GetAt(crystal - 1); }
    const Int_t GetLab(Int_t crystal) { return fLab->GetAt(crystal - 1); }
    const Int_t GetMrccModule(Int_t crystal) { return fMrcc_module->GetAt(crystal - 1); }
    const Int_t GetMrccBus(Int_t crystal) { return fMrcc_bus->GetAt(crystal - 1); }
    const Int_t GetMrccPreamp(Int_t crystal) { return fMrcc_preamp->GetAt(crystal - 1); }
    const Int_t GetInUse(Int_t crystal) { return fIn_use->GetAt(crystal - 1); }

    inline void SetNumCrystals(int numberCry) { fNumCrystals = numberCry; }
    inline void SetHalf(Int_t value, Int_t crystal) { fHalf->AddAt(value, crystal - 1); }
    inline void SetRing(Int_t value, Int_t crystal) { fRing->AddAt(value, crystal - 1); }
    inline void SetPreamp(Int_t value, Int_t crystal) { fPreamp->AddAt(value, crystal - 1); }
    inline void SetChannel(Int_t value, Int_t crystal) { fChannel->AddAt(value, crystal - 1); }
    inline void SetCrystalType(Int_t value, Int_t crystal) { fCrystal_type->AddAt(value, crystal - 1); }
    inline void SetApdNumber(Int_t value, Int_t crystal) { fApd_number->AddAt(value, crystal - 1); }
    inline void SetVoltage(Float_t value, Int_t crystal) { fVoltage->AddAt(value, crystal - 1); }
    inline void SetFebexPC(Int_t value, Int_t crystal) { fFebex_pc->AddAt(value, crystal - 1); }
    inline void SetFebexSlot(Int_t value, Int_t crystal) { fFebex_slot->AddAt(value, crystal - 1); }
    inline void SetFebexMod(Int_t value, Int_t crystal) { fFebex_mod->AddAt(value, crystal - 1); }
    inline void SetFebexChannel(Int_t value, Int_t crystal) { fFebex_channel->AddAt(value, crystal - 1); }
    inline void SetLab(Int_t value, Int_t crystal) { fLab->AddAt(value, crystal - 1); }
    inline void SetMrccModule(Int_t value, Int_t crystal) { fMrcc_module->AddAt(value, crystal - 1); }
    inline void SetMrccBus(Int_t value, Int_t crystal) { fMrcc_bus->AddAt(value, crystal - 1); }
    inline void SetMrccPreamp(Int_t value, Int_t crystal) { fMrcc_preamp->AddAt(value, crystal - 1); }
    inline void SetInUse(Int_t value, Int_t crystal) { fIn_use->AddAt(value, crystal - 1); }

  private:
    int fNumCrystals = 5088; // number of crystals
    TArrayI* fHalf;          // half (1 Wixhausen side or RIGHT side, 2 Messel side or LEFT side)
    TArrayI* fRing;          // ring from 1 to 5 (5 is iPhos)
    TArrayI* fPreamp;        // preamp from 1 to 16
    TArrayI* fChannel;       // preamp channel from 1 to 16
    TArrayI* fCrystal_type;  // crystal type (crystal engineering model)
    TArrayI* fApd_number;    // APD number attached to the crystal
    TArrayF* fVoltage;       // nominal voltage for gain=50 (volts)
    TArrayI* fFebex_pc;      // FEBEX pc from 0 to 1
    TArrayI* fFebex_slot;    // FEBEX spf from 0 to 3
    TArrayI* fFebex_mod;     // FEBEX module from 0 to 15 (0 to f)
    TArrayI* fFebex_channel; // FEBEX channel from 0 to 15
    TArrayI* fLab;           // 1:LUND 2:TUD 3:USC
    TArrayI* fMrcc_module;   // Slow Control MRCC module 0 to 1
    TArrayI* fMrcc_bus;      // Slow Control MRCC bus  (0 or 1)
    TArrayI* fMrcc_preamp;   // Slow Control MRCC bus from 0 to 15 (0 to f)
    TArrayI* fIn_use;        // 1: crystal installed and ready  0:otherwise

    const R3BCalifaMappingPar& operator=(const R3BCalifaMappingPar&); /*< an assignment operator>*/
    R3BCalifaMappingPar(const R3BCalifaMappingPar&);                  /*< a copy constructor >*/

  public:
    ClassDef(R3BCalifaMappingPar, 1);
};
