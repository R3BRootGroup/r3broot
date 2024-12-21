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

#include "R3BNeulandMultiplicityCalorimetricPar.h"
#include "FairLogger.h"
#include "TObjString.h"

R3BNeulandMultiplicityCalorimetricPar::R3BNeulandMultiplicityCalorimetricPar(const char* name,
                                                                             const char* title,
                                                                             const char* context)
    : FairParGenericSet(name, title, context)
    , fNeutronCuts(nullptr)
{
}

R3BNeulandMultiplicityCalorimetricPar::~R3BNeulandMultiplicityCalorimetricPar()
{
    // Note: Deleting stuff here or in clear() causes segfaults?
}

void R3BNeulandMultiplicityCalorimetricPar::clear() {}

void R3BNeulandMultiplicityCalorimetricPar::putParams(FairParamList* l)
{
    if (!l)
    {
        return;
    }
    l->addObject("NeulandNeutronCuts", fNeutronCuts);
}

Bool_t R3BNeulandMultiplicityCalorimetricPar::getParams(FairParamList* l)
{
    if (!l)
    {
        return kFALSE;
    }
    if (!l->fillObject("NeulandNeutronCuts", fNeutronCuts))
    {
        return kFALSE;
    }
    return kTRUE;
}

void R3BNeulandMultiplicityCalorimetricPar::printParams()
{
    LOG(info) << "R3BNeulandMultiplicityCalorimetricPar: Neuland Neutron Cuts ...";
    for (const auto& nc : GetNeutronCuts())
    {
        LOG(info) << nc.first;
        nc.second->Print();
    }
}

void R3BNeulandMultiplicityCalorimetricPar::SetNeutronCuts(const std::map<UInt_t, TCutG*>& cuts)
{
    fNeutronCuts = new TMap();
    fNeutronCuts->SetOwner(kTRUE);
    for (const auto& nc : cuts)
    {
        auto key = new TObjString(TString::Itoa(nc.first, 10));
        fNeutronCuts->Add(key, nc.second->Clone());
    }
}

std::map<UInt_t, TCutG*> R3BNeulandMultiplicityCalorimetricPar::GetNeutronCuts() const
{
    if (fNeutronCuts == nullptr)
    {
        LOG(fatal) << "R3BNeulandMultiplicityCalorimetricPar: NeutronCuts not set!";
    }

    std::map<UInt_t, TCutG*> map;

    TObjString* key;
    TIterator* nextobj = fNeutronCuts->MakeIterator();
    while ((key = dynamic_cast<TObjString*>(nextobj->Next())))
    {
        UInt_t nNeutrons = key->GetString().Atoi();
        map[nNeutrons] = dynamic_cast<TCutG*>(fNeutronCuts->GetValue(key)->Clone());
    }

    return map;
}

TCutG* R3BNeulandMultiplicityCalorimetricPar::GetNeutronCut(const Int_t n) const { return GetNeutronCuts().at(n); }

UInt_t R3BNeulandMultiplicityCalorimetricPar::GetNeutronMultiplicity(const Double_t energy,
                                                                     const Double_t nClusters) const
{
    // Note: it might be better to implement std::map as a member and sync between the tmap for this type of usage.
    if (fNeutronCuts == nullptr)
    {
        LOG(fatal) << "R3BNeulandMultiplicityCalorimetricPar: NeutronCuts not set!";
    }

    if (nClusters < 1)
    {
        return 0;
    }

    TObjString* key;
    TIterator* nextobj = fNeutronCuts->MakeIterator();
    while ((key = dynamic_cast<TObjString*>(nextobj->Next())))
    {
        // UInt_t nNeutrons = key->GetString().Atoi();
        if ((dynamic_cast<TCutG*>(fNeutronCuts->GetValue(key)))->IsInside(energy, nClusters))
        {
            return (UInt_t)key->GetString().Atoi();
        }
    }
    // Assume if no match is found, the neutron multiplicity must be higher than the highest saved cut.
    return GetNeutronCuts().rbegin()->first + 1;
}

ClassImp(R3BNeulandMultiplicityCalorimetricPar);
