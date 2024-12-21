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

#include <FairLogger.h>
#include <FairRootManager.h>
#include <FairRunAna.h>
#include <FairRuntimeDb.h>

#include "R3BCalifaCrystalCal2TotCalPar.h"
#include "R3BCalifaCrystalCalData.h"
#include "R3BCalifaMappingPar.h"
#include "R3BCalifaTotCalPar.h"

#include <TClonesArray.h>
#include <TF1.h>
#include <TGraph.h>
#include <TH1F.h>
#include <TMath.h>
#include <TObjArray.h>
#include <TRandom.h>
#include <TSpectrum.h>
#include <TVector3.h>
#include <iostream>
#include <stdlib.h>

R3BCalifaCrystalCal2TotCalPar::R3BCalifaCrystalCal2TotCalPar()
    : R3BCalifaCrystalCal2TotCalPar("R3BCalifaCrystalCal2TotCalPar", 1)
{
}

R3BCalifaCrystalCal2TotCalPar::R3BCalifaCrystalCal2TotCalPar(const char* name, Int_t iVerbose)
    : FairTask(name, iVerbose)
    , fMap_Par(NULL)
    , fTotCal_Par(NULL)
    , fCrystalCalDataCA(NULL)
    , fNumCrystals(1)
    , fNumParam(2)
    , fMinStadistics(100)
    , fThreshold(0)
    , fLeft(500.)
    , fRight(1400.)
    , fDebugMode(false)
{
}

R3BCalifaCrystalCal2TotCalPar::~R3BCalifaCrystalCal2TotCalPar()
{
    if (fCrystalCalDataCA)
        delete fCrystalCalDataCA;
}

void R3BCalifaCrystalCal2TotCalPar::SetParContainers()
{
    // Parameter Container
    // Reading califaMappingPar from FairRuntimeDb
    FairRuntimeDb* rtdb = FairRuntimeDb::instance();
    if (!rtdb)
    {
        LOG(error) << "FairRuntimeDb not opened!";
    }
    fMap_Par = dynamic_cast<R3BCalifaMappingPar*>(rtdb->getContainer("califaMappingPar"));
    if (!fMap_Par)
    {
        LOG(error) << "R3BCalifaCrystalCal2TotCalPar::Init() Couldn't get handle on califamappingpar container";
    }
    else
    {
        LOG(info) << "R3BCalifaCrystalCal2TotCalPar:: califamappingpar container open";
    }
}

void R3BCalifaCrystalCal2TotCalPar::SetParameter()
{
    if (!fMap_Par)
    {
        LOG(error) << "R3BCalifaMapped2CrystalCalPar::Container califaMappingPar not found.";
    }
    //--- Parameter Container ---
    fNumCrystals = fMap_Par->GetNumCrystals(); // Number of crystals x 2
    LOG(info) << "R3BCalifaMapped2CrystalCalPar::NumCry " << fNumCrystals;
    // fMap_Par->printParams();
}

InitStatus R3BCalifaCrystalCal2TotCalPar::Init()
{
    LOG(info) << "R3BCalifaCrystalCal2TotCalPar::Init()";

    FairRootManager* rootManager = FairRootManager::Instance();
    if (!rootManager)
    {
        LOG(error) << "R3BCalifaCrystalCal2TotCalPar::Init() FairRootManager not found";
        return kFATAL;
    }

    fCrystalCalDataCA = dynamic_cast<TClonesArray*>(rootManager->GetObject("CalifaCrystalCalData"));
    if (!fCrystalCalDataCA)
    {
        LOG(error) << "R3BCalifaCrystalCal2TotCalPar::Init() CalifaCrystalCalData not found";
        return kFATAL;
    }

    FairRuntimeDb* rtdb = FairRuntimeDb::instance();
    if (!rtdb)
    {
        LOG(error) << "R3BCalifaCrystalCal2TotCalPar::Init() FairRuntimeDb not found";
        return kFATAL;
    }

    fTotCal_Par = dynamic_cast<R3BCalifaTotCalPar*>(rtdb->getContainer("CalifaTotCalPar"));
    if (!fTotCal_Par)

    {
        LOG(error) << "R3BCalifaCrystalCal2TotCalPar::Init() Couldn't get handle on CalifaTotCalPar container";
        return kFATAL;
    }

    // Set container with mapping parameters
    SetParameter();

    // Create histograms for crystal calibration
    char name[100];
    energy_vs_tot_crystal = new TGraph*[fNumCrystals];
    LOG(debug) << "fNumCrystals = " << fNumCrystals;
    for (Int_t i = 0; i < fNumCrystals; i++)
        if (fMap_Par->GetInUse(i + 1) == 1)
        {
            sprintf(name, "energy_vs_tot_crystal_%d", i + 1);
            energy_vs_tot_crystal[i] = new TGraph();
            energy_vs_tot_crystal[i]->SetName(name);
        }

    return kSUCCESS;
}

InitStatus R3BCalifaCrystalCal2TotCalPar::ReInit()
{
    SetParContainers();
    SetParameter();
    return kSUCCESS;
}

void R3BCalifaCrystalCal2TotCalPar::Reset() {}

void R3BCalifaCrystalCal2TotCalPar::FinishEvent() {}

void R3BCalifaCrystalCal2TotCalPar::FinishTask()
{
    Search_TotParams();
    if (fDebugMode)
        for (Int_t i = 0; i < fNumCrystals; i++)
            if (fMap_Par->GetInUse(i + 1) == 1)
                energy_vs_tot_crystal[i]->Write();
}

void R3BCalifaCrystalCal2TotCalPar::Exec(Option_t* opt)
{
    Int_t nHits = fCrystalCalDataCA->GetEntries();
    if (!nHits)
        return;

    R3BCalifaCrystalCalData** CalHit = new R3BCalifaCrystalCalData*[nHits];
    Int_t crystalId = 0;

    for (Int_t i = 0; i < nHits; i++)
    {
        CalHit[i] = dynamic_cast<R3BCalifaCrystalCalData*>(fCrystalCalDataCA->At(i));
        crystalId = CalHit[i]->GetCrystalId();
        // Fill histograms
        if (CalHit[i]->GetEnergy() > fThreshold)
            energy_vs_tot_crystal[crystalId - 1]->SetPoint(
                energy_vs_tot_crystal[crystalId - 1]->GetN(), CalHit[i]->GetToTEnergy(), CalHit[i]->GetEnergy());
        LOG(debug) << "Tot = " << CalHit[i]->GetToTEnergy();
        LOG(debug) << "Energy = " << CalHit[i]->GetEnergy();
    }

    if (CalHit)
        delete[] CalHit;
    return;
}

void R3BCalifaCrystalCal2TotCalPar::Search_TotParams()
{

    fTotCal_Par->SetNumCrystals(fNumCrystals);
    fTotCal_Par->SetNumParametersFit(fNumParam);
    fTotCal_Par->GetCryCalParams()->Set(fNumParam * fNumCrystals);

    TF1* f1 = nullptr;
    if (fNumParam == 2 || fNumParam == 3)
    {
        if (fNumParam == 2)
        {
            f1 = new TF1("f1", "[0]*TMath::Exp(x/[1])", fLeft, fRight);
            f1->SetParameter(0, 7700);
            f1->SetParameter(1, 1000);
        }
        if (fNumParam == 3)
        {
            f1 = new TF1("f1", "[0]*TMath::Exp(x/[1]) + [2]", fLeft, fRight);
            f1->SetParameter(0, 7700);
            f1->SetParameter(1, 1000);
            f1->SetParameter(2, 100);
        }
    }
    else
    {
        LOG(info) << "R3BCalifaCrystalCal2TotCalPar::Search_TotParams() Wrong number of fit parameters, therefore, by "
                     "default "
                     "NumberParameters=2";
        f1 = new TF1("f1", "[0]*TMath::Exp(x/[1])", fLeft, fRight);
        f1->SetParameter(0, 7700);
        f1->SetParameter(1, 1000);
    }

    for (Int_t i = 0; i < fNumCrystals; i++)
        if (fMap_Par->GetInUse(i + 1) == 1)
        {
            energy_vs_tot_crystal[i]->Fit("f1", "R0");
            for (Int_t h = 0; h < fNumParam; h++)
            {
                fTotCal_Par->SetCryCalParams(f1->GetParameter(h), fNumParam * i + h + 1); // 1-base
            }
        }

        else
        {
            LOG(warn) << "R3BCalifaCrystalCal2TotCalPar::Histogram number " << i + 1 << " not Fitted";
        }

    fTotCal_Par->setChanged();
    return;
}

ClassImp(R3BCalifaCrystalCal2TotCalPar)
