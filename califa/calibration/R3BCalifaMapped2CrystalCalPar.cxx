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

#include "TClonesArray.h"
#include "TF1.h"
#include "TGraph.h"
#include "TH1F.h"
#include "TMath.h"
#include "TObjArray.h"
#include "TRandom.h"
#include "TSpectrum.h"
#include "TVector3.h"

#include "FairLogger.h"
#include "FairRootManager.h"
#include "FairRunAna.h"
#include "FairRuntimeDb.h"

#include "R3BCalifaCrystalCalPar.h"
#include "R3BCalifaMapped2CrystalCalPar.h"
#include "R3BCalifaMappedData.h"
#include "R3BCalifaMappingPar.h"

#include <iostream>
#include <stdlib.h>

R3BCalifaMapped2CrystalCalPar::R3BCalifaMapped2CrystalCalPar()
    : R3BCalifaMapped2CrystalCalPar("R3B CALIFA Calibration Parameters Finder ", 1)
{
}

R3BCalifaMapped2CrystalCalPar::R3BCalifaMapped2CrystalCalPar(const char* name, Int_t iVerbose)
    : FairTask(name, iVerbose)
    , fMap_Par(NULL)
    , fCal_Par(NULL)
    , fCalifaMappedDataCA(NULL)
    , fNumCrystals(1)
    , fNumParam(0)
    , fMinStadistics(100)
    , fMapHistos_left(0)
    , fMapHistos_right(0)
    , fMapHistos_bins(0)
    , fMapHistos_leftp(0)
    , fMapHistos_rightp(0)
    , fMapHistos_binsp(0)
    , fNumPeaks(0)
    , fSigma(0)
    , fThreshold(0)
    , fEnergyPeaks(NULL)
    , fDebugMode(0)
{
}

R3BCalifaMapped2CrystalCalPar::~R3BCalifaMapped2CrystalCalPar()
{
    LOG(info) << "R3BCalifaMapped2CrystalCalPar: Delete instance";
    if (fCalifaMappedDataCA)
        delete fCalifaMappedDataCA;
    if (fEnergyPeaks)
        delete fEnergyPeaks;
}

void R3BCalifaMapped2CrystalCalPar::SetParContainers()
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
        LOG(error) << "R3BCalifaMapped2CrystalCalPar::Init() Couldn't get handle on califaMappingPar container";
    }
    else
    {
        LOG(info) << "R3BCalifaMapped2CrystalCalPar:: califaMappingPar container open";
    }
}

void R3BCalifaMapped2CrystalCalPar::SetParameter()
{
    if (!fMap_Par)
    {
        LOG(warn) << "R3BCalifaMapped2CrystalCalPar::Container califaMappingPar not found.";
    }
    //--- Parameter Container ---
    fNumCrystals = fMap_Par->GetNumCrystals(); // Number of crystals x 2
    LOG(info) << "R3BCalifaMapped2CrystalCalPar::NumCry " << fNumCrystals;
    // fMap_Par->printParams();
}

InitStatus R3BCalifaMapped2CrystalCalPar::Init()
{
    LOG(info) << "R3BCalifaMapped2CrystalCalPar::Init()";

    if (!fEnergyPeaks)
    {
        fEnergyPeaks = new TArrayF;
        fEnergyPeaks->Set(fNumPeaks);
    }

    FairRootManager* rootManager = FairRootManager::Instance();
    if (!rootManager)
    {
        LOG(error) << "R3BCalifaMapped2CrystalCalPar::Init() FairRootManager not found";
        return kFATAL;
    }

    fCalifaMappedDataCA = dynamic_cast<TClonesArray*>(rootManager->GetObject("CalifaMappedData"));
    if (!fCalifaMappedDataCA)
    {
        LOG(error) << "R3BCalifaMapped2CrystalCalPar::Init() CalifaMappedData not found";
        return kFATAL;
    }

    FairRuntimeDb* rtdb = FairRuntimeDb::instance();
    if (!rtdb)
    {
        LOG(error) << "R3BCalifaMapped2CrystalCalPar::Init() FairRuntimeDb not found";
        return kFATAL;
    }

    fCal_Par = dynamic_cast<R3BCalifaCrystalCalPar*>(rtdb->getContainer("califaCrystalCalPar"));
    if (!fCal_Par)
    {
        LOG(error) << "R3BCalifaMapped2CrystalCalPar::Init() Couldn't get handle on califaCrystalCalPar container";
        return kFATAL;
    }

    // Set container with mapping parameters
    SetParameter();

    // Create histograms for crystal calibration
    char name[100];
    Int_t fright, fleft, fbins;
    fh_Map_energy_crystal = new TH1F*[fNumCrystals];
    for (Int_t i = 0; i < fNumCrystals; i++)
        if (fMap_Par->GetInUse(i + 1) == 1)
        {

            sprintf(name, "fh_Map_energy_crystal_%i", i + 1);
            if (i < fMap_Par->GetNumCrystals() / 2)
            {
                fright = fMapHistos_right;
                fleft = fMapHistos_left;
                fbins = fMapHistos_bins;
            }
            else
            {
                fright = fMapHistos_rightp;
                fleft = fMapHistos_leftp;
                fbins = fMapHistos_binsp;
            }
            fh_Map_energy_crystal[i] = new TH1F(name, name, fbins, fleft, fright);
        }

    return kSUCCESS;
}

InitStatus R3BCalifaMapped2CrystalCalPar::ReInit()
{
    SetParContainers();
    SetParameter();
    return kSUCCESS;
}

void R3BCalifaMapped2CrystalCalPar::Exec(Option_t* opt)
{
    Int_t nHits = fCalifaMappedDataCA->GetEntries();
    if (!nHits)
        return;

    R3BCalifaMappedData** MapHit = new R3BCalifaMappedData*[nHits];
    Int_t crystalId = 0;

    for (Int_t i = 0; i < nHits; i++)
    {
        MapHit[i] = dynamic_cast<R3BCalifaMappedData*>(fCalifaMappedDataCA->At(i));
        crystalId = MapHit[i]->GetCrystalId();
        // Fill histograms
        if (fMap_Par->GetInUse(crystalId) == 1)
            fh_Map_energy_crystal[crystalId - 1]->Fill(MapHit[i]->GetEnergy());
    }

    if (MapHit)
        delete[] MapHit;
    return;
}

void R3BCalifaMapped2CrystalCalPar::Reset() {}

void R3BCalifaMapped2CrystalCalPar::FinishEvent() {}

void R3BCalifaMapped2CrystalCalPar::FinishTask()
{
    // Look for the peaks in each spectrum
    SearchPeaks();
    // fCal_Par->printParams();
    if (fDebugMode)
        for (Int_t i = 0; i < fNumCrystals; i++)
            if (fMap_Par->GetInUse(i + 1) == 1 && fh_Map_energy_crystal[i]->GetEntries() > fMinStadistics)
                fh_Map_energy_crystal[i]->Write();
}

void R3BCalifaMapped2CrystalCalPar::SearchPeaks()
{
    Int_t nfound = 0;
    Int_t numPars = 2; // Number of parameters=2 by default
    if (fNumParam)
    {
        numPars = fNumParam;
    }

    fCal_Par->SetNumCrystals(fNumCrystals);
    fCal_Par->SetNumParametersFit(fNumParam);
    fCal_Par->GetCryCalParams()->Set(numPars * fNumCrystals);

    TSpectrum* ss = new TSpectrum(fNumPeaks);

    Int_t fright, fleft;

    for (Int_t i = 0; i < fNumCrystals; i++)
        if (fMap_Par->GetInUse(i + 1) == 1)
        {

            if (fh_Map_energy_crystal[i]->GetEntries() > fMinStadistics)
            {

                if (fDebugMode)
                    nfound = ss->Search(fh_Map_energy_crystal[i], fSigma, "", fThreshold); // number of peaks
                else
                    nfound = ss->Search(fh_Map_energy_crystal[i], fSigma, "goff", fThreshold);
                LOG(debug) << "CrystalId=" << i + 1 << " " << nfound << " " << fThreshold;
                fChannelPeaks = (Double_t*)ss->GetPositionX();

                Int_t idx[nfound];
                TMath::Sort(nfound, fChannelPeaks, idx, kTRUE);

                // Calibrated Spectrum
                Double_t X[nfound + 1];
                Double_t Y[nfound + 1];

                for (Int_t j = 0; j < nfound; j++)
                {
                    X[j] = fChannelPeaks[idx[nfound - j - 1]];
                    Y[j] = fEnergyPeaks->GetAt(nfound - j - 1);
                    LOG(debug) << "CrystalId=" << i + 1 << " " << j + 1 << " " << X[j + 1];
                }
                X[nfound] = 0.;
                Y[nfound] = 0.;

                if (i < fMap_Par->GetNumCrystals() / 2)
                {
                    fright = fMapHistos_right;
                    fleft = fMapHistos_left;
                }
                else
                {
                    fright = fMapHistos_rightp;
                    fleft = fMapHistos_leftp;
                }

                TF1* f1 = nullptr;
                if (fNumParam)
                {

                    if (fNumParam == 1)
                    {
                        f1 = new TF1("f1", "[0]*x", fleft, fright);
                    }
                    if (fNumParam == 2)
                    {
                        f1 = new TF1("f1", "[0]+[1]*x", fleft, fright);
                    }
                    if (fNumParam == 3)
                    {
                        f1 = new TF1("f1", "[0]+[1]*x+[2]*pow(x,2)", fleft, fright);
                    }
                    if (fNumParam == 4)
                    {
                        f1 = new TF1("f1", "[0]+[1]*x+[2]*pow(x,2)+[3]*pow(x,3)", fleft, fright);
                    }
                    if (fNumParam == 5)
                    {
                        f1 = new TF1("f1", "[0]+[1]*x+[2]*pow(x,2)+[3]*pow(x,3)+[4]*pow(x,4)", fleft, fright);
                    }
                    if (fNumParam > 5)
                    {
                        LOG(warn)
                            << "R3BCalifaMapped2CrystalCalPar:: The number of fit parameters can not be higher than 5";
                    }
                }
                else
                {
                    LOG(warn)
                        << "R3BCalifaMapped2CrystalCalPar:: No imput number of fit parameters, therefore, by default "
                           "NumberParameters=2";
                    f1 = new TF1("f1", "[0]+[1]*x", fleft, fright);
                }

                TGraph* graph = new TGraph(fNumPeaks + 1, X, Y);
                graph->Fit("f1", "Q"); // Quiet mode (minimum printing)

                for (Int_t h = 0; h < numPars; h++)
                {
                    fCal_Par->SetCryCalParams(f1->GetParameter(h), numPars * i + h + 1); // 1-base
                }
            }
            else
            {
                LOG(warn) << "R3BCalifaMapped2CrystalCalPar::Histogram number " << i + 1 << "not Fitted";
            }
        }

    delete ss;
    fCal_Par->setChanged();
    return;
}

ClassImp(R3BCalifaMapped2CrystalCalPar)
