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

#include <FairRootManager.h>

#include "R3BCalifa.h"
#include "R3BCalifaCrystalCalData.h"
#include "R3BCalifaDigitizer.h"
#include "R3BCalifaGeometry.h"
#include "R3BCalifaPoint.h"
#include "R3BLogger.h"

#include <TArrayD.h>
#include <TClonesArray.h>
#include <TMath.h>
#include <TRandom.h>
#include <TVector3.h>
#include <cmath>
#include <iostream>
#include <stdlib.h>

auto MakeCompFun = [](std::array<double, 6> par)
{
    return [par](double val)
    { return (val > 0.0) ? par[5] / (par[0] + par[1] * pow(val, par[2]) + par[3] / pow(val, par[4])) : 0.0; };
};

R3BCalifaDigitizer::R3BCalifaDigitizer()
    : FairTask("R3BCalifaDigitizer")
{
}

R3BCalifaDigitizer::~R3BCalifaDigitizer()
{
    if (fCalifaCryCalDataCA)
    {
        fCalifaCryCalDataCA->Delete();
        delete fCalifaCryCalDataCA;
    }
}

void R3BCalifaDigitizer::SetParContainers()
{
    if (fRealConfig)
    {
        FairRuntimeDb* rtdb = FairRuntimeDb::instance();
        R3BLOG_IF(error, !rtdb, "FairRuntimeDb not found");

        fSim_Par = dynamic_cast<R3BCalifaCrystalPars4Sim*>(rtdb->getContainer("califaCrystalPars4Sim"));
        if (!fSim_Par)
        {
            R3BLOG(error, "Couldn't get handle on califaCrystalPars4Sim container");
        }
        else
        {
            R3BLOG(info, "califaCrystalPars4Sim container opened");
        }
    }
}

void R3BCalifaDigitizer::SetParameter()
{
    if (fRealConfig == true)
    {
        //--- Parameter Container ---
        fNumCrystals = fSim_Par->GetNumCrystals() / 2;      // Number of Crystals
        fNumberOfParams = fSim_Par->GetNumParameters4Sim(); // Number of Parameters

        const auto numcryroot = R3BCalifaGeometry::Instance()->GetNbCrystals() / 2;
        R3BLOG_IF(fatal,
                  fNumCrystals != numcryroot,
                  "The number of crystals in the container califaCrystalPars4Sim: "
                      << fNumCrystals << " and the Root geometry: " << numcryroot << " are different");

        fSim_Par->printParams();

        R3BLOG(info, "Max number of crystals: " << fNumCrystals);
        R3BLOG(info, "Nb of parameters used in simulation: " << fNumberOfParams);
    }
    else
    {
        // R3BCalifaGeometry returns the total number of crystals in proton and gamma range
        fNumCrystals = R3BCalifaGeometry::Instance()->GetNbCrystals() / 2;
        R3BLOG(info, "Max number of crystals: " << fNumCrystals);
    }
}

InitStatus R3BCalifaDigitizer::Init()
{
    R3BLOG(info, "");
    // Get a handle from the IO manager
    auto* rootman = FairRootManager::Instance();
    if (rootman == nullptr)
    {
        R3BLOG(error, "FairRootManager not found");
        return kFATAL;
    }

    fCalifaPointDataCA = dynamic_cast<TClonesArray*>(rootman->GetObject("CalifaPoint"));
    if (!fCalifaPointDataCA)
    {
        R3BLOG(error, "CalifaPoint not found");
        return kFATAL;
    }

    fCalifaCryCalDataCA = new TClonesArray("R3BCalifaCrystalCalData");
    rootman->Register("CalifaCrystalCalData", "CALIFA Crystal Cal", fCalifaCryCalDataCA, true);
    // output_data.init();

    SetParameter();
    return kSUCCESS;
}

void R3BCalifaDigitizer::Exec(Option_t*)
{
    // Reset entries for output
    // output_data.clear();
    if (fCalifaCryCalDataCA)
        fCalifaCryCalDataCA->Clear();

    // Reading the Input -- Point data --
    const int nHits = fCalifaPointDataCA->GetEntriesFast();
    if (!nHits)
        return;

    auto tf_dNf_dE = MakeCompFun({ -1.79, 1.36e-2, 7.84e-1, 4.97, 1.75e-1, 0.78 });
    auto tf_dNs_dE = MakeCompFun({ -1.24e2, 6.3e-3, 1.27, 1.262e2, 2.3e-3, 1.32 });

    R3BCalifaPoint** pointData = nullptr;
    pointData = new R3BCalifaPoint*[nHits];
    for (int i = 0; i < nHits; i++)
        pointData[i] = dynamic_cast<R3BCalifaPoint*>(fCalifaPointDataCA->At(i));

    int crystalId = 0;
    double Nf = 0.;
    double Ns = 0.;
    uint64_t time = 0;
    double energy = 0.;
    double pathlength = 0.;

    for (int i = 0; i < nHits; i++)
    {
        crystalId = pointData[i]->GetCrystalId();
        time = pointData[i]->GetTime();
        energy = NUSmearing(pointData[i]->GetEnergyLoss());
        pathlength = pointData[i]->GetLength();
        auto dE_dx = energy / pathlength;
        if (fComponentRes > 0)
        {
            auto dE_dx_ns = CompSmearing(dE_dx);
            auto dE_dx_nf = CompSmearing(dE_dx);
            auto eneres = ExpResSmearing(energy);
            Nf = tf_dNf_dE(dE_dx_nf) * eneres;
            Ns = tf_dNs_dE(dE_dx_ns) * eneres;
        }
        else
        {
            Nf = tf_dNf_dE(dE_dx) * energy;
            Ns = tf_dNs_dE(dE_dx) * energy;
        }

        int nCrystalCals = 0;
        if (fCalifaCryCalDataCA != nullptr)
            nCrystalCals = fCalifaCryCalDataCA->GetEntriesFast();
        bool existHit = false;
        if (nCrystalCals == 0)
            AddCrystalCal(crystalId, energy, Nf, Ns, time, energy);
        else
        {
            for (int j = 0; j < nCrystalCals; j++)
            {
                if ((dynamic_cast<R3BCalifaCrystalCalData*>(fCalifaCryCalDataCA->At(j)))->GetCrystalId() == crystalId)
                {
                    (dynamic_cast<R3BCalifaCrystalCalData*>(fCalifaCryCalDataCA->At(j)))->AddMoreEnergy(energy);
                    (dynamic_cast<R3BCalifaCrystalCalData*>(fCalifaCryCalDataCA->At(j)))->AddMoreTot(energy);
                    (dynamic_cast<R3BCalifaCrystalCalData*>(fCalifaCryCalDataCA->At(j)))->AddMoreNf(Nf);
                    (dynamic_cast<R3BCalifaCrystalCalData*>(fCalifaCryCalDataCA->At(j)))->AddMoreNs(Ns);
                    if ((dynamic_cast<R3BCalifaCrystalCalData*>(fCalifaCryCalDataCA->At(j)))->GetTime() > time)
                    {
                        (dynamic_cast<R3BCalifaCrystalCalData*>(fCalifaCryCalDataCA->At(j)))->SetTime(time);
                    }
                    existHit = true; // to avoid the creation of a new CrystalHit
                    break;
                }
            }
            if (!existHit)
                AddCrystalCal(crystalId, energy, Nf, Ns, time, energy);
        }
    }

    if (pointData)
        delete[] pointData;

    int nCrystalCals = 0;
    if (fCalifaCryCalDataCA != nullptr)
        nCrystalCals = fCalifaCryCalDataCA->GetEntriesFast();

    if (nCrystalCals == 0)
    {
        return;
    }

    if (fRealConfig == true)
    {
        FillRealConfig(nCrystalCals);
    }
    else if (fFullProtonRange == true)
    {
        FillIdealConfigWithProtonRange(nCrystalCals);
    }
    else
    {
        FillIdealConfig(nCrystalCals);
    }

    /*
            int nCrysCals = fCalifaCryCalDataCA->GetEntriesFast();
            for (int i = 0; i < nCrysCals; i++)
            {

                auto tempCal = dynamic_cast<R3BCalifaCrystalCalData*>(fCalifaCryCalDataCA->At(i));

                R3BCalifaCrystalCalData hit(tempCal->GetCrystalId(),
                                            tempCal->GetEnergy(),
                                            tempCal->GetNf(),
                                            tempCal->GetNs(),
                                            tempCal->GetTime(),
                                            tempCal->GetToTEnergy());
                output_data.get().push_back(hit);
            }*/
}

void R3BCalifaDigitizer::FillRealConfig(int nCrystalCals)
{
    double tempE = 0;
    double tempNf = 0;
    double tempNs = 0;
    int tempCryID = 0;
    int parThres = 0;
    bool inUse = false;

    for (int i = 0; i < nCrystalCals; i++)
    {
        tempCryID = (dynamic_cast<R3BCalifaCrystalCalData*>(fCalifaCryCalDataCA->At(i)))->GetCrystalId();
        tempE = (dynamic_cast<R3BCalifaCrystalCalData*>(fCalifaCryCalDataCA->At(i)))->GetEnergy();

        inUse = fSim_Par->GetInUse(tempCryID - 1);
        fResolution = fSim_Par->GetResolution(tempCryID - 1);
        parThres = fSim_Par->GetThreshold(tempCryID - 1);

        if (inUse && parThres < tempE * 1000.)
        { // Thresholds are in KeV!!

            (dynamic_cast<R3BCalifaCrystalCalData*>(fCalifaCryCalDataCA->At(i)))->SetEnergy(ExpResSmearing(tempE));
            if (fComponentRes > 0)
            {
                tempNf = (dynamic_cast<R3BCalifaCrystalCalData*>(fCalifaCryCalDataCA->At(i)))->GetNf();
                (dynamic_cast<R3BCalifaCrystalCalData*>(fCalifaCryCalDataCA->At(i)))->SetNf(CompSmearing(tempNf));
                tempNs = (dynamic_cast<R3BCalifaCrystalCalData*>(fCalifaCryCalDataCA->At(i)))->GetNs();
                (dynamic_cast<R3BCalifaCrystalCalData*>(fCalifaCryCalDataCA->At(i)))->SetNs(CompSmearing(tempNs));
            }
        }
        else
        {
            fCalifaCryCalDataCA->RemoveAt(i); // remove from CalData those below threshold
            fCalifaCryCalDataCA->Compress();
            nCrystalCals--;
            i--;
            continue;
        }
    }
}

void R3BCalifaDigitizer::FillIdealConfig(int nCrystalCals)
{
    double tempE = 0;
    double realE = 0;
    double tempNf = 0;
    double realNf = 0;
    double tempNs = 0;
    double realNs = 0;
    int tempCryID = 0;

    for (int i = 0; i < nCrystalCals; i++)
    {
        tempCryID = (dynamic_cast<R3BCalifaCrystalCalData*>(fCalifaCryCalDataCA->At(i)))->GetCrystalId();

        tempE = (dynamic_cast<R3BCalifaCrystalCalData*>(fCalifaCryCalDataCA->At(i)))->GetEnergy();
        tempNf = (dynamic_cast<R3BCalifaCrystalCalData*>(fCalifaCryCalDataCA->At(i)))->GetNf();
        tempNs = (dynamic_cast<R3BCalifaCrystalCalData*>(fCalifaCryCalDataCA->At(i)))->GetNs();

        realE = tempE;
        realNf = tempNf;
        realNs = tempNs;

        if (tempE < fThreshold)
        {
            fCalifaCryCalDataCA->RemoveAt(i);
            fCalifaCryCalDataCA->Compress();
            nCrystalCals--; // remove from CalData those below threshold
            i--;
            continue;
        }

        if (tempCryID <= fNumCrystals)
        {
            // Gamma branch
            if (fGammaResolution > 0)
            {
                fResolution = fGammaResolution;
                realE = ExpResSmearing(tempE);
            }

            if (tempE >= fGammaSaturation && fGammaSaturation > 0 && fGammaResolution > 0)
            {
                (dynamic_cast<R3BCalifaCrystalCalData*>(fCalifaCryCalDataCA->At(i)))
                    ->SetEnergy(ExpResSmearing(fGammaSaturation));
                (dynamic_cast<R3BCalifaCrystalCalData*>(fCalifaCryCalDataCA->At(i)))
                    ->SetNf(ExpResSmearing(fGammaSaturation));
                (dynamic_cast<R3BCalifaCrystalCalData*>(fCalifaCryCalDataCA->At(i)))
                    ->SetNs(ExpResSmearing(fGammaSaturation));
                (dynamic_cast<R3BCalifaCrystalCalData*>(fCalifaCryCalDataCA->At(i)))->SetToTEnergy(realE);
            }
            else if (tempE >= fGammaSaturation && fGammaSaturation > 0 && fGammaResolution == 0)
            {
                (dynamic_cast<R3BCalifaCrystalCalData*>(fCalifaCryCalDataCA->At(i)))->SetEnergy(fGammaSaturation);
                (dynamic_cast<R3BCalifaCrystalCalData*>(fCalifaCryCalDataCA->At(i)))->SetNf(fGammaSaturation);
                (dynamic_cast<R3BCalifaCrystalCalData*>(fCalifaCryCalDataCA->At(i)))->SetNs(fGammaSaturation);
                (dynamic_cast<R3BCalifaCrystalCalData*>(fCalifaCryCalDataCA->At(i)))->SetToTEnergy(realE);
            }
            else if (fGammaResolution > 0)
            {
                (dynamic_cast<R3BCalifaCrystalCalData*>(fCalifaCryCalDataCA->At(i)))->SetEnergy(realE);
                (dynamic_cast<R3BCalifaCrystalCalData*>(fCalifaCryCalDataCA->At(i)))->SetNf(realNf);
                (dynamic_cast<R3BCalifaCrystalCalData*>(fCalifaCryCalDataCA->At(i)))->SetNs(realNs);
                (dynamic_cast<R3BCalifaCrystalCalData*>(fCalifaCryCalDataCA->At(i)))->SetToTEnergy(realE);
            }
            else if (fGammaResolution == 0)
            {
                (dynamic_cast<R3BCalifaCrystalCalData*>(fCalifaCryCalDataCA->At(i)))->SetEnergy(realE);
                (dynamic_cast<R3BCalifaCrystalCalData*>(fCalifaCryCalDataCA->At(i)))->SetNf(realNf);
                (dynamic_cast<R3BCalifaCrystalCalData*>(fCalifaCryCalDataCA->At(i)))->SetNs(realNs);
                (dynamic_cast<R3BCalifaCrystalCalData*>(fCalifaCryCalDataCA->At(i)))->SetToTEnergy(realE);
            }

            // Proton branch
            if (tempCryID > barrelCrystals)
            {
                if (fProtonResolution > 0)
                {
                    fResolution = fProtonResolution;
                    realE = ExpResSmearing(tempE);
                }
                auto time = (dynamic_cast<R3BCalifaCrystalCalData*>(fCalifaCryCalDataCA->At(i)))->GetTime();

                AddCrystalCal(tempCryID + fNumCrystals,
                              ((realE >= fProtonSaturation && fProtonSaturation > 0) ? fProtonSaturation : realE),
                              realNf,
                              realNs,
                              time);
            }
        }
        else
        {
            R3BLOG(error, "CrystalID: " << tempCryID << " not valid.");
        }
    }
}

void R3BCalifaDigitizer::FillIdealConfigWithProtonRange(int nCrystalCals)
{
    double tempE = 0;
    double realE = 0;
    double tempNf = 0;
    double realNf = 0;
    double tempNs = 0;
    double realNs = 0;
    int tempCryID = 0;

    for (int i = 0; i < nCrystalCals; i++)
    {
        tempCryID = (dynamic_cast<R3BCalifaCrystalCalData*>(fCalifaCryCalDataCA->At(i)))->GetCrystalId();

        tempE = (dynamic_cast<R3BCalifaCrystalCalData*>(fCalifaCryCalDataCA->At(i)))->GetEnergy();
        tempNf = (dynamic_cast<R3BCalifaCrystalCalData*>(fCalifaCryCalDataCA->At(i)))->GetNf();
        tempNs = (dynamic_cast<R3BCalifaCrystalCalData*>(fCalifaCryCalDataCA->At(i)))->GetNs();

        realE = tempE;
        realNf = tempNf;
        realNs = tempNs;

        if (tempE < fThreshold)
        {
            fCalifaCryCalDataCA->RemoveAt(i);
            fCalifaCryCalDataCA->Compress();
            nCrystalCals--; // remove from CalData those below threshold
            i--;
            continue;
        }

        if (tempCryID <= fNumCrystals)
        {
            // Gamma branch
            if (tempCryID > barrelCrystals)
            {
                if (fGammaResolution > 0)
                {
                    fResolution = fGammaResolution;
                    realE = ExpResSmearing(tempE);
                }

                if (tempE >= fGammaSaturation && fGammaSaturation > 0 && fGammaResolution > 0)
                {
                    (dynamic_cast<R3BCalifaCrystalCalData*>(fCalifaCryCalDataCA->At(i)))
                        ->SetEnergy(ExpResSmearing(fGammaSaturation));
                    (dynamic_cast<R3BCalifaCrystalCalData*>(fCalifaCryCalDataCA->At(i)))
                        ->SetNf(ExpResSmearing(fGammaSaturation));
                    (dynamic_cast<R3BCalifaCrystalCalData*>(fCalifaCryCalDataCA->At(i)))
                        ->SetNs(ExpResSmearing(fGammaSaturation));
                    (dynamic_cast<R3BCalifaCrystalCalData*>(fCalifaCryCalDataCA->At(i)))->SetToTEnergy(realE);
                }
                else if (tempE >= fGammaSaturation && fGammaSaturation > 0 && fGammaResolution == 0)
                {
                    (dynamic_cast<R3BCalifaCrystalCalData*>(fCalifaCryCalDataCA->At(i)))->SetEnergy(fGammaSaturation);
                    (dynamic_cast<R3BCalifaCrystalCalData*>(fCalifaCryCalDataCA->At(i)))->SetNf(fGammaSaturation);
                    (dynamic_cast<R3BCalifaCrystalCalData*>(fCalifaCryCalDataCA->At(i)))->SetNs(fGammaSaturation);
                    (dynamic_cast<R3BCalifaCrystalCalData*>(fCalifaCryCalDataCA->At(i)))->SetToTEnergy(realE);
                }
                else if (fGammaResolution > 0)
                {
                    (dynamic_cast<R3BCalifaCrystalCalData*>(fCalifaCryCalDataCA->At(i)))->SetEnergy(realE);
                    (dynamic_cast<R3BCalifaCrystalCalData*>(fCalifaCryCalDataCA->At(i)))->SetNf(realNf);
                    (dynamic_cast<R3BCalifaCrystalCalData*>(fCalifaCryCalDataCA->At(i)))->SetNs(realNs);
                    (dynamic_cast<R3BCalifaCrystalCalData*>(fCalifaCryCalDataCA->At(i)))->SetToTEnergy(realE);
                }
                else if (fGammaResolution == 0)
                {
                    (dynamic_cast<R3BCalifaCrystalCalData*>(fCalifaCryCalDataCA->At(i)))->SetEnergy(realE);
                    (dynamic_cast<R3BCalifaCrystalCalData*>(fCalifaCryCalDataCA->At(i)))->SetNf(realNf);
                    (dynamic_cast<R3BCalifaCrystalCalData*>(fCalifaCryCalDataCA->At(i)))->SetNs(realNs);
                    (dynamic_cast<R3BCalifaCrystalCalData*>(fCalifaCryCalDataCA->At(i)))->SetToTEnergy(realE);
                }
            }
            // Proton branch
            if (fProtonResolution > 0)
            {
                fResolution = fProtonResolution;
                realE = ExpResSmearing(tempE);
            }
            auto time = (dynamic_cast<R3BCalifaCrystalCalData*>(fCalifaCryCalDataCA->At(i)))->GetTime();

            AddCrystalCal(tempCryID + fNumCrystals,
                          ((realE >= fProtonSaturation && fProtonSaturation > 0) ? fProtonSaturation : realE),
                          realNf,
                          realNs,
                          time);

            if (tempCryID <= barrelCrystals)
            {
                fCalifaCryCalDataCA->RemoveAt(i);
                fCalifaCryCalDataCA->Compress();
                nCrystalCals--; // remove from CalData those below threshold
                i--;
            }
        }
        else
        {
            R3BLOG(error, "CrystalID: " << tempCryID << " not valid.");
        }
    }
}

void R3BCalifaDigitizer::SetDetectionThreshold(double thresholdEne)
{
    fThreshold = thresholdEne;
    R3BLOG(info, "To " << fThreshold << " MeV");
}

void R3BCalifaDigitizer::SetRealConfig(bool isRealSet)
{
    fRealConfig = isRealSet;
    R3BLOG(info, "To " << isRealSet);
}

R3BCalifaCrystalCalData* R3BCalifaDigitizer::AddCrystalCal(int cryid,
                                                           double energy,
                                                           double Nf,
                                                           double Ns,
                                                           uint64_t time,
                                                           double tot_energy)
{
    if (fVerbose > 1)
    {
        R3BLOG(info,
               "Adding CrystalCalData with crystalID " << cryid << " entering with " << energy << " (MeV), Nf=" << Nf
                                                       << ", Ns=" << Ns << ", Time=" << time
                                                       << ", ToT energy=" << tot_energy);
    }
    if (fCalifaCryCalDataCA != nullptr)
    {
        return new ((*fCalifaCryCalDataCA)[fCalifaCryCalDataCA->GetEntriesFast()])
            R3BCalifaCrystalCalData(cryid, energy, Nf, Ns, time, tot_energy);
    }
    else
    {
        return nullptr;
    }
}

void R3BCalifaDigitizer::SetExpGammaEnergyRes(double crysRes)
{
    fGammaResolution = crysRes;
    R3BLOG(info, "To " << fGammaResolution << "% @ 1 MeV");
}

void R3BCalifaDigitizer::SetExpProtonEnergyRes(double crysRes)
{
    fProtonResolution = crysRes;
    R3BLOG(info, "To " << fProtonResolution << "% for protons");
}

void R3BCalifaDigitizer::SetComponentRes(double componentRes)
{
    fComponentRes = componentRes;
    R3BLOG(info, "To " << fComponentRes);
}

double R3BCalifaDigitizer::NUSmearing(double inputEnergy)
{
    // Very simple preliminary scheme where the NU is introduced as a flat random
    // distribution with limits fNonUniformity (%) of the energy value.
    //
    return gRandom->Uniform(inputEnergy - inputEnergy * fNonUniformity / 100.,
                            inputEnergy + inputEnergy * fNonUniformity / 100.);
}

void R3BCalifaDigitizer::SetNonUniformity(Double_t nonU)
{
    fNonUniformity = nonU;
    R3BLOG(info, "To " << fNonUniformity << " %");
}

double R3BCalifaDigitizer::ExpResSmearing(double inputEnergy)
{
    // Smears the energy according to some Experimental Resolution distribution
    // Very simple  scheme where the Experimental Resolution
    // is introduced as a gaus random distribution with a width given by the
    // parameter fResolution(in % @ MeV). Scales according to 1/sqrt(E)
    //
    // The formula is   TF1("name","0.058*x/sqrt(x)",0,10) for 3% at 1MeV (3.687 @
    // 662keV)
    //  ( % * energy ) / sqrt( energy )
    // and then the % is given at 1 MeV!!
    //
    if (inputEnergy < 0)
    {
        R3BLOG(error, "Negative value for energy: " << inputEnergy);
        return NAN;
    }
    else if (fResolution > 0)
    {
        // Energy in MeV
        double randomIs = gRandom->Gaus(0., inputEnergy * fResolution * 1. / (235.5 * sqrt(inputEnergy)));
        return inputEnergy + randomIs;
    }
    else
    {
        return inputEnergy;
    }
}

double R3BCalifaDigitizer::CompSmearing(double inputComponent)
{
    // Smears the components Ns and Nf according to fComponentRes
    //
    if (inputComponent < 0)
    {
        R3BLOG(error, "Negative value for Ns or Nf: " << inputComponent);
        return NAN;
    }
    else if (fComponentRes > 0)
    {
        // Energy in MeV
        double randomIs = gRandom->Gaus(0, inputComponent * fComponentRes * 1. / (235.5 * sqrt(inputComponent)));
        return inputComponent + randomIs;
    }
    else
    {
        return inputComponent;
    }
}

ClassImp(R3BCalifaDigitizer)
