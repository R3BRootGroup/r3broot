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

#include "R3BNeulandTimeRes.h"
#include "R3BNeulandCalData.h"
#include "TF1.h"
#include "TF2.h"
#include "TFile.h"
#include "TMath.h"
#include "TSpectrum.h"
#include <FairRootManager.h>
#include <TClonesArray.h>

#include <fstream>
#include <signal.h>
#include <sstream>

R3BNeulandTimeRes::R3BNeulandTimeRes()
    : FairTask("NeulandTimeRes", 1)
    , fUpdateRate(1000000)
    , fNEventsNeeded(10000)
    , fTrigger(-1)
{
}

R3BNeulandTimeRes::R3BNeulandTimeRes(const char* name, Int_t iVerbose)
    : FairTask(name, iVerbose)
    , fUpdateRate(1000000)
    , fNEventsNeeded(10000)
    , fTrigger(-1)
{
}

R3BNeulandTimeRes::~R3BNeulandTimeRes() {}

InitStatus R3BNeulandTimeRes::Init()
{
    std::cout << "init init" << std::endl;

    FairRootManager* rm = FairRootManager::Instance();
    if (!rm)
    {
        return kFATAL;
    }

    header = (R3BEventHeader*)rm->GetObject("EventHeader.");
    if (!header)
    {
        return kFATAL;
    }

    fPmt = (TClonesArray*)rm->GetObject("NeulandCalData");
    if (!fPmt)
    {
        return kFATAL;
    }

    for (Int_t pln = 0; pln < fNofPlanes; pln++)
    {
        for (Int_t bar = 0; bar < fNofBarsPerPlane; bar++)
        {

            bar_done[pln][bar] = 0;

            std::ostringstream hss;
            hss << "h_p" << pln + 1 << "b" << bar + 1;
            TString hname = hss.str();

            std::ostringstream hssq;
            hssq << "hq_p" << pln + 1 << "b" << bar + 1;
            TString hnameq = hssq.str();

            hTimeRes[pln][bar] = new TH1F(hname, hname, 1000, -20., 20.);
            hTimeResQ[pln][bar] = new TH2F(hnameq, hnameq, 400, 0, 1000, 400, -20., 20.);
        }
    }

    fNEvents = 0;
    bars_done = 0;
    finished = 0;

    return kSUCCESS;
}

void R3BNeulandTimeRes::Exec(Option_t* option)
{

    fNEvents++;

    // std::cout << fNEvents << std::endl;

    if (finished)
    {
        // if (fNEvents == fNEventsNeeded) {
        std::cout << std::endl;
        std::cout << "done 1" << std::endl;
        raise(SIGINT);
    }

    if (fTrigger >= 0)
    {
        if (header->GetTrigger() != fTrigger)
        {
            return;
        }
    }

    Int_t nHits = fPmt->GetEntries();

    // if (nHits < 4) return; // better cosmic peaks

    Double_t fCharge1[300] = { 0. }, fCharge2[300] = { 0. };
    Double_t fTime1[300] = { 0. }, fTime2[300] = { 0. };

    Int_t mult[fNofPlanes];

    for (Int_t pl = 0; pl < fNofPlanes; pl++)
        mult[pl] = 0;

    // Loop over mapped hits
    for (Int_t i = 0; i < nHits; i++)
    {

        R3BNeulandCalData* hit = (R3BNeulandCalData*)fPmt->At(i);

        if (!hit)
            continue; // should not happen

        Int_t iPlane = (hit->GetBarId() - 1) / 50;
        Int_t iBar = (hit->GetBarId() - 1) % 50;
        Int_t iSide = hit->GetSide() - 1;

        if (hit->GetTime() > 0)
            mult[iPlane]++;

        if (0 == iSide)
        {
            fCharge1[iPlane * 50 + iBar] = hit->GetQdc();
            fTime1[iPlane * 50 + iBar] = hit->GetTime() - wlk(fCharge1[iPlane * 50 + iBar]);
            // fTime1[iPlane*50+iBar] = hit->GetTime();
        }
        else
        {
            fCharge2[iPlane * 50 + iBar] = hit->GetQdc();
            fTime2[iPlane * 50 + iBar] = hit->GetTime() - wlk(fCharge2[iPlane * 50 + iBar]);
            // fTime2[iPlane*50+iBar] = hit->GetTime();
        }
    }

    for (Int_t pl = 0; pl < fNofPlanes; pl++)
    {
        for (Int_t bar = 1; bar < fNofBarsPerPlane; bar++)
        {

            if (bar_done[pl][bar] == 1)
                continue;

            Double_t t1 = fTime1[pl * 50 + bar - 1];
            Double_t t2 = fTime2[pl * 50 + bar - 1];
            Double_t t3 = fTime1[pl * 50 + bar];
            Double_t t4 = fTime2[pl * 50 + bar];

            Double_t qdc12 = sqrt(fCharge1[pl * 50 + bar - 1] * fCharge2[pl * 50 + bar - 1]);
            Double_t qdc34 = sqrt(fCharge1[pl * 50 + bar] * fCharge2[pl * 50 + bar]);

            // fill time diff histos
            if (t1 > 0 && t2 > 0 && t3 > 0 && t4 > 0 && qdc12 > 75. && qdc34 > 75.)
            {
                // if (t1 > 0 && t2 > 0 && t3 > 0 && t4 > 0) {
                if (pl % 2 == 0 && mult[pl] > 10)
                {
                    hTimeRes[pl][bar]->Fill((t3 + t4) / 2. - (t1 + t2) / 2.);
                    hTimeResQ[pl][bar]->Fill(qdc34, (t3 + t4) / 2. - (t1 + t2) / 2.);
                }
                if (pl % 2 == 1 && mult[pl] > 10)
                {
                    hTimeRes[pl][bar]->Fill((t3 + t4) / 2. - (t1 + t2) / 2.);
                    hTimeResQ[pl][bar]->Fill(qdc34, (t3 + t4) / 2. - (t1 + t2) / 2.);
                }
            }

            if (hTimeRes[pl][bar]->GetEntries() > 3000)
            {

                bar_done[pl][bar] = 1;

                bars_done++;

                std::cout << std::endl;
                std::cout << "bars done: " << bars_done << std::endl;

                // do fitting
                Int_t a = hTimeRes[pl][bar]->GetEntries();

                Double_t ampli = hTimeRes[pl][bar]->GetMaximum();
                Int_t binmax = hTimeRes[pl][bar]->GetMaximumBin();
                Double_t posi = hTimeRes[pl][bar]->GetXaxis()->GetBinCenter(binmax);

                Double_t sigma = 0.150;
                Double_t min = posi - 1.0;
                Double_t max = posi + 1.0;

                std::cout << "position of peak: " << posi << "  " << ampli << std::endl;

                std::ostringstream fss;
                fss << "f_p" << pl + 1 << "b" << bar + 1;
                TString fname = fss.str();

                if (pl % 2 == 0)
                {
                    fitfunc[pl][bar] = new TF1(fname, "gaus");
                    fitfunc[pl][bar]->SetParameters(ampli, posi, sigma);
                }
                else
                {
                    fitfunc[pl][bar] = new TF1(fname, "gaus(0)+gaus(3)");
                    fitfunc[pl][bar]->SetParameters(ampli, posi - 0.3, sigma, ampli, posi + 0.3, sigma);
                }

                hTimeRes[pl][bar]->Fit(fitfunc[pl][bar], "", "", min, max);

                Double_t nmax = fitfunc[pl][bar]->GetParameter(0);
                Double_t mean = fitfunc[pl][bar]->GetParameter(1);
                Double_t tres = fitfunc[pl][bar]->GetParameter(2);

                delete fitfunc[pl][bar];

                std::ofstream timeres_file("timeres.txt", std::ios::app);

                timeres_file << "time res p" << pl + 1 << "b" << bar + 1 << "  " << nmax << "   " << mean << "   "
                             << tres << std::endl;

                timeres_file.close();
            }

            if (bars_done == fNofPlanes * (fNofBarsPerPlane - 1))
                finished = 1;
        }
    }
}

void R3BNeulandTimeRes::FinishEvent() {}

void R3BNeulandTimeRes::FinishTask()
{

    std::cout << "done 1.5" << std::endl;

    TFile* fff = new TFile("histos_timeres.root", "recreate");

    for (Int_t pln = 0; pln < fNofPlanes; pln++)
    {
        for (Int_t bar = 0; bar < fNofBarsPerPlane; bar++)
        {
            hTimeRes[pln][bar]->Write();
            hTimeResQ[pln][bar]->Write();
        }
    }
    fff->Close();
}

Double_t R3BNeulandTimeRes::wlk(Double_t x)
{
    Double_t y = 0;
    Double_t par1 = 3.34290e+01; // not considered
    Double_t par2 = -4.48725e+05;
    Double_t par3 = -2.07310e+02;
    Double_t par4 = 2.48974e+02;
    Double_t par5 = 6.36355e-03;
    Double_t par6 = -2.93865e-05;

    if (x > 0 && x < 400)
        y = par2 * TMath::Power(x, par3) + par4 / x + par5 * x + par6 * x * x;
    else
        y = 0;

    return y;
    // return 0.;
}
