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

#include "R3BNeulandGainMatching.h"
#include "R3BNeulandCalData.h"
#include "TF1.h"
#include "TFile.h"
#include "TSpectrum.h"
#include <FairRootManager.h>
#include <TClonesArray.h>

#include <fstream>
#include <sstream>

namespace
{
    Float_t const MAX_STEPS = 10;
    Float_t const Kp = 1.0;   // Adjust PID here !!
    Float_t const Ki = 0.075; // Adjust PID here !!
    Float_t const Kd = 0.025; // Adjust PID here !!
    Float_t const Ta = 0.5;   // Adjust PID here !!
    Float_t const xtargetp1 = 95.;
    Float_t const accuracy = 0.01;
    Float_t const starthv = 1100.;

    Int_t const LSEARCHNB = 0;
    Int_t const HSEARCHNB = 200;
    Int_t const LGETCOS = 81;
    Int_t const HGETCOS = 200;
} // namespace

R3BNeulandGainMatching::R3BNeulandGainMatching()
    : FairTask("NeulandGainMatching", 1)
    , fUpdateRate(1000000)
    , fNEventsNeeded(10000)
    , fTrigger(-1)
{
}

R3BNeulandGainMatching::R3BNeulandGainMatching(const char* name, Int_t iVerbose)
    : FairTask(name, iVerbose)
    , fUpdateRate(1000000)
    , fNEventsNeeded(10000)
    , fTrigger(-1)
{
}

R3BNeulandGainMatching::~R3BNeulandGainMatching() {}

InitStatus R3BNeulandGainMatching::Init()
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

    // epics Neuland HV channels
    for (Int_t pln = 0; pln < fNofPlanes; pln++)
    {
        for (Int_t bar = 0; bar < fNofBarsPerPlane; bar++)
        {
            for (Int_t pmt = 0; pmt < 2; pmt++)
            {

                std::ostringstream oss;
                oss << "r3b:nl:hv:p" << pln + 1 << "b" << bar + 1 << "t" << pmt + 1;

                auto vmon = oss.str() + ":vmon";
                auto vtarget = oss.str() + ":vtarget";

                auto& entry = ca[pln][bar][pmt];
                entry.group = epics.CreateGroup();

                entry.vmon = entry.group->CreateChannel(vmon);
                entry.vtarget = entry.group->CreateChannel(vtarget);

                // std::cout << "r3b:nl:hv:p" << pln+1 << "b" << bar+1 << "t" << pmt+1 << std::endl;
                entry.group->Fetch();
                hv[pln][bar][pmt] = entry.vmon->Get();
                // std::cout << entry.vmon->Get() << std::endl;

                std::ostringstream hss;
                hss << "h_p" << pln + 1 << "b" << bar + 1 << "t" << pmt + 1;

                TString hname = hss.str();
                hCosmicPeak[pln][bar][pmt] = new TH1F(hname, hname, 300, 0, 600);

                iteration[pln][bar][pmt] = 0;
                esum[pln][bar][pmt] = 0.;
                ealt[pln][bar][pmt] = 0.;
            }
        }
    }

    // std::cout << ca[1][23][0].vmon->Get() << std::endl;

    finished = false;

    return kSUCCESS;
}

void R3BNeulandGainMatching::Exec(Option_t* option)
{

    Double_t maxhv = 1400;
    // check high voltage

    if (finished)
        FinishTask();

    if (fTrigger >= 0)
    {
        if (header->GetTrigger() != fTrigger)
        {
            return;
        }
    }

    Int_t nHits = fPmt->GetEntries();

    Int_t checkiteration = 0;

    for (Int_t pln = 0; pln < fNofPlanes; pln++)
    {
        for (Int_t bar = 0; bar < fNofBarsPerPlane; bar++)
        {
            for (Int_t pmt = 0; pmt < 2; pmt++)
            {
                if (iteration[pln][bar][pmt] == MAX_STEPS)
                {
                    checkiteration++;
                    // if (checkiteration > 0)
                    //  std::cout << "Done: p" << pln+1 << "b" << bar+1 << "t" << pmt+1 << std::endl;
                }
            }
        }
    }

    if (checkiteration == fNofPMTs)
        finished = 1;
    // std::cout << "Channels done: " << checkiteration << std::endl;

    if (nHits < 10)
        return; // better cosmic peaks

    // Loop over mapped hits
    for (Int_t i = 0; i < nHits; i++)
    {

        R3BNeulandCalData* hit = (R3BNeulandCalData*)fPmt->At(i);

        if (!hit)
            continue; // should not happen

        Int_t iPlane = (hit->GetBarId() - 1) / 50;
        Int_t iBar = (hit->GetBarId() - 1) % 50;
        Int_t iSide = hit->GetSide() - 1;
        Double_t fCharge = hit->GetQdc();

        // fill histos with cosmics peaks
        hCosmicPeak[iPlane][iBar][iSide]->Fill(fCharge);

        // check entries in histograms
        if (hCosmicPeak[iPlane][iBar][iSide]->GetEntries() > fNEventsNeeded &&
            iteration[iPlane][iBar][iSide] < MAX_STEPS)
        {
            Int_t a = hCosmicPeak[iPlane][iBar][iSide]->GetEntries();
            Int_t b = iteration[iPlane][iBar][iSide];
            std::cout << "Iteration: " << iteration[iPlane][iBar][iSide] << std::endl;

            std::ofstream hv_file("gainmatched_hv.txt", std::ios::app);

            Double_t newhv = 0;
            peakmethod = 0;

            auto& hventry = ca[iPlane][iBar][iSide];

            Double_t xtarget;
            if (iPlane % 2 == 1)
            {
                xtarget = xtargetp1 * 1.16;
            }
            else
            {
                xtarget = xtargetp1; // more energy deposit for cosmics in vertical paddles
            }
            Double_t e = xtarget - getcosmicpeak(hCosmicPeak[iPlane][iBar][iSide]);
            if (e == xtarget - 10000)
            {
                // Peak finding failed. Set HV back to last or start value.
                if (hv[iPlane][iBar][iSide] < starthv && hv[iPlane][iBar][iSide] >= starthv - 200)
                {
                    std::cout << "failed, try again with same hv" << std::endl;
                }
                else
                {
                    hv[iPlane][iBar][iSide] = starthv;
                    std::cout << "failed, back to start hv: " << hv[iPlane][iBar][iSide] << std::endl;
                }

                maxhv = 1400;
                if (iPlane < 2)
                    maxhv = 1800;
                if (hv[iPlane][iBar][iSide] >= 0 && hv[iPlane][iBar][iSide] <= maxhv)
                {
                    hventry.vtarget->Set(hv[iPlane][iBar][iSide]);
                    hventry.group->Commit();
                }
                hCosmicPeak[iPlane][iBar][iSide]->Reset();
                iteration[iPlane][iBar][iSide]++;
            }
            else if ((e <= xtarget * accuracy) && (e >= -xtarget * accuracy))
            {
                // Peak finding worked. Matching finished
                newhv = hv[iPlane][iBar][iSide];
                std::cout << "Matching for this pm finished. New HV: " << newhv << std::endl;
                iteration[iPlane][iBar][iSide] = MAX_STEPS;
            }
            else
            {
                // Peak finding worked.
                esum[iPlane][iBar][iSide] = esum[iPlane][iBar][iSide] + e;
                hventry.group->Fetch();
                std::cout << "old hv: " << hventry.vmon->Get() << std::endl;
                std::cout << "p" << iPlane + 1 << "b" << iBar + 1 << "t" << iSide + 1 << "  e: " << e
                          << "  esum: " << esum[iPlane][iBar][iSide] << "  ealt: " << ealt[iPlane][iBar][iSide]
                          << std::endl;

                hv[iPlane][iBar][iSide] = hv[iPlane][iBar][iSide] + Kp * e + Ki * Ta * esum[iPlane][iBar][iSide] +
                                          Kd * (e - ealt[iPlane][iBar][iSide]) / Ta;
                ealt[iPlane][iBar][iSide] = e;
                std::cout << "new hv: " << hv[iPlane][iBar][iSide] << std::endl;
                maxhv = 1400;
                if (iPlane < 2)
                    maxhv = 1800;
                if (hv[iPlane][iBar][iSide] >= 0 && hv[iPlane][iBar][iSide] <= maxhv)
                {
                    hventry.vtarget->Set(hv[iPlane][iBar][iSide]);
                    hventry.group->Commit();
                }
                hCosmicPeak[iPlane][iBar][iSide]->Reset();
                iteration[iPlane][iBar][iSide]++;
            }
            if (iteration[iPlane][iBar][iSide] == MAX_STEPS)
            {
                // Reached last Iteration -> abort
                if (newhv == 0)
                {
                    if ((e <= xtarget * accuracy * 5) && (e >= -xtarget * accuracy * 5))
                    {
                        std::cout
                            << "Reached max steps but last check was already inside 5% accuracy. Taking last hv value."
                            << std::endl;
                        newhv = hv[iPlane][iBar][iSide];
                    }
                    else
                    {
                        std::cout << "Error, check matching by hand. Voltage set to default value" << std::endl;
                        newhv = 1000;
                        hventry.vtarget->Set(newhv);
                        hventry.group->Commit();
                    }
                }
                hv_file << "caput r3b:nl:hv:p" << iPlane + 1 << "b" << iBar + 1 << "t" << iSide + 1 << ":vtarget "
                        << newhv << " " << peakmethod << " " << b << std::endl;
            }
            hv_file.close();
            // End gain matching
        }
    }
}

Double_t R3BNeulandGainMatching::getcosmicpeak(TH1* hin)
{

    Double_t xp = 10000;
    Double_t xpf = 0;
    Double_t sig = 0;
    Double_t max = 0;
    TF1* fitgaus;
    // std::cout << "xtarget: " << xtarget << std::endl;
    for (Int_t q = 0; q < 5; q++)
    {
        xpf = searchcosmicpeaknb(hin, 15. + q * 5.);
        if (xpf != 0)
        {
            xp = xpf;
            peakmethod = 1;
            break;
        }
        xpf = searchcosmicpeak(hin, 20. + q * 5.);
        if (xpf != 0)
        {
            xp = xpf;
            peakmethod = 2;
            break;
        }
    }

    if (xp == 10000)
    {
        hin->GetXaxis()->SetRange(LGETCOS, HGETCOS); // ig 81<  <200
        max = hin->GetMaximumBin();
        for (int k = 0; k < 5; k++)
        {
            fitgaus = new TF1("fitgaus", "gaus(0)", max - (10 + k * 3), max + (40 + k * 5));
            fitgaus->SetParameters(200, max, 20);
            hin->Fit("fitgaus", "QR0");
            xpf = fitgaus->GetParameter(1);
            sig = fitgaus->GetParameter(2);
            delete fitgaus;
            if (70 < xpf && xpf < 180 && 10 < sig && sig < 50)
            {
                xp = xpf;
                std::cout << "Found cosmic peak by gaus fitting at x = " << xp << std::endl;
                peakmethod = 3;
                break;
            }
        }
        if (xp == 10000)
        {
            if (70 < max && max < 180)
            {
                xp = max;
                std::cout << "Using maximum bin as cosmic peak position at x = " << xp << std::endl;
                peakmethod = 4;
            }
        }
        hin->GetXaxis()->SetRange(1, 600);
        // hin->GetXaxis()->SetDefaults();
    }
    return xp;
}

Double_t R3BNeulandGainMatching::searchcosmicpeak(TH1* hin, Double_t width)
{

    Double_t xp = 0;
    TSpectrum* s = new TSpectrum(3);
    Int_t nfound = s->Search(hin, width, "", 0.05);
    std::cout << "searchcosmicpeak: Found " << nfound << " candidate peaks" << std::endl;

    Double_t* xpeaks = s->GetPositionX();
    for (Int_t p = 0; p < nfound; p++)
    {
        // std::cout << "Peak " << p << " at x=" << xpeaks[p] << std::endl;
        if (75 < xpeaks[p] && xpeaks[p] < 180)
        {
            xp = xpeaks[p];
            std::cout << "searchcosmicpeak: Found cosmic peak at x = " << xp << std::endl;
            break;
        }
    }
    delete s;
    return xp;
}

Double_t R3BNeulandGainMatching::searchcosmicpeaknb(TH1* hin, Double_t width)
{

    Int_t max = hin->GetMaximumBin();
    hin->GetXaxis()->SetRange(LSEARCHNB, HSEARCHNB); // ig 60<  <200
    Double_t xp = 0;
    TSpectrum* s = new TSpectrum(3);
    Int_t nfound = s->Search(hin, width, "nobackground new", 0.05);
    Int_t diff = 0;
    Int_t difflast = 0;
    Int_t icall = 0;
    std::cout << "searchcosmicpeakNB: Found " << nfound << " candidate peaks" << std::endl;
    std::cout << "searchcosmicpeakNB: maximum bin content at = " << max << std::endl;
    Double_t* xpeaks = s->GetPositionX();
    for (Int_t p = 0; p < nfound; p++)
    {
        std::cout << "Peak " << p << " at x=" << xpeaks[p] << std::endl;
        if (LSEARCHNB < xpeaks[p] && xpeaks[p] < HSEARCHNB)
        { // ig 60<  <200
            // xp = xpeaks[p];
            diff = abs(xpeaks[p] - max);
            std::cout << "searchcosmicpeakNB: diff = " << diff << std::endl;
            if (0 == icall)
            {
                xp = xpeaks[p];
                difflast = diff;
                icall = 1;
            }
            if (diff < difflast)
            {
                std::cout << "searchcosmicpeakNB: diff = " << diff << std::endl;
                std::cout << "searchcosmicpeakNB: difflast = " << difflast << std::endl;
                xp = xpeaks[p];
                difflast = diff;
            }
        }
        //}
        // std::cout << "searchcosmicpeakNB: Found cosmic peak at x = " << xp << std::endl;
        // break;
    }
    std::cout << "searchcosmicpeakNB: Chose cosmic peak at x = " << xp << std::endl;
    //}
    delete s;
    hin->GetXaxis()->SetRange(1, 600);
    return xp;
}

void R3BNeulandGainMatching::FinishEvent() {}

void R3BNeulandGainMatching::FinishTask()
{

    TFile* fff = new TFile("histos.root", "recreate");

    for (Int_t pln = 0; pln < fNofPlanes; pln++)
    {
        for (Int_t bar = 0; bar < fNofBarsPerPlane; bar++)
        {
            for (Int_t pmt = 0; pmt < 2; pmt++)
            {
                hCosmicPeak[pln][bar][pmt]->Write();
            }
        }
    }
    fff->Close();
}
