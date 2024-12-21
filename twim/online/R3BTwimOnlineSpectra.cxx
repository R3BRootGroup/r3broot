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

// ------------------------------------------------------------
// -----               R3BTwimOnlineSpectra               -----
// -----    Created 29/09/19 by J.L. Rodriguez-Sanchez    -----
// -----           Fill TWIM online histograms            -----
// ------------------------------------------------------------

/*
 * This task should fill histograms with TWIM online data
 */

#include "R3BTwimOnlineSpectra.h"
#include "R3BEventHeader.h"
#include "R3BLogger.h"
#include "R3BMwpcHitData.h"
#include "R3BShared.h"
#include "R3BTwimCalData.h"
#include "R3BTwimHitData.h"
#include "R3BTwimMappedData.h"
#ifdef DSOFIADEFINED
#include "R3BSofTofWHitData.h"
#endif

#include "FairLogger.h"
#include "FairRootManager.h"
#include "FairRunAna.h"
#include "FairRunOnline.h"
#include "FairRuntimeDb.h"
#include "TCanvas.h"
#include "TClonesArray.h"
#include "TFolder.h"
#include "TH1F.h"
#include "TH2F.h"
#include "THttpServer.h"
#include "TMath.h"
#include "TRandom.h"
#include "TVector3.h"
#include <fstream>
#include <iostream>
#include <sstream>

R3BTwimOnlineSpectra::R3BTwimOnlineSpectra()
    : R3BTwimOnlineSpectra("TwimOnlineSpectra", 1)
{
}

R3BTwimOnlineSpectra::R3BTwimOnlineSpectra(const TString& name, Int_t iVerbose)
    : FairTask(name, iVerbose)
{
}

InitStatus R3BTwimOnlineSpectra::Init()
{
    R3BLOG(info, "");

    // try to get a handle on the EventHeader. EventHeader may not be
    // present though and hence may be null. Take care when using.

    FairRootManager* rootManager = FairRootManager::Instance();
    R3BLOG_IF(fatal, rootManager == nullptr, "FairRootManager not found");

    FairRunOnline* run = FairRunOnline::Instance();
    run->GetHttpServer()->Register("", this);

    auto header = dynamic_cast<R3BEventHeader*>(rootManager->GetObject("EventHeader."));
    R3BLOG_IF(warn, header == nullptr, "EventHeader. not found");

    if (fExpId == 0) // Obtain global ExpId if it's not set locally.
    {
        fExpId = header->GetExpId();
        R3BLOG(info, "fExpId: " << fExpId);
    }

    // get access to mapped data of the TWIM
    fMappedItemsTwim = dynamic_cast<TClonesArray*>(rootManager->GetObject("TwimMappedData"));
    R3BLOG_IF(fatal, fMappedItemsTwim == nullptr, "TwimMappedData not found");

    // get access to cal data of the TWIM
    fCalItemsTwim = dynamic_cast<TClonesArray*>(rootManager->GetObject("TwimCalData"));
    R3BLOG_IF(warn, fCalItemsTwim == nullptr, "TwimCalData not found");

    // get access to hit data of the TWIM
    fHitItemsTwim = dynamic_cast<TClonesArray*>(rootManager->GetObject("TwimHitData"));
    R3BLOG_IF(warn, fHitItemsTwim == nullptr, "TwimHitData not found");

    // get access to hit data of the MWPC3
    fHitItemsMwpc3 = dynamic_cast<TClonesArray*>(rootManager->GetObject("Mwpc3HitData"));
    R3BLOG_IF(warn, fHitItemsMwpc3 == nullptr, "Mwpc3HitData not found");

    // get access to hit data of the Tof-Wall
    fHitItemsTofW = dynamic_cast<TClonesArray*>(rootManager->GetObject("TofWHitData"));
    R3BLOG_IF(warn, fHitItemsTofW == nullptr, "TofWHitData not found");

    if (fExpId == 444 || fExpId == 467)
    {
        fNbSections = 1;
        fNbAnodes = 16;
        fNbTref = 2;
        fNbTrig = 2;
    }
    else if (fExpId == 509 || fExpId == 522 || fExpId == 91 || fExpId == 118)
    {
        fNbSections = 1;
        fNbAnodes = 16;
        fNbTref = 1;
        fNbTrig = 1;
    }
    else if (fExpId == 455)
    {
        fNbSections = 4;
        fNbAnodes = 16;
        fNbTref = 1;
        fNbTrig = 1;
    }

    // Create histograms for detectors
    char Name1[255];
    char Name2[255];

    cTwimMap_E = new TCanvas*[fNbSections];
    cTwimMap_T = new TCanvas*[fNbSections];
    cTwimMult = new TCanvas*[fNbSections];
    // TWIM: Map data for E and T
    for (Int_t i = 0; i < fNbSections; i++)
    {
        sprintf(Name1, "Twim_Emap_Sec_%d", i + 1);
        sprintf(Name2, "Section %d", i + 1);
        cTwimMap_E[i] = new TCanvas(Name1, Name2, 10, 10, 800, 700);
        cTwimMap_E[i]->Divide(4, 4);

        sprintf(Name1, "Twim_Tmap_Sec_%d", i + 1);
        sprintf(Name2, "Section %d", i + 1);
        cTwimMap_T[i] = new TCanvas(Name1, Name2, 10, 10, 800, 700);
        cTwimMap_T[i]->Divide(4, 4);

        sprintf(Name1, "Twim_Mult_Sec_%d", i + 1);
        sprintf(Name2, "Section %d", i + 1);
        cTwimMult[i] = new TCanvas(Name1, Name2, 10, 10, 800, 700);
        cTwimMult[i]->Divide(4, 4);
    }

    fh1_twimmap_E = new TH1F*[fNbSections * fNbAnodes];
    fh1_twimmap_T = new TH1F*[fNbSections * fNbAnodes];
    fh1_twimmult = new TH1F*[fNbSections * fNbAnodes];
    for (Int_t i = 0; i < fNbSections; i++)
    {
        for (Int_t j = 0; j < fNbAnodes; j++)
        {
            sprintf(Name1, "fh1_twim_Emap_sec%d_a%d", i + 1, j + 1);
            sprintf(Name2, "Sec %d:Anode %d", i + 1, j + 1);
            if (fExpId == 455)
                fh1_twimmap_E[i * fNbAnodes + j] = R3B::root_owned<TH1F>(Name1, Name2, 64000, 0, 64000);
            else
            {
                fh1_twimmap_E[i * fNbAnodes + j] = R3B::root_owned<TH1F>(Name1, Name2, 8192, 0, 8192);
            }
            fh1_twimmap_E[i * fNbAnodes + j]->GetXaxis()->SetTitle("Raw Energy [channels], if mult==1");
            fh1_twimmap_E[i * fNbAnodes + j]->GetYaxis()->SetTitle("Counts");
            fh1_twimmap_E[i * fNbAnodes + j]->GetYaxis()->SetTitleOffset(1.1);
            fh1_twimmap_E[i * fNbAnodes + j]->GetXaxis()->CenterTitle(true);
            fh1_twimmap_E[i * fNbAnodes + j]->GetYaxis()->CenterTitle(true);
            fh1_twimmap_E[i * fNbAnodes + j]->GetXaxis()->SetLabelSize(0.045);
            fh1_twimmap_E[i * fNbAnodes + j]->GetXaxis()->SetTitleSize(0.045);
            fh1_twimmap_E[i * fNbAnodes + j]->GetYaxis()->SetLabelSize(0.045);
            fh1_twimmap_E[i * fNbAnodes + j]->GetYaxis()->SetTitleSize(0.045);
            fh1_twimmap_E[i * fNbAnodes + j]->SetFillColor(31);
            cTwimMap_E[i]->cd(j + 1);
            fh1_twimmap_E[i * fNbAnodes + j]->Draw("");

            sprintf(Name1, "fh1_twim_Tmap_sec%d_a%d", i + 1, j + 1);
            sprintf(Name2, "Sec %d:Anode %d - T-Ttrig", i + 1, j + 1);
            fh1_twimmap_T[i * fNbAnodes + j] = R3B::root_owned<TH1F>(Name1, Name2, 25000, -50000, 50000);
            fh1_twimmap_T[i * fNbAnodes + j]->GetXaxis()->SetTitle("Raw Time - Trigger Time [channels]");
            fh1_twimmap_T[i * fNbAnodes + j]->GetYaxis()->SetTitle("Counts");
            fh1_twimmap_T[i * fNbAnodes + j]->GetYaxis()->SetTitleOffset(1.1);
            fh1_twimmap_T[i * fNbAnodes + j]->GetXaxis()->CenterTitle(true);
            fh1_twimmap_T[i * fNbAnodes + j]->GetYaxis()->CenterTitle(true);
            fh1_twimmap_T[i * fNbAnodes + j]->GetXaxis()->SetLabelSize(0.045);
            fh1_twimmap_T[i * fNbAnodes + j]->GetXaxis()->SetTitleSize(0.045);
            fh1_twimmap_T[i * fNbAnodes + j]->GetYaxis()->SetLabelSize(0.045);
            fh1_twimmap_T[i * fNbAnodes + j]->GetYaxis()->SetTitleSize(0.045);
            fh1_twimmap_T[i * fNbAnodes + j]->SetFillColor(31);
            cTwimMap_T[i]->cd(j + 1);
            fh1_twimmap_T[i * fNbAnodes + j]->Draw("");

            sprintf(Name1, "fh1_twim_mult_sec%d_a%d", i + 1, j + 1);
            sprintf(Name2, "Sec %d:Anode %d - Mult per event", i + 1, j + 1);
            fh1_twimmult[i * fNbAnodes + j] = R3B::root_owned<TH1F>(Name1, Name2, 10, -0.5, 9.5);
            fh1_twimmult[i * fNbAnodes + j]->GetXaxis()->SetTitle("mult per event");
            fh1_twimmult[i * fNbAnodes + j]->GetYaxis()->SetTitle("Counts");
            fh1_twimmult[i * fNbAnodes + j]->GetYaxis()->SetTitleOffset(1.1);
            fh1_twimmult[i * fNbAnodes + j]->GetXaxis()->CenterTitle(true);
            fh1_twimmult[i * fNbAnodes + j]->GetYaxis()->CenterTitle(true);
            fh1_twimmult[i * fNbAnodes + j]->GetXaxis()->SetLabelSize(0.045);
            fh1_twimmult[i * fNbAnodes + j]->GetXaxis()->SetTitleSize(0.045);
            fh1_twimmult[i * fNbAnodes + j]->GetYaxis()->SetLabelSize(0.045);
            fh1_twimmult[i * fNbAnodes + j]->GetYaxis()->SetTitleSize(0.045);
            fh1_twimmult[i * fNbAnodes + j]->SetFillColor(31);
            cTwimMult[i]->cd(j + 1);
            fh1_twimmult[i * fNbAnodes + j]->Draw("");
        }
    }

    // TWIM: Map data for E vs DT
    cTwimMap_EvsDT = new TCanvas*[fNbSections];
    cTwim_DTvsDT = new TCanvas*[fNbSections];
    for (Int_t i = 0; i < fNbSections; i++)
    {
        sprintf(Name1, "Twim_EvsDT_Sec_%d", i + 1);
        sprintf(Name2, "Section %d", i + 1);
        cTwimMap_EvsDT[i] = new TCanvas(Name1, Name2, 10, 10, 800, 700);
        cTwimMap_EvsDT[i]->Divide(4, 4);
        sprintf(Name1, "Twim_DTvsDT_Sec_%d", i + 1);
        sprintf(Name2, "DTvsDT Section %d", i + 1);
        cTwim_DTvsDT[i] = new TCanvas(Name1, Name2, 10, 10, 800, 700);
        cTwim_DTvsDT[i]->Divide(4, 4);
    }

    char Name3[255];
    char Name4[255];
    cTwimMap_EsumvsDT = new TCanvas("Twim_EsumvsDT", "Twim_EsumvsDT", 10, 10, 800, 700);
    if (fNbSections > 1)
        cTwimMap_EsumvsDT->Divide(2, 2);

    fh2_twim_EneRawVsDriftTime = new TH2F*[fNbSections * fNbAnodes];
    fh2_twim_DTvsDT = new TH2F*[fNbSections * fNbAnodes];
    fh2_twim_EneRawSumVsDriftTime = new TH2F*[fNbSections];
    for (Int_t i = 0; i < fNbSections; i++)
    {
        for (Int_t j = 0; j < fNbAnodes; j++)
        {
            sprintf(Name1, "fh1_twim_EvsDT_sec%d_a%d", i + 1, j + 1);
            sprintf(Name2, "Sec %d:Anode %d", i + 1, j + 1);
            if (fExpId == 455)
                fh2_twim_EneRawVsDriftTime[i * fNbAnodes + j] =
                    R3B::root_owned<TH2F>(Name1, Name2, 1000, 0, 64000, 400, 1, 40001);
            else
                fh2_twim_EneRawVsDriftTime[i * fNbAnodes + j] =
                    R3B::root_owned<TH2F>(Name1, Name2, 500, 0, 8192, 800, 1, 40001);
            fh2_twim_EneRawVsDriftTime[i * fNbAnodes + j]->GetXaxis()->SetTitle(
                "Energy [channels, 8k ADC resolution], if mult==1");
            fh2_twim_EneRawVsDriftTime[i * fNbAnodes + j]->GetYaxis()->SetTitle(
                "Drift time [channels, 100ps TDC resolution], if mult==1");
            fh2_twim_EneRawVsDriftTime[i * fNbAnodes + j]->GetYaxis()->SetTitleOffset(1.1);
            fh2_twim_EneRawVsDriftTime[i * fNbAnodes + j]->GetXaxis()->CenterTitle(true);
            fh2_twim_EneRawVsDriftTime[i * fNbAnodes + j]->GetYaxis()->CenterTitle(true);
            fh2_twim_EneRawVsDriftTime[i * fNbAnodes + j]->GetXaxis()->SetLabelSize(0.045);
            fh2_twim_EneRawVsDriftTime[i * fNbAnodes + j]->GetXaxis()->SetTitleSize(0.045);
            fh2_twim_EneRawVsDriftTime[i * fNbAnodes + j]->GetYaxis()->SetLabelSize(0.045);
            fh2_twim_EneRawVsDriftTime[i * fNbAnodes + j]->GetYaxis()->SetTitleSize(0.045);
            cTwimMap_EvsDT[i]->cd(j + 1);
            fh2_twim_EneRawVsDriftTime[i * fNbAnodes + j]->Draw("col");
        }

        for (Int_t j = 0; j < fNbAnodes - 1; j++)
        {
            sprintf(Name1, "fh2_twim_DTvsDT_sec%d_a%d", i + 1, j);
            sprintf(Name2, "Sec %d: DT_%d vs DT_%d", i + 1, j + 1, j + 2);
            fh2_twim_DTvsDT[i * fNbAnodes + j] = R3B::root_owned<TH2F>(Name1, Name2, 800, 1, 40001, 200, -200, 200);
            sprintf(Name1, "Drift time %d [channels, 100ps TDC resolution]", j + 1);
            sprintf(Name2, "DT %d - DT %d [channels, 100ps TDC res.]", j + 1, j + 2);
            fh2_twim_DTvsDT[i * fNbAnodes + j]->GetXaxis()->SetTitle(Name1);
            fh2_twim_DTvsDT[i * fNbAnodes + j]->GetYaxis()->SetTitle(Name2);
            fh2_twim_DTvsDT[i * fNbAnodes + j]->GetYaxis()->SetTitleOffset(1.1);
            fh2_twim_DTvsDT[i * fNbAnodes + j]->GetXaxis()->CenterTitle(true);
            fh2_twim_DTvsDT[i * fNbAnodes + j]->GetYaxis()->CenterTitle(true);
            fh2_twim_DTvsDT[i * fNbAnodes + j]->GetXaxis()->SetLabelSize(0.045);
            fh2_twim_DTvsDT[i * fNbAnodes + j]->GetXaxis()->SetTitleSize(0.045);
            fh2_twim_DTvsDT[i * fNbAnodes + j]->GetYaxis()->SetLabelSize(0.045);
            fh2_twim_DTvsDT[i * fNbAnodes + j]->GetYaxis()->SetTitleSize(0.045);
            cTwim_DTvsDT[i]->cd(j + 1);
            fh2_twim_DTvsDT[i * fNbAnodes + j]->Draw("col");
        }

        sprintf(Name3, "fh2_Twim_EsumvsDT_sec%d", i + 1);
        sprintf(Name4, "Twim: Esum vs DT for sec%d", i + 1);
        if (fExpId == 455)
            fh2_twim_EneRawSumVsDriftTime[i] = R3B::root_owned<TH2F>(Name3, Name4, 1000, 1, 25000, 1000, 0, 64000);
        else
            fh2_twim_EneRawSumVsDriftTime[i] = R3B::root_owned<TH2F>(Name3, Name4, 1000, 1, 25000, 1000, 0, 8192);
        fh2_twim_EneRawSumVsDriftTime[i]->GetXaxis()->SetTitle("Drift time anode 5 [channels, 100ps TDC res.]");
        fh2_twim_EneRawSumVsDriftTime[i]->GetYaxis()->SetTitle("Energy Sum [channels]");
        fh2_twim_EneRawSumVsDriftTime[i]->GetYaxis()->SetTitleOffset(1.1);
        fh2_twim_EneRawSumVsDriftTime[i]->GetXaxis()->CenterTitle(true);
        fh2_twim_EneRawSumVsDriftTime[i]->GetYaxis()->CenterTitle(true);
        fh2_twim_EneRawSumVsDriftTime[i]->GetXaxis()->SetLabelSize(0.045);
        fh2_twim_EneRawSumVsDriftTime[i]->GetXaxis()->SetTitleSize(0.045);
        fh2_twim_EneRawSumVsDriftTime[i]->GetYaxis()->SetLabelSize(0.045);
        fh2_twim_EneRawSumVsDriftTime[i]->GetYaxis()->SetTitleSize(0.045);
        if (fNbSections > 1)
            cTwimMap_EsumvsDT->cd(i + 1);
        else
            cTwimMap_EsumvsDT->cd();
        fh2_twim_EneRawSumVsDriftTime[i]->Draw("col");
    }

    // TWIM: Map data for DT
    cTwimMap_DT = new TCanvas*[fNbSections];
    for (Int_t i = 0; i < fNbSections; i++)
    {
        sprintf(Name1, "Twim_DTmap_Sec_%d", i + 1);
        sprintf(Name2, "Section %d", i + 1);
        cTwimMap_DT[i] = new TCanvas(Name1, Name2, 10, 10, 800, 700);
        cTwimMap_DT[i]->Divide(4, 4);
    }

    fh1_twimmap_DT = new TH1F*[fNbSections * fNbAnodes];
    for (Int_t i = 0; i < fNbSections; i++)
    {
        for (Int_t j = 0; j < fNbAnodes; j++)
        {
            sprintf(Name1, "fh1_twim_DTmap_sec%d_a%d", i + 1, j + 1);
            sprintf(Name2, "Sec %d:Anode %d", i + 1, j + 1);
            fh1_twimmap_DT[i * fNbAnodes + j] =
                R3B::root_owned<TH1F>(Name1, Name2, 26000, -1000, 25000); // keep negative values to have the FG
            fh1_twimmap_DT[i * fNbAnodes + j]->GetXaxis()->SetTitle("Drift time [channels, 100ps TDC resolution]");
            fh1_twimmap_DT[i * fNbAnodes + j]->GetYaxis()->SetTitle("Counts");
            fh1_twimmap_DT[i * fNbAnodes + j]->GetYaxis()->SetTitleOffset(1.1);
            fh1_twimmap_DT[i * fNbAnodes + j]->GetXaxis()->CenterTitle(true);
            fh1_twimmap_DT[i * fNbAnodes + j]->GetYaxis()->CenterTitle(true);
            fh1_twimmap_DT[i * fNbAnodes + j]->GetXaxis()->SetLabelSize(0.045);
            fh1_twimmap_DT[i * fNbAnodes + j]->GetXaxis()->SetTitleSize(0.045);
            fh1_twimmap_DT[i * fNbAnodes + j]->GetYaxis()->SetLabelSize(0.045);
            fh1_twimmap_DT[i * fNbAnodes + j]->GetYaxis()->SetTitleSize(0.045);
            fh1_twimmap_DT[i * fNbAnodes + j]->SetFillColor(31);
            cTwimMap_DT[i]->cd(j + 1);
            fh1_twimmap_DT[i * fNbAnodes + j]->Draw("");
        }
    }

    // TWIM: Map data for multiplicities
    TCanvas* cTwim_Mult[fNbSections];
    fh1_Twimmap_mult = new TH1F*[fNbSections];
    for (Int_t i = 0; i < fNbSections; i++)
    {
        sprintf(Name1, "Twim_multiplicity_Sec_%d", i + 1);
        sprintf(Name2, "TWIM: Multiplicity per anode in section %d", i + 1);
        cTwim_Mult[i] = new TCanvas(Name1, Name2, 10, 10, 800, 700);
        sprintf(Name1, "fh1_Twim_mult_sec_%d", i + 1);
        sprintf(Name2, "Twim: Multiplicity section %d", i + 1);
        fh1_Twimmap_mult[i] =
            R3B::root_owned<TH1F>(Name1, Name2, fNbAnodes + fNbTref + fNbTrig, 1, fNbAnodes + fNbTref + fNbTrig + 1);
        fh1_Twimmap_mult[i]->GetXaxis()->SetTitle("Anode");
        fh1_Twimmap_mult[i]->GetYaxis()->SetTitle("Counts");
        fh1_Twimmap_mult[i]->GetYaxis()->SetTitleOffset(1.1);
        fh1_Twimmap_mult[i]->GetXaxis()->CenterTitle(true);
        fh1_Twimmap_mult[i]->GetYaxis()->CenterTitle(true);
        fh1_Twimmap_mult[i]->GetXaxis()->SetLabelSize(0.045);
        fh1_Twimmap_mult[i]->GetXaxis()->SetTitleSize(0.045);
        fh1_Twimmap_mult[i]->GetYaxis()->SetLabelSize(0.045);
        fh1_Twimmap_mult[i]->GetYaxis()->SetTitleSize(0.045);
        fh1_Twimmap_mult[i]->SetFillColor(31);
        fh1_Twimmap_mult[i]->Draw("");
    }

    // TWIM: Map data for E in the first (anodes from 0 to 7) and final sections (anodes from 8 to 15)
    cTwimMap_ESum = new TCanvas("twim_ESum_1,2", "twim_ESum_1,2", 10, 10, 800, 700);
    cTwimMap_ESum->Divide(1, 2);
    cTwimMap_ESum->cd(1);
    fh1_twim_ESum[0] = R3B::root_owned<TH1F>("fh1_twim_ESum1", "twim:ESum:first", 8192, 0, 8192);
    fh1_twim_ESum[0]->GetXaxis()->SetTitle("Energy [channels]");
    fh1_twim_ESum[0]->GetYaxis()->SetTitle("Counts");
    fh1_twim_ESum[0]->GetYaxis()->SetTitleOffset(1.1);
    fh1_twim_ESum[0]->GetXaxis()->CenterTitle(true);
    fh1_twim_ESum[0]->GetYaxis()->CenterTitle(true);
    fh1_twim_ESum[0]->GetXaxis()->SetLabelSize(0.045);
    fh1_twim_ESum[0]->GetXaxis()->SetTitleSize(0.045);
    fh1_twim_ESum[0]->GetYaxis()->SetLabelSize(0.045);
    fh1_twim_ESum[0]->GetYaxis()->SetTitleSize(0.045);
    fh1_twim_ESum[0]->SetFillColor(31);
    fh1_twim_ESum[0]->Draw("");
    cTwimMap_ESum->cd(2);
    fh1_twim_ESum[1] = R3B::root_owned<TH1F>("fh1_twim_ESum2", "twim:ESum:Second", 8192, 0, 8192);
    fh1_twim_ESum[1]->GetXaxis()->SetTitle("Energy [channels]");
    fh1_twim_ESum[1]->GetYaxis()->SetTitle("Counts");
    fh1_twim_ESum[1]->GetYaxis()->SetTitleOffset(1.1);
    fh1_twim_ESum[1]->GetXaxis()->CenterTitle(true);
    fh1_twim_ESum[1]->GetYaxis()->CenterTitle(true);
    fh1_twim_ESum[1]->GetXaxis()->SetLabelSize(0.045);
    fh1_twim_ESum[1]->GetXaxis()->SetTitleSize(0.045);
    fh1_twim_ESum[1]->GetYaxis()->SetLabelSize(0.045);
    fh1_twim_ESum[1]->GetYaxis()->SetTitleSize(0.045);
    fh1_twim_ESum[1]->SetFillColor(31);
    fh1_twim_ESum[1]->Draw("");

    // TWIM: Map data for Esum
    cTwimMap_ESum1 = new TCanvas("twim_ESum", "twim_ESum", 10, 10, 800, 700);
    fh1_twim_ESum[2] = R3B::root_owned<TH1F>("fh1_twim_ESum", "twim:ESum", 8192, 0, 8192);
    fh1_twim_ESum[2]->GetXaxis()->SetTitle("Energy [channels]");
    fh1_twim_ESum[2]->GetYaxis()->SetTitle("Counts");
    fh1_twim_ESum[2]->GetYaxis()->SetTitleOffset(1.1);
    fh1_twim_ESum[2]->GetXaxis()->CenterTitle(true);
    fh1_twim_ESum[2]->GetYaxis()->CenterTitle(true);
    fh1_twim_ESum[2]->GetXaxis()->SetLabelSize(0.045);
    fh1_twim_ESum[2]->GetXaxis()->SetTitleSize(0.045);
    fh1_twim_ESum[2]->GetYaxis()->SetLabelSize(0.045);
    fh1_twim_ESum[2]->GetYaxis()->SetTitleSize(0.045);
    fh1_twim_ESum[2]->SetFillColor(31);
    fh1_twim_ESum[2]->Draw("");

    // TWIM: Map data for Esum1 vs Esum2
    cTwimMap_ESum2 = new TCanvas("twim_E1vsE2", "twim_E1vsE2", 10, 10, 800, 700);
    fh2_twim_ESum = R3B::root_owned<TH2F>("fh2_twim_ESum", "twim: ESum1 vs Esum2", 1000, 0, 8192, 1000, 0, 8192);
    fh2_twim_ESum->GetXaxis()->SetTitle("Energy1 [channels]");
    fh2_twim_ESum->GetYaxis()->SetTitle("Energy2 [channels]");
    fh2_twim_ESum->GetYaxis()->SetTitleOffset(1.1);
    fh2_twim_ESum->GetXaxis()->CenterTitle(true);
    fh2_twim_ESum->GetYaxis()->CenterTitle(true);
    fh2_twim_ESum->GetXaxis()->SetLabelSize(0.045);
    fh2_twim_ESum->GetXaxis()->SetTitleSize(0.045);
    fh2_twim_ESum->GetYaxis()->SetLabelSize(0.045);
    fh2_twim_ESum->GetYaxis()->SetTitleSize(0.045);
    fh2_twim_ESum->Draw("col");

    // TWIM: Esum versus DT16-DT1
    cTwimMap_ESum_vs_diffDT = new TCanvas("twim_ESum_vs_diffDT", "twim_ESum_vs_diffDT", 10, 10, 800, 700);
    if (fNbSections > 1)
        cTwimMap_ESum_vs_diffDT->Divide(2, 2);

    fh2_twim_ESum_vs_diffDT = new TH2F*[fNbSections];
    for (Int_t i = 0; i < fNbSections; i++)
    {
        sprintf(Name1, "fh2_twim_ESum_vs_diffDT_Sec_%d", i + 1);
        sprintf(Name2, "TWIM: ESum vs diffDT in section %d", i + 1);
        if (fExpId == 455)
            fh2_twim_ESum_vs_diffDT[i] = R3B::root_owned<TH2F>(Name1, Name2, 1000, -2000, 2000, 1000, 0, 64000);
        else
            fh2_twim_ESum_vs_diffDT[i] = R3B::root_owned<TH2F>(Name1, Name2, 1000, -2000, 2000, 1000, 0, 8192);
        fh2_twim_ESum_vs_diffDT[i]->GetXaxis()->SetTitle("DT difference [channels]");
        fh2_twim_ESum_vs_diffDT[i]->GetYaxis()->SetTitle("Energy Sum [channels]");
        fh2_twim_ESum_vs_diffDT[i]->GetYaxis()->SetTitleOffset(1.1);
        fh2_twim_ESum_vs_diffDT[i]->GetXaxis()->CenterTitle(true);
        fh2_twim_ESum_vs_diffDT[i]->GetYaxis()->CenterTitle(true);
        fh2_twim_ESum_vs_diffDT[i]->GetXaxis()->SetLabelSize(0.045);
        fh2_twim_ESum_vs_diffDT[i]->GetXaxis()->SetTitleSize(0.045);
        fh2_twim_ESum_vs_diffDT[i]->GetYaxis()->SetLabelSize(0.045);
        fh2_twim_ESum_vs_diffDT[i]->GetYaxis()->SetTitleSize(0.045);
        if (fNbSections > 1)
            cTwimMap_ESum_vs_diffDT->cd(i + 1);
        else
            cTwimMap_ESum_vs_diffDT->cd();
        fh2_twim_ESum_vs_diffDT[i]->Draw("col");
    }

    // Twin_MUSIC: Map for delta t between Tref and Trig
    cTwimMap_DeltaTrefTrig = new TCanvas*[fNbSections];
    fh1_twimmap_DeltaTrefTrig = new TH1F*[fNbSections * fNbTref];
    for (Int_t i = 0; i < fNbSections; i++)
    {
        sprintf(Name1, "Twim_TrefTrig_Sec_%d", i + 1);
        sprintf(Name2, "Delta T (Tref-Trig) section %d", i + 1);
        cTwimMap_DeltaTrefTrig[i] = new TCanvas(Name1, Name2, 10, 10, 800, 700);
        cTwimMap_DeltaTrefTrig[i]->Divide(1, fNbTref);
        for (Int_t j = 0; j < fNbTref; j++)
        {
            sprintf(Name1, "fh1_Twim_Sec%d_DeltaTrefTrig%d", i + 1, j + 1);
            sprintf(Name2,
                    "Delta T (Tref%d-Trig%d) in section %d [channels, 100ps TDC resolution] 1ns/bin",
                    j + 1,
                    j + 1,
                    i + 1);
            fh1_twimmap_DeltaTrefTrig[i * fNbTref + j] = R3B::root_owned<TH1F>(Name1, Name2, 8000, -40000, 40000);
            cTwimMap_DeltaTrefTrig[i]->cd(j + 1);
            fh1_twimmap_DeltaTrefTrig[i * fNbTref + j]->Draw("");
        }
    }

    // Cal data, position in mm for each anode
    // Section 1: LEFT - DOWN,  x increasing from cathode to FG
    // Section 2: LEFT - UP,    x increasing from cathode to FG
    // Section 3: RIGHT - UP,   x increasing from FG to cathode
    // Section 4: RIGHT - DOWN, x increasing from FG to cathode
    // the distance from the FG to the cathode is 110 mm
    cTwimCal_Pos = new TCanvas*[fNbSections];
    fh1_Twimcal_Pos = new TH1F*[fNbSections * fNbAnodes];
    for (Int_t i = 0; i < fNbSections; i++)
    {
        sprintf(Name1, "Twim_Pos-X_Sec_%d", i + 1);
        sprintf(Name2, "Pos-X in mm section %d", i + 1);
        cTwimCal_Pos[i] = new TCanvas(Name1, Name2, 10, 10, 800, 700);
        cTwimCal_Pos[i]->Divide(4, 4);
        for (Int_t j = 0; j < fNbAnodes; j++)
        {
            sprintf(Name1, "fh1_Twim_Sec%d_Pos_a%d", i + 1, j + 1);
            sprintf(Name2, "Sec %d, Anode %d", i + 1, j + 1);
            if (i == 0 || i == 1)
            {
                fh1_Twimcal_Pos[i * fNbAnodes + j] = R3B::root_owned<TH1F>(Name1, Name2, 1400, 120, -10.);
                fh1_Twimcal_Pos[i * fNbAnodes + j]->GetXaxis()->SetTitle("LEFT / MESSEL Position-X [mm]");
            }
            else
            {
                fh1_Twimcal_Pos[i * fNbAnodes + j] = R3B::root_owned<TH1F>(Name1, Name2, 1400, -120., 10.);
                fh1_Twimcal_Pos[i * fNbAnodes + j]->GetXaxis()->SetTitle("RIGHT / WIXHAUSEN Position-X [mm]");
            }
            fh1_Twimcal_Pos[i * fNbAnodes + j]->GetYaxis()->SetTitle("Counts");
            fh1_Twimcal_Pos[i * fNbAnodes + j]->GetYaxis()->SetTitleOffset(1.1);
            fh1_Twimcal_Pos[i * fNbAnodes + j]->GetXaxis()->CenterTitle(true);
            fh1_Twimcal_Pos[i * fNbAnodes + j]->GetYaxis()->CenterTitle(true);
            fh1_Twimcal_Pos[i * fNbAnodes + j]->GetXaxis()->SetLabelSize(0.045);
            fh1_Twimcal_Pos[i * fNbAnodes + j]->GetXaxis()->SetTitleSize(0.045);
            fh1_Twimcal_Pos[i * fNbAnodes + j]->GetYaxis()->SetLabelSize(0.045);
            fh1_Twimcal_Pos[i * fNbAnodes + j]->GetYaxis()->SetTitleSize(0.045);
            fh1_Twimcal_Pos[i * fNbAnodes + j]->SetFillColor(31);
            fh1_Twimcal_Pos[i * fNbAnodes + j]->SetLineColor(1);
            cTwimCal_Pos[i]->cd(j + 1);
            fh1_Twimcal_Pos[i * fNbAnodes + j]->Draw("");
        }
    }

    // Hit data
    // s444 and s467 experiment: init of histograms ------
    TCanvas* cTwim_Z = new TCanvas("Twim_charge_z", "Twim: Charge Z", 10, 10, 800, 700);
    fh1_Twimhit_z = R3B::root_owned<TH1F>("fh1_Twim_charge_z", "Twim: Charge Z", 1200, 6, 38);
    fh1_Twimhit_z->GetXaxis()->SetTitle("Charge (Z)");
    fh1_Twimhit_z->GetYaxis()->SetTitle("Counts");
    fh1_Twimhit_z->GetYaxis()->SetTitleOffset(1.1);
    fh1_Twimhit_z->GetXaxis()->CenterTitle(true);
    fh1_Twimhit_z->GetYaxis()->CenterTitle(true);
    fh1_Twimhit_z->GetXaxis()->SetLabelSize(0.045);
    fh1_Twimhit_z->GetXaxis()->SetTitleSize(0.045);
    fh1_Twimhit_z->GetYaxis()->SetLabelSize(0.045);
    fh1_Twimhit_z->GetYaxis()->SetTitleSize(0.045);
    fh1_Twimhit_z->SetFillColor(29);
    fh1_Twimhit_z->SetLineColor(1);
    fh1_Twimhit_z->Draw("");

    TCanvas* cTwim_theta = new TCanvas("Twim_theta", "Twim: #theta_{XZ}", 10, 10, 800, 700);
    fh1_Twimhit_theta = R3B::root_owned<TH1F>("fh1_Twim_theta", "Twim: #theta_{XZ}", 900, -60, 60);
    fh1_Twimhit_theta->GetXaxis()->SetTitle("#theta_{XZ} [mrad]");
    fh1_Twimhit_theta->GetYaxis()->SetTitle("Counts");
    fh1_Twimhit_theta->GetYaxis()->SetTitleOffset(1.1);
    fh1_Twimhit_theta->GetXaxis()->CenterTitle(true);
    fh1_Twimhit_theta->GetYaxis()->CenterTitle(true);
    fh1_Twimhit_theta->GetXaxis()->SetLabelSize(0.045);
    fh1_Twimhit_theta->GetXaxis()->SetTitleSize(0.045);
    fh1_Twimhit_theta->GetYaxis()->SetLabelSize(0.045);
    fh1_Twimhit_theta->GetYaxis()->SetTitleSize(0.045);
    fh1_Twimhit_theta->SetFillColor(29);
    fh1_Twimhit_theta->SetLineColor(1);
    fh1_Twimhit_theta->Draw("");

    TCanvas* cTwim_zvstheta = new TCanvas("Twim_charge_vs_theta", "Twim: Charge Z vs #theta_{XZ}", 10, 10, 800, 700);
    fh2_Twimhit_zvstheta =
        R3B::root_owned<TH2F>("fh2_Twim_Charge_Z_vs_theta", "Twim: Charge Z vs #theta_{XZ}", 900, -40, 40, 1000, 6, 38);
    fh2_Twimhit_zvstheta->GetXaxis()->SetTitle("#theta_{XZ} [mrad]");
    fh2_Twimhit_zvstheta->GetYaxis()->SetTitle("Charge (Z)");
    fh2_Twimhit_zvstheta->GetYaxis()->SetTitleOffset(1.1);
    fh2_Twimhit_zvstheta->GetXaxis()->CenterTitle(true);
    fh2_Twimhit_zvstheta->GetYaxis()->CenterTitle(true);
    fh2_Twimhit_zvstheta->GetXaxis()->SetLabelSize(0.045);
    fh2_Twimhit_zvstheta->GetXaxis()->SetTitleSize(0.045);
    fh2_Twimhit_zvstheta->GetYaxis()->SetLabelSize(0.045);
    fh2_Twimhit_zvstheta->GetYaxis()->SetTitleSize(0.045);
    fh2_Twimhit_zvstheta->Draw("colz");

    cTwimTheta_vs_mwpc3x = new TCanvas("Twim_theta_vs_mwpc3x", "Twim: Theta vs Mwpc3-X", 10, 10, 800, 700);
    fh2_TwimTheta_vs_mwpc3x =
        R3B::root_owned<TH2F>("fh2_TwiwTheta_vs_mwpc3", "Twim: #theta_{xz} vs MWPC3X", 900, -40, 40, 800, -400, 400);
    fh2_TwimTheta_vs_mwpc3x->GetXaxis()->SetTitle("#theta_{XZ} [mrad]");
    fh2_TwimTheta_vs_mwpc3x->GetYaxis()->SetTitle("(Wixhausen)<---  Mwpc3-X [mm]  ---> (Messel)");
    fh2_TwimTheta_vs_mwpc3x->GetYaxis()->SetTitleOffset(1.1);
    fh2_TwimTheta_vs_mwpc3x->GetXaxis()->CenterTitle(true);
    fh2_TwimTheta_vs_mwpc3x->GetYaxis()->CenterTitle(true);
    fh2_TwimTheta_vs_mwpc3x->GetXaxis()->SetLabelSize(0.045);
    fh2_TwimTheta_vs_mwpc3x->GetXaxis()->SetTitleSize(0.045);
    fh2_TwimTheta_vs_mwpc3x->GetYaxis()->SetLabelSize(0.045);
    fh2_TwimTheta_vs_mwpc3x->GetYaxis()->SetTitleSize(0.045);
    fh2_TwimTheta_vs_mwpc3x->Draw("colz");

    cTwimZ_vs_mwpc3x = new TCanvas("Twim_Z_vs_mwpc3x", "Twim: Z vs Mwpc3-X", 10, 10, 800, 700);
    fh2_TwimZ_vs_mwpc3x = R3B::root_owned<TH2F>("fh2_TwiwZ_vs_mwpc3", "Twim: Z vs Mwpc3-X", 800, -400, 400, 700, 6, 38);
    fh2_TwimZ_vs_mwpc3x->GetXaxis()->SetTitle("(Wixhausen)<---  Mwpc3-X [mm]  ---> (Messel)");
    fh2_TwimZ_vs_mwpc3x->GetYaxis()->SetTitle("Charge Z");
    fh2_TwimZ_vs_mwpc3x->GetYaxis()->SetTitleOffset(1.1);
    fh2_TwimZ_vs_mwpc3x->GetXaxis()->CenterTitle(true);
    fh2_TwimZ_vs_mwpc3x->GetYaxis()->CenterTitle(true);
    fh2_TwimZ_vs_mwpc3x->GetXaxis()->SetLabelSize(0.045);
    fh2_TwimZ_vs_mwpc3x->GetXaxis()->SetTitleSize(0.045);
    fh2_TwimZ_vs_mwpc3x->GetYaxis()->SetLabelSize(0.045);
    fh2_TwimZ_vs_mwpc3x->GetYaxis()->SetTitleSize(0.045);
    fh2_TwimZ_vs_mwpc3x->Draw("colz");
    // s444 and s467 experiment: end of histograms ------

    //-----------------------------------------------------------
    // Hit data
    // s455 experiment: init of histograms ----------------------
    cTwimZs[0] = new TCanvas("Twim_ZL1_vs_ZR3", "Twim: Zs", 10, 10, 800, 700);
    cTwimZs[0]->Divide(2, 1);
    cTwimZs[0]->cd(1);
    fh1_Twimhit_Zl[0] = R3B::root_owned<TH1F>("fh1_Twim_zl1", "Twim: ZL1(blue) and ZR3(red)", 1200, 5, 100);
    fh1_Twimhit_Zl[0]->GetXaxis()->SetTitle("Z [atomic number]");
    fh1_Twimhit_Zl[0]->GetYaxis()->SetTitle("Counts");
    fh1_Twimhit_Zl[0]->GetYaxis()->SetTitleOffset(1.1);
    fh1_Twimhit_Zl[0]->GetXaxis()->CenterTitle(true);
    fh1_Twimhit_Zl[0]->GetYaxis()->CenterTitle(true);
    fh1_Twimhit_Zl[0]->GetXaxis()->SetLabelSize(0.045);
    fh1_Twimhit_Zl[0]->GetXaxis()->SetTitleSize(0.045);
    fh1_Twimhit_Zl[0]->GetYaxis()->SetLabelSize(0.045);
    fh1_Twimhit_Zl[0]->GetYaxis()->SetTitleSize(0.045);
    fh1_Twimhit_Zl[0]->SetLineColor(4);
    fh1_Twimhit_Zl[0]->Draw("");
    fh1_Twimhit_Zr[0] = R3B::root_owned<TH1F>("fh1_Twim_zr1", "Twim: ZR sec 3", 1200, 5, 100);
    fh1_Twimhit_Zr[0]->SetLineColor(2);
    fh1_Twimhit_Zr[0]->Draw("same");
    cTwimZs[0]->cd(2);
    fh2_Twimhit_ZrZl[0] =
        R3B::root_owned<TH2F>("fh2_Twim_zlzr1", "Twim: ZL sec 1 vs ZR sec 3", 1200, 5, 100, 1200, 5, 100);
    fh2_Twimhit_ZrZl[0]->GetXaxis()->SetTitle("Charge ZL");
    fh2_Twimhit_ZrZl[0]->GetYaxis()->SetTitle("Charge ZR");
    fh2_Twimhit_ZrZl[0]->GetYaxis()->SetTitleOffset(1.1);
    fh2_Twimhit_ZrZl[0]->GetXaxis()->CenterTitle(true);
    fh2_Twimhit_ZrZl[0]->GetYaxis()->CenterTitle(true);
    fh2_Twimhit_ZrZl[0]->GetXaxis()->SetLabelSize(0.045);
    fh2_Twimhit_ZrZl[0]->GetXaxis()->SetTitleSize(0.045);
    fh2_Twimhit_ZrZl[0]->GetYaxis()->SetLabelSize(0.045);
    fh2_Twimhit_ZrZl[0]->GetYaxis()->SetTitleSize(0.045);
    fh2_Twimhit_ZrZl[0]->Draw("colz");

    cTwimZs[1] = new TCanvas("Twim_ZL2_vs_ZR4", "Twim: Zs", 10, 10, 800, 700);
    cTwimZs[1]->Divide(2, 1);
    cTwimZs[1]->cd(1);
    fh1_Twimhit_Zl[1] = R3B::root_owned<TH1F>("fh1_Twim_zl2", "Twim: ZL2(blue) and ZR4(red)", 1200, 5, 100);
    fh1_Twimhit_Zl[1]->GetXaxis()->SetTitle("Z [atomic number]");
    fh1_Twimhit_Zl[1]->GetYaxis()->SetTitle("Counts");
    fh1_Twimhit_Zl[1]->GetYaxis()->SetTitleOffset(1.1);
    fh1_Twimhit_Zl[1]->GetXaxis()->CenterTitle(true);
    fh1_Twimhit_Zl[1]->GetYaxis()->CenterTitle(true);
    fh1_Twimhit_Zl[1]->GetXaxis()->SetLabelSize(0.045);
    fh1_Twimhit_Zl[1]->GetXaxis()->SetTitleSize(0.045);
    fh1_Twimhit_Zl[1]->GetYaxis()->SetLabelSize(0.045);
    fh1_Twimhit_Zl[1]->GetYaxis()->SetTitleSize(0.045);
    fh1_Twimhit_Zl[1]->SetLineColor(4);
    fh1_Twimhit_Zl[1]->Draw("");
    fh1_Twimhit_Zr[1] = R3B::root_owned<TH1F>("fh1_Twim_zr2", "Twim: ZR", 1200, 5, 100);
    fh1_Twimhit_Zr[1]->SetLineColor(2);
    fh1_Twimhit_Zr[1]->Draw("same");
    cTwimZs[1]->cd(2);
    fh2_Twimhit_ZrZl[1] =
        R3B::root_owned<TH2F>("fh2_Twim_zlzr2", "Twim: ZL sec 2 vs ZR sec 4", 1200, 5, 100, 1200, 5, 100);
    fh2_Twimhit_ZrZl[1]->GetXaxis()->SetTitle("Charge ZL");
    fh2_Twimhit_ZrZl[1]->GetYaxis()->SetTitle("Charge ZR");
    fh2_Twimhit_ZrZl[1]->GetYaxis()->SetTitleOffset(1.1);
    fh2_Twimhit_ZrZl[1]->GetXaxis()->CenterTitle(true);
    fh2_Twimhit_ZrZl[1]->GetYaxis()->CenterTitle(true);
    fh2_Twimhit_ZrZl[1]->GetXaxis()->SetLabelSize(0.045);
    fh2_Twimhit_ZrZl[1]->GetXaxis()->SetTitleSize(0.045);
    fh2_Twimhit_ZrZl[1]->GetYaxis()->SetLabelSize(0.045);
    fh2_Twimhit_ZrZl[1]->GetYaxis()->SetTitleSize(0.045);
    fh2_Twimhit_ZrZl[1]->Draw("colz");

    cTwimZsvsTof = new TCanvas("Twim_Zs_vs_Tof", "Twim: Zs vs Tof", 10, 10, 800, 700);
    cTwimZsvsTof->Divide(2, 2);
    cTwimZsvsTof->cd(1);
    fh2_Twimhit_ZlvsTof[0] =
        R3B::root_owned<TH2F>("fh1_Twim_zlvstof_sec1", "Twim: ZL vs ToF section 1", 1000, 29., 40., 1000, 5, 100);
    fh2_Twimhit_ZlvsTof[0]->GetYaxis()->SetTitle("Charge ZL");
    fh2_Twimhit_ZlvsTof[0]->GetXaxis()->SetTitle("ToF [ns]");
    fh2_Twimhit_ZlvsTof[0]->GetYaxis()->SetTitleOffset(1.1);
    fh2_Twimhit_ZlvsTof[0]->GetXaxis()->CenterTitle(true);
    fh2_Twimhit_ZlvsTof[0]->GetYaxis()->CenterTitle(true);
    fh2_Twimhit_ZlvsTof[0]->GetXaxis()->SetLabelSize(0.045);
    fh2_Twimhit_ZlvsTof[0]->GetXaxis()->SetTitleSize(0.045);
    fh2_Twimhit_ZlvsTof[0]->GetYaxis()->SetLabelSize(0.045);
    fh2_Twimhit_ZlvsTof[0]->GetYaxis()->SetTitleSize(0.045);
    fh2_Twimhit_ZlvsTof[0]->SetLineColor(4);
    fh2_Twimhit_ZlvsTof[0]->Draw("colz");
    cTwimZsvsTof->cd(2);
    fh2_Twimhit_ZlvsTof[1] =
        R3B::root_owned<TH2F>("fh1_Twim_zlvstof_sec2", "Twim: ZL vs ToF section 2", 1000, 29., 40., 1000, 5, 100);
    fh2_Twimhit_ZlvsTof[1]->GetYaxis()->SetTitle("Charge ZL");
    fh2_Twimhit_ZlvsTof[1]->GetXaxis()->SetTitle("ToF [ns]");
    fh2_Twimhit_ZlvsTof[1]->GetYaxis()->SetTitleOffset(1.1);
    fh2_Twimhit_ZlvsTof[1]->GetXaxis()->CenterTitle(true);
    fh2_Twimhit_ZlvsTof[1]->GetYaxis()->CenterTitle(true);
    fh2_Twimhit_ZlvsTof[1]->GetXaxis()->SetLabelSize(0.045);
    fh2_Twimhit_ZlvsTof[1]->GetXaxis()->SetTitleSize(0.045);
    fh2_Twimhit_ZlvsTof[1]->GetYaxis()->SetLabelSize(0.045);
    fh2_Twimhit_ZlvsTof[1]->GetYaxis()->SetTitleSize(0.045);
    fh2_Twimhit_ZlvsTof[1]->SetLineColor(4);
    fh2_Twimhit_ZlvsTof[1]->Draw("colz");
    cTwimZsvsTof->cd(3);
    fh2_Twimhit_ZrvsTof[0] =
        R3B::root_owned<TH2F>("fh2_Twim_zrvstof_sec3", "Twim: ZR vs ToF section 3", 1000, 29., 40., 1000, 5, 100);
    fh2_Twimhit_ZrvsTof[0]->GetYaxis()->SetTitle("Charge ZR");
    fh2_Twimhit_ZrvsTof[0]->GetXaxis()->SetTitle("ToF [ns]");
    fh2_Twimhit_ZrvsTof[0]->GetYaxis()->SetTitleOffset(1.1);
    fh2_Twimhit_ZrvsTof[0]->GetXaxis()->CenterTitle(true);
    fh2_Twimhit_ZrvsTof[0]->GetYaxis()->CenterTitle(true);
    fh2_Twimhit_ZrvsTof[0]->GetXaxis()->SetLabelSize(0.045);
    fh2_Twimhit_ZrvsTof[0]->GetXaxis()->SetTitleSize(0.045);
    fh2_Twimhit_ZrvsTof[0]->GetYaxis()->SetLabelSize(0.045);
    fh2_Twimhit_ZrvsTof[0]->GetYaxis()->SetTitleSize(0.045);
    fh2_Twimhit_ZrvsTof[0]->Draw("colz");
    cTwimZsvsTof->cd(4);
    fh2_Twimhit_ZrvsTof[1] =
        R3B::root_owned<TH2F>("fh2_Twim_zrvstof_sec4", "Twim: ZR vs ToF section 4", 1000, 29., 40., 1000, 5, 100);
    fh2_Twimhit_ZrvsTof[1]->GetYaxis()->SetTitle("Charge ZR");
    fh2_Twimhit_ZrvsTof[1]->GetXaxis()->SetTitle("ToF [ns]");
    fh2_Twimhit_ZrvsTof[1]->GetYaxis()->SetTitleOffset(1.1);
    fh2_Twimhit_ZrvsTof[1]->GetXaxis()->CenterTitle(true);
    fh2_Twimhit_ZrvsTof[1]->GetYaxis()->CenterTitle(true);
    fh2_Twimhit_ZrvsTof[1]->GetXaxis()->SetLabelSize(0.045);
    fh2_Twimhit_ZrvsTof[1]->GetXaxis()->SetTitleSize(0.045);
    fh2_Twimhit_ZrvsTof[1]->GetYaxis()->SetLabelSize(0.045);
    fh2_Twimhit_ZrvsTof[1]->GetYaxis()->SetTitleSize(0.045);
    fh2_Twimhit_ZrvsTof[1]->Draw("colz");

    for (Int_t i = 0; i < 4; i++)
    {
        sprintf(Name1, "ZL_vs_TofL%d", i + 1);
        sprintf(Name2, "Section %d", i + 1);
        cTwimZsvsTofL[i] = new TCanvas(Name1, Name2, 10, 10, 800, 700);
        cTwimZsvsTofL[i]->Divide(3, 2);
        for (Int_t j = 0; j < 6; j++)
        {
            sprintf(Name1, "fh2_ZL_vs_TofL%d", i * 6 + j + 1);
            sprintf(Name2, "ZL vs TofL for Sci %d", i * 6 + j + 1);
            fh2_Twimhit_ZlvsTofl[i * 6 + j] = R3B::root_owned<TH2F>(Name1, Name2, 1000, 29., 40., 1000, 5, 100);
            fh2_Twimhit_ZlvsTofl[i * 6 + j]->GetYaxis()->SetTitle("Charge ZL");
            fh2_Twimhit_ZlvsTofl[i * 6 + j]->GetXaxis()->SetTitle("ToF [ns]");
            fh2_Twimhit_ZlvsTofl[i * 6 + j]->GetYaxis()->SetTitleOffset(1.1);
            fh2_Twimhit_ZlvsTofl[i * 6 + j]->GetXaxis()->CenterTitle(true);
            fh2_Twimhit_ZlvsTofl[i * 6 + j]->GetYaxis()->CenterTitle(true);
            fh2_Twimhit_ZlvsTofl[i * 6 + j]->GetXaxis()->SetLabelSize(0.045);
            fh2_Twimhit_ZlvsTofl[i * 6 + j]->GetXaxis()->SetTitleSize(0.045);
            fh2_Twimhit_ZlvsTofl[i * 6 + j]->GetYaxis()->SetLabelSize(0.045);
            fh2_Twimhit_ZlvsTofl[i * 6 + j]->GetYaxis()->SetTitleSize(0.045);
            cTwimZsvsTofL[i]->cd(j + 1);
            fh2_Twimhit_ZlvsTofl[i * 6 + j]->Draw("colz");
        }

        sprintf(Name1, "ZR_vs_TofR%d", i + 1);
        sprintf(Name2, "Section %d", i + 1);
        cTwimZsvsTofR[i] = new TCanvas(Name1, Name2, 10, 10, 800, 700);
        cTwimZsvsTofR[i]->Divide(3, 2);
        for (Int_t j = 0; j < 6; j++)
        {
            sprintf(Name1, "fh2_ZR_vs_TofR%d", i * 6 + j + 5);
            sprintf(Name2, "ZR vs TofR for Sci %d", i * 6 + j + 5);
            fh2_Twimhit_ZrvsTofr[i * 6 + j] = R3B::root_owned<TH2F>(Name1, Name2, 1000, 29., 40., 1000, 5, 100);
            fh2_Twimhit_ZrvsTofr[i * 6 + j]->GetYaxis()->SetTitle("Charge ZR");
            fh2_Twimhit_ZrvsTofr[i * 6 + j]->GetXaxis()->SetTitle("ToF [ns]");
            fh2_Twimhit_ZrvsTofr[i * 6 + j]->GetYaxis()->SetTitleOffset(1.1);
            fh2_Twimhit_ZrvsTofr[i * 6 + j]->GetXaxis()->CenterTitle(true);
            fh2_Twimhit_ZrvsTofr[i * 6 + j]->GetYaxis()->CenterTitle(true);
            fh2_Twimhit_ZrvsTofr[i * 6 + j]->GetXaxis()->SetLabelSize(0.045);
            fh2_Twimhit_ZrvsTofr[i * 6 + j]->GetXaxis()->SetTitleSize(0.045);
            fh2_Twimhit_ZrvsTofr[i * 6 + j]->GetYaxis()->SetLabelSize(0.045);
            fh2_Twimhit_ZrvsTofr[i * 6 + j]->GetYaxis()->SetTitleSize(0.045);
            cTwimZsvsTofR[i]->cd(j + 1);
            fh2_Twimhit_ZrvsTofr[i * 6 + j]->Draw("colz");
        }
    }

    cTwimZsum = new TCanvas("Twim_Zsum", "Twim: ZL+ZR", 10, 10, 800, 700);
    fh1_twim_ZSum[0] = R3B::root_owned<TH1F>("fh1_Twim_ZSum13", "Twim: ZL1+ZR3(blue) and ZL2+ZR4(red)", 800, 60, 100);
    fh1_twim_ZSum[0]->GetXaxis()->SetTitle("Fissioning system --> Z [atomic number]");
    fh1_twim_ZSum[0]->GetYaxis()->SetTitle("Counts");
    fh1_twim_ZSum[0]->GetYaxis()->SetTitleOffset(1.1);
    fh1_twim_ZSum[0]->GetXaxis()->CenterTitle(true);
    fh1_twim_ZSum[0]->GetYaxis()->CenterTitle(true);
    fh1_twim_ZSum[0]->GetXaxis()->SetLabelSize(0.045);
    fh1_twim_ZSum[0]->GetXaxis()->SetTitleSize(0.045);
    fh1_twim_ZSum[0]->GetYaxis()->SetLabelSize(0.045);
    fh1_twim_ZSum[0]->GetYaxis()->SetTitleSize(0.045);
    fh1_twim_ZSum[0]->SetLineColor(4);
    fh1_twim_ZSum[0]->Draw("");
    fh1_twim_ZSum[1] = R3B::root_owned<TH1F>("fh1_Twim_ZSum24", "Twim: ZL+ZR", 800, 60, 100);
    fh1_twim_ZSum[1]->SetLineColor(2);
    fh1_twim_ZSum[1]->Draw("same");
    // s455 experiment: end of histograms -----------------------

    // MAIN FOLDER-Twim
    TFolder* mainfolTwim = new TFolder("TWIM", "TWIM info");

    // Folder for mapped data
    TFolder* mapfolTwim = new TFolder("Map", "Map TWIM info");
    for (Int_t i = 0; i < fNbSections; i++)
        mapfolTwim->Add(cTwimMap_E[i]);
    for (Int_t i = 0; i < fNbSections; i++)
        mapfolTwim->Add(cTwimMap_T[i]);
    for (Int_t i = 0; i < fNbSections; i++)
        mapfolTwim->Add(cTwimMult[i]);
    for (Int_t i = 0; i < fNbSections; i++)
        mapfolTwim->Add(cTwimMap_DT[i]);
    for (Int_t i = 0; i < fNbSections; i++)
        mapfolTwim->Add(cTwimMap_DeltaTrefTrig[i]);
    for (Int_t i = 0; i < fNbSections; i++)
        mapfolTwim->Add(cTwim_Mult[i]);
    if (fExpId == 91 || fExpId == 118 || fExpId == 444 || fExpId == 467 || fExpId == 509 || fExpId == 522)
    {
        mapfolTwim->Add(cTwimMap_ESum);
        mapfolTwim->Add(cTwimMap_ESum1);
        mapfolTwim->Add(cTwimMap_ESum2);
    }
    mapfolTwim->Add(cTwimMap_ESum_vs_diffDT);
    for (Int_t i = 0; i < fNbSections; i++)
    {
        mapfolTwim->Add(cTwimMap_EvsDT[i]);
        mapfolTwim->Add(cTwim_DTvsDT[i]);
    }
    mainfolTwim->Add(mapfolTwim);

    // Cal data --------------------------------------------------
    if (fCalItemsTwim)
    {
        // Folder for cal data
        TFolder* calfolTwim = new TFolder("Cal", "Cal TWIM info");
        for (Int_t i = 0; i < fNbSections; i++)
            calfolTwim->Add(cTwimCal_Pos[i]);

        mainfolTwim->Add(calfolTwim);
    }

    // Folder for hit data
    TFolder* hitfolTwim = new TFolder("Hit", "Hit TWIM info");

    // Hit data --------------------------------------------------
    if (fHitItemsTwim &&
        (fExpId == 91 || fExpId == 118 || fExpId == 444 || fExpId == 467 || fExpId == 509 || fExpId == 522))
    {
        hitfolTwim->Add(cTwim_Z);
        hitfolTwim->Add(cTwim_theta);
        hitfolTwim->Add(cTwim_zvstheta);
        if (fHitItemsMwpc3)
        {
            hitfolTwim->Add(cTwimTheta_vs_mwpc3x);
            hitfolTwim->Add(cTwimZ_vs_mwpc3x);
        }

        mainfolTwim->Add(hitfolTwim);
    }
    else if (fHitItemsTwim && fExpId == 455)
    {
        hitfolTwim->Add(cTwimZs[0]);
        hitfolTwim->Add(cTwimZs[1]);
        hitfolTwim->Add(cTwimZsum);
        hitfolTwim->Add(cTwimZsvsTof);
        for (Int_t i = 0; i < 4; i++)
            hitfolTwim->Add(cTwimZsvsTofL[i]);

        for (Int_t i = 0; i < 4; i++)
            hitfolTwim->Add(cTwimZsvsTofR[i]);

        mainfolTwim->Add(hitfolTwim);
    }

    run->AddObject(mainfolTwim);

    // Register command to reset histograms
    run->GetHttpServer()->RegisterCommand("Reset_TWIM_HIST", Form("/Objects/%s/->Reset_Histo()", GetName()));

    return kSUCCESS;
}

void R3BTwimOnlineSpectra::Reset_Histo()
{
    R3BLOG(info, "Reset_Histo");
    // Map data
    for (Int_t i = 0; i < fNbSections; i++)
    {
        fh1_Twimmap_mult[i]->Reset();
        fh2_twim_EneRawSumVsDriftTime[i]->Reset();
        fh2_twim_ESum_vs_diffDT[i]->Reset();
        for (Int_t j = 0; j < fNbAnodes; j++)
        {
            fh1_twimmap_E[i * fNbAnodes + j]->Reset();
            fh1_twimmap_T[i * fNbAnodes + j]->Reset();
            fh1_twimmap_DT[i * fNbAnodes + j]->Reset();
            fh2_twim_EneRawVsDriftTime[i * fNbAnodes + j]->Reset();
        }
        for (Int_t j = 0; j < fNbAnodes - 1; j++)
        {
            fh2_twim_DTvsDT[i * fNbAnodes + j]->Reset();
        }
        for (Int_t j = 0; j < fNbTref; j++)
        {
            fh1_twimmap_DeltaTrefTrig[i * fNbTref + j]->Reset();
        }
    }
    fh1_twim_ESum[0]->Reset();
    fh1_twim_ESum[1]->Reset();
    fh1_twim_ESum[2]->Reset();
    fh2_twim_ESum->Reset();

    // Cal data
    if (fCalItemsTwim)
    {
        for (Int_t i = 0; i < fNbSections; i++)
            for (Int_t j = 0; j < fNbAnodes; j++)
                fh1_Twimcal_Pos[i * fNbAnodes + j]->Reset();
    }

    // Hit data
    if (fHitItemsTwim && (fExpId == 444 || fExpId == 467 || fExpId == 509 || fExpId == 522))
    {
        fh1_Twimhit_z->Reset();
        fh1_Twimhit_theta->Reset();
        fh2_Twimhit_zvstheta->Reset();
        if (fHitItemsMwpc3)
        {
            fh2_TwimTheta_vs_mwpc3x->Reset();
            fh2_TwimZ_vs_mwpc3x->Reset();
        }
    }
    else if (fHitItemsTwim && fExpId == 455)
    {
        fh1_Twimhit_Zl[0]->Reset();
        fh1_Twimhit_Zr[0]->Reset();
        fh1_Twimhit_Zl[1]->Reset();
        fh1_Twimhit_Zr[1]->Reset();
        fh2_Twimhit_ZrZl[0]->Reset();
        fh2_Twimhit_ZrZl[1]->Reset();
        fh1_twim_ZSum[0]->Reset();
        fh1_twim_ZSum[1]->Reset();
        fh2_Twimhit_ZrvsTof[0]->Reset();
        fh2_Twimhit_ZlvsTof[1]->Reset();
        fh2_Twimhit_ZrvsTof[0]->Reset();
        fh2_Twimhit_ZlvsTof[1]->Reset();
        for (Int_t i = 0; i < 24; i++)
        {
            fh2_Twimhit_ZrvsTofr[i]->Reset();
            fh2_Twimhit_ZlvsTofl[i]->Reset();
        }
    }
}

void R3BTwimOnlineSpectra::Exec(Option_t* /*option*/)
{
    if (fExpId == 91 || fExpId == 118 || fExpId == 444 || fExpId == 467 || fExpId == 509 || fExpId == 522)
    {
        s444_s467();
    }
    else if (fExpId == 455)
    {
        s455();
    }
}

void R3BTwimOnlineSpectra::s455()
{
    // Fill mapped data
    if (fMappedItemsTwim && fMappedItemsTwim->GetEntriesFast() > 0)
    {
        Double_t Eraw[fNbSections][fNbAnodes];
        Double_t e1 = 0., e2 = 0.;
        Double_t Traw[fNbSections][fNbAnodes + fNbTref + fNbTrig];
        Double_t n1 = 0., n2 = 0.;
        UInt_t mult[fNbSections][fNbAnodes + fNbTref + fNbTrig];
        for (Int_t j = 0; j < fNbSections; j++)
            for (Int_t i = 0; i < fNbAnodes; i++)
                Eraw[j][i] = 0.; // mult=1 !!!
        for (Int_t j = 0; j < fNbSections; j++)
            for (Int_t i = 0; i < fNbAnodes + fNbTref + fNbTrig; i++)
            {
                Traw[j][i] = 0.; // mult=1 !!!
                mult[j][i] = 0;
            }

        Int_t nHits = fMappedItemsTwim->GetEntriesFast();
        // std::cout << "Event:\n";
        for (Int_t ihit = 0; ihit < nHits; ihit++)
        {
            R3BTwimMappedData* hit = dynamic_cast<R3BTwimMappedData*>(fMappedItemsTwim->At(ihit));
            if (!hit)
                continue;
            fh1_Twimmap_mult[hit->GetSecID() - 1]->Fill(hit->GetAnodeID());
            mult[hit->GetSecID() - 1][hit->GetAnodeID() - 1]++;
            if (Traw[hit->GetSecID() - 1][hit->GetAnodeID() - 1] == 0)
            {
                if (hit->GetEnergy() < 65535 && hit->GetEnergy() > 0 &&
                    Eraw[hit->GetSecID() - 1][hit->GetAnodeID() - 1] == 0)
                    Eraw[hit->GetSecID() - 1][hit->GetAnodeID() - 1] = hit->GetEnergy(); // mult=1 !!!
                Traw[hit->GetSecID() - 1][hit->GetAnodeID() - 1] = hit->GetTime();       // mult=1 !!!
            }
        } // end of loop over the Mapped data

        Int_t idTref = fNbAnodes;
        for (Int_t j = 0; j < fNbSections; j++)
        {
            e1 = 0.;
            n1 = 0.;
            for (Int_t i = 0; i < fNbTref; i++)
            {
                if (mult[j][16 + i] == 1 && mult[j][17 + i] == 1) // FIXME
                {
                    fh1_twimmap_DeltaTrefTrig[j * fNbTref + i]->Fill(Traw[j][16 + i] - Traw[j][17 + i]);
                }
            }
            for (Int_t i = 0; i < fNbAnodes; i++)
            {
                fh1_twimmult[j * fNbAnodes + i]->Fill(mult[j][i]);
                if ((mult[j][i] == 1) && (mult[j][fNbAnodes + fNbTref] == 1))
                {
                    // to check the searching window of the MDPP16 fNbAnodes+fNbTref = idTtri=g
                    fh1_twimmap_T[j * fNbAnodes + i]->Fill(Traw[j][i] - Traw[j][fNbAnodes + fNbTref]);
                }
                if ((mult[j][i] == 1) && (mult[j][idTref] == 1))
                {
                    fh1_twimmap_E[j * fNbAnodes + i]->Fill(Eraw[j][i]);
                    fh1_twimmap_DT[j * fNbAnodes + i]->Fill(Traw[j][i] - Traw[j][idTref]);
                    fh2_twim_EneRawVsDriftTime[j * fNbAnodes + i]->Fill(Eraw[j][i], Traw[j][i] - Traw[j][idTref]);

                    e1 = e1 + Eraw[j][i];
                    n1++;
                }
            }
            for (Int_t i = 0; i < fNbAnodes - 1; i++)
            {
                // Tref = 16 for ch 0 to 15
                if ((mult[j][i] == 1) && (mult[j][i + 1] == 1) && (mult[j][idTref] == 1))
                {
                    fh2_twim_DTvsDT[j * fNbAnodes + i]->Fill(
                        Traw[j][i] - Traw[j][idTref],
                        Traw[j][i] - Traw[j][idTref] - (Traw[j][i + 1] - Traw[j][idTref]));
                }
            }
            if ((mult[j][15] == 1) && (mult[j][0] == 1) && (mult[j][16] == 1))
            {
                fh2_twim_ESum_vs_diffDT[j]->Fill(Traw[j][15] - Traw[j][0], e1 / n1);
            }
            if ((mult[j][5] == 1) && (mult[j][16] == 1))
            {
                fh2_twim_EneRawSumVsDriftTime[j]->Fill(Traw[j][5] - Traw[j][16], e1 / n1);
            }
        }
    }

    // Fill cal data
    if (fCalItemsTwim && fCalItemsTwim->GetEntriesFast() > 0)
    {
        Int_t nHits = fCalItemsTwim->GetEntriesFast();
        for (Int_t ihit = 0; ihit < nHits; ihit++)
        {
            R3BTwimCalData* hit = dynamic_cast<R3BTwimCalData*>(fCalItemsTwim->At(ihit));
            if (!hit)
                continue;
            fh1_Twimcal_Pos[(hit->GetSecID() - 1) * fNbAnodes + (hit->GetAnodeID() - 1)]->Fill(hit->GetDTime());
        }
    }

    // Fill TofW Hit data
    Double_t tof[2] = { 0., 0. };
    int padid[2] = { 0, 0 };
#ifdef DSOFIADEFINED
    if (fHitItemsTofW && fHitItemsTofW->GetEntriesFast() > 0)
    {
        Int_t nHits = fHitItemsTofW->GetEntriesFast();
        for (Int_t ihit = 0; ihit < nHits; ihit++)
        {
            R3BSofTofWHitData* hit = (R3BSofTofWHitData*)fHitItemsTofW->At(ihit);
            if (!hit)
                continue;
            if (padid[0] == 0)
            {
                padid[0] = hit->GetPaddle();
                tof[0] = hit->GetTof();
            }
            else
            {
                if (hit->GetPaddle() > padid[0] && hit->GetPaddle() - padid[0] > 1)
                {
                    padid[1] = hit->GetPaddle();
                    tof[1] = hit->GetTof(); // right
                }
                else if (hit->GetPaddle() - padid[0] < -1)
                {
                    tof[1] = tof[0]; // right
                    padid[1] = padid[0];
                    tof[0] = hit->GetTof(); // new left
                    padid[0] = hit->GetPaddle();
                }
            }
        }
    }
#endif

    // Fill hit data
    if (fHitItemsTwim && fHitItemsTwim->GetEntriesFast() > 0)
    {
        Int_t nHits = fHitItemsTwim->GetEntriesFast();
        Float_t zr[2] = { 0., 0. }, zl[2] = { 0., 0. };
        for (Int_t ihit = 0; ihit < nHits; ihit++)
        {
            R3BTwimHitData* hit = dynamic_cast<R3BTwimHitData*>(fHitItemsTwim->At(ihit));
            if (!hit)
                continue;
            // FIXME: this is defined only for the experiment 4-march-2021
            if (hit->GetSecID() == 1 || hit->GetSecID() == 2)
                zl[hit->GetSecID() - 1] = hit->GetZcharge();
            else if (hit->GetSecID() == 3 || hit->GetSecID() == 4)
                zr[hit->GetSecID() - 3] = hit->GetZcharge();
        }
        // FIXME in the future for complex fission events!
        if (zl[0] > 0. && tof[0] > 0.)
        {
            fh2_Twimhit_ZlvsTof[0]->Fill(tof[0], zl[0]);
            if (padid[0] < 24)
                fh2_Twimhit_ZlvsTofl[padid[0] - 1]->Fill(tof[0], zl[0]);
        }
        else if (zl[1] > 0. && tof[0] > 0.)
        {
            fh2_Twimhit_ZlvsTof[1]->Fill(tof[0], zl[1]);
            if (padid[0] < 24)
                fh2_Twimhit_ZlvsTofl[padid[0] - 1]->Fill(tof[0], zl[1]);
        }

        if (zr[0] > 0. && tof[1] > 0.)
        {
            fh2_Twimhit_ZrvsTof[0]->Fill(tof[1], zr[0]);
            if (padid[1] > 4)
                fh2_Twimhit_ZrvsTofr[padid[1] - 5]->Fill(tof[1], zr[0]);
        }
        else if (zr[1] > 0. && tof[1] > 0.)
        {

            fh2_Twimhit_ZrvsTof[1]->Fill(tof[1], zr[1]);
            if (padid[1] > 4)
                fh2_Twimhit_ZrvsTofr[padid[1] - 5]->Fill(tof[1], zr[1]);
        }

        if (zl[0] > 0.)
            fh1_Twimhit_Zl[0]->Fill(zl[0]);
        else if (zl[1] > 0.)
            fh1_Twimhit_Zl[1]->Fill(zl[1]);

        if (zr[0] > 0.)
            fh1_Twimhit_Zr[0]->Fill(zr[0]);
        else if (zr[1] > 0.)
            fh1_Twimhit_Zr[1]->Fill(zr[1]);

        if (zr[0] > 0. && zl[0] > 0. /*&& tof[0] > 0. && tof[1] > 0.*/)
        {
            fh2_Twimhit_ZrZl[0]->Fill(zl[0], zr[0]);
            fh1_twim_ZSum[0]->Fill(zl[0] + zr[0]);
        }
        else if (zr[1] > 0. && zl[1] > 0. /*&& tof[0] > 0. && tof[1] > 0.*/)
        {
            fh2_Twimhit_ZrZl[1]->Fill(zl[1], zr[1]);
            fh1_twim_ZSum[1]->Fill(zl[1] + zr[1]);
        }
    }
}

void R3BTwimOnlineSpectra::s444_s467()
{
    // Fill mapped data
    if (fMappedItemsTwim && fMappedItemsTwim->GetEntriesFast() > 0)
    {
        Double_t Eraw[fNbSections][fNbAnodes];
        Double_t e1 = 0., e2 = 0.;
        Double_t Traw[fNbSections][fNbAnodes + fNbTref + fNbTrig];
        Double_t n1 = 0., n2 = 0.;
        UInt_t mult[fNbSections][fNbAnodes + fNbTref + fNbTrig];
        for (Int_t j = 0; j < fNbSections; j++)
            for (Int_t i = 0; i < fNbAnodes; i++)
                Eraw[j][i] = 0.; // mult=1 !!!
        for (Int_t j = 0; j < fNbSections; j++)
            for (Int_t i = 0; i < fNbAnodes + fNbTref + fNbTrig; i++)
            {
                Traw[j][i] = 0.; // mult=1 !!!
                mult[j][i] = 0;
            }
        Int_t nHits = fMappedItemsTwim->GetEntriesFast();
        // std::cout << "Event:\n";
        for (Int_t ihit = 0; ihit < nHits; ihit++)
        {
            R3BTwimMappedData* hit = dynamic_cast<R3BTwimMappedData*>(fMappedItemsTwim->At(ihit));
            if (!hit)
                continue;
            fh1_Twimmap_mult[hit->GetSecID() - 1]->Fill(hit->GetAnodeID());
            mult[hit->GetSecID() - 1][hit->GetAnodeID() - 1]++;
            if (Traw[hit->GetSecID() - 1][hit->GetAnodeID() - 1] == 0)
            {
                if (hit->GetEnergy() < 8192 && hit->GetEnergy() > 0 &&
                    Eraw[hit->GetSecID() - 1][hit->GetAnodeID() - 1] == 0)
                    Eraw[hit->GetSecID() - 1][hit->GetAnodeID() - 1] = hit->GetEnergy(); // mult=1 !!!
                Traw[hit->GetSecID() - 1][hit->GetAnodeID() - 1] = hit->GetTime();       // mult=1 !!!
            }
        }
        Int_t idTref;
        for (Int_t j = 0; j < fNbSections; j++)
        {
            for (Int_t i = 0; i < fNbTref; i++)
            {
                if (mult[j][16 + i] == 1 && mult[j][18 + i] == 1) // NOLINT
                {
                    fh1_twimmap_DeltaTrefTrig[j * fNbTref + i]->Fill(Traw[j][16 + i] - Traw[j][18 + i]);
                }
            }
            for (Int_t i = 0; i < fNbAnodes; i++)
            {
                // Tref = 16 for ch 0 to 7 and Tref = 17 for ch 8 to 15
                idTref = fNbAnodes + i / 8;
                if ((mult[j][i] == 1) && (mult[j][idTref] == 1)) // NOLINT
                {
                    fh1_twimmap_E[j * fNbAnodes + i]->Fill(Eraw[j][i]);
                    fh1_twimmap_DT[j * fNbAnodes + i]->Fill(Traw[j][i] - Traw[j][idTref]);
                    fh2_twim_EneRawVsDriftTime[j * fNbAnodes + i]->Fill(Eraw[j][i], Traw[j][i] - Traw[j][idTref]);

                    if (i < fNbAnodes / 2)
                    {
                        e1 = e1 + Eraw[j][i];
                        n1++;
                    }
                    else if (i >= fNbAnodes / 2 && i < fNbAnodes)
                    {
                        e2 = e2 + Eraw[j][i];
                        n2++;
                    }
                }
            }
            for (Int_t i = 0; i < fNbAnodes - 1; i++)
            {
                idTref = fNbAnodes + i / 8;
                if ((mult[j][i] == 1) && (mult[j][i + 1] == 1) && (mult[j][idTref] == 1))
                {
                    fh2_twim_DTvsDT[j * fNbAnodes + i]->Fill(Traw[j][i] - Traw[j][idTref],
                                                             Traw[j][i + 1] - Traw[j][idTref]);
                }
            }
            if ((mult[j][15] == 1) && (mult[j][0] == 1) && (mult[j][16] == 1) && (mult[j][17] == 1)) // NOLINT
            {
                fh2_twim_ESum_vs_diffDT[j]->Fill((Traw[j][15] - Traw[j][17]) - (Traw[j][0] - Traw[j][16]),
                                                 (e1 + e2) / (n1 + n2));
            }
            if ((mult[j][5] == 1) && (mult[j][16] == 1)) // NOLINT
            {
                fh2_twim_EneRawSumVsDriftTime[j]->Fill(Traw[j][5] - Traw[j][16], (e1 + e2) / (n1 + n2));
            }
        }
        fh1_twim_ESum[0]->Fill(e1 / n1);
        fh1_twim_ESum[1]->Fill(e2 / n2);
        fh1_twim_ESum[2]->Fill((e1 + e2) / (n1 + n2));
        fh2_twim_ESum->Fill(e1 / n1, e2 / n2);
    }

    // Fill cal data
    if (fCalItemsTwim && fCalItemsTwim->GetEntriesFast() > 0)
    {
        Int_t nHits = fCalItemsTwim->GetEntriesFast();
        for (Int_t ihit = 0; ihit < nHits; ihit++)
        {
            R3BTwimCalData* hit = dynamic_cast<R3BTwimCalData*>(fCalItemsTwim->At(ihit));
            if (!hit)
                continue;
            fh1_Twimcal_Pos[(hit->GetSecID() - 1) * fNbAnodes + (hit->GetAnodeID() - 1)]->Fill(hit->GetDTime());
        }
    }

    Double_t mwpc3x = -1500;
    // Fill MWPC3 Hit data
    if (fHitItemsMwpc3 && fHitItemsMwpc3->GetEntriesFast() > 0)
    {
        Int_t nHits = fHitItemsMwpc3->GetEntriesFast();
        for (Int_t ihit = 0; ihit < nHits; ihit++)
        {
            R3BMwpcHitData* hit = dynamic_cast<R3BMwpcHitData*>(fHitItemsMwpc3->At(ihit));
            if (!hit)
                continue;
            mwpc3x = hit->GetX();
        }
    }

    // Fill hit data
    if (fHitItemsTwim && fHitItemsTwim->GetEntriesFast() > 0 && mwpc3x > -1500)
    {
        Int_t nHits = fHitItemsTwim->GetEntriesFast();
        for (Int_t ihit = 0; ihit < nHits; ihit++)
        {
            R3BTwimHitData* hit = dynamic_cast<R3BTwimHitData*>(fHitItemsTwim->At(ihit));
            if (!hit)
                continue;
            fh1_Twimhit_z->Fill(hit->GetZcharge());
            fh1_Twimhit_theta->Fill(hit->GetTheta() * 1000.);
            fh2_Twimhit_zvstheta->Fill(hit->GetTheta() * 1000., hit->GetZcharge());
            if (mwpc3x > -500)
            {
                fh2_TwimTheta_vs_mwpc3x->Fill(hit->GetTheta() * 1000., mwpc3x);
                fh2_TwimZ_vs_mwpc3x->Fill(mwpc3x, hit->GetZcharge());
            }
        }
    }

    fNEvents += 1;
}

void R3BTwimOnlineSpectra::FinishEvent()
{
    if (fMappedItemsTwim)
    {
        fMappedItemsTwim->Clear();
    }
    if (fCalItemsTwim)
    {
        fCalItemsTwim->Clear();
    }
    if (fHitItemsTwim)
    {
        fHitItemsTwim->Clear();
    }
    if (fHitItemsMwpc3)
    {
        fHitItemsMwpc3->Clear();
    }
    if (fHitItemsTofW)
    {
        fHitItemsTofW->Clear();
    }
}

void R3BTwimOnlineSpectra::FinishTask()
{
    if (fMappedItemsTwim)
    {
        for (Int_t i = 0; i < fNbSections; i++)
        {
            fh1_Twimmap_mult[i]->Write();
            cTwimMap_E[i]->Write();
            cTwimMap_T[i]->Write();
            cTwimMult[i]->Write();
            cTwimMap_DT[i]->Write();
            cTwimMap_EvsDT[i]->Write();
            cTwim_DTvsDT[i]->Write();
            cTwimMap_DeltaTrefTrig[i]->Write();
        }
        cTwimMap_EsumvsDT->Write();
        cTwimMap_ESum_vs_diffDT->Write();
        fh1_twim_ESum[0]->Write();
        fh1_twim_ESum[1]->Write();
        fh1_twim_ESum[2]->Write();
        fh2_twim_ESum->Write();
    }
    if (fCalItemsTwim)
    {
        for (Int_t i = 0; i < fNbSections; i++)
            cTwimCal_Pos[i]->Write();
    }
    if (fHitItemsTwim)
    {
        if (fExpId == 91 || fExpId == 118 || fExpId == 444 || fExpId == 467 || fExpId == 509 || fExpId == 522)
        {
            fh1_Twimhit_z->Write();
            fh1_Twimhit_theta->Write();
            fh2_Twimhit_zvstheta->Write();
            if (fHitItemsMwpc3)
            {
                cTwimTheta_vs_mwpc3x->Write();
                cTwimZ_vs_mwpc3x->Write();
            }
        }
        else if (fExpId == 455)
        {
            cTwimZs[0]->Write();
            cTwimZs[1]->Write();
            cTwimZsum->Write();
            cTwimZsvsTof->Write();
            for (Int_t i = 0; i < 4; i++)
            {
                cTwimZsvsTofL[i]->Write();
                cTwimZsvsTofR[i]->Write();
            }
        }
    }
}

ClassImp(R3BTwimOnlineSpectra)
