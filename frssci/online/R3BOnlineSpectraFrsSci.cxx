/******************************************************************************
 *   Copyright (C) 2019 GSI Helmholtzzentrum für Schwerionenforschung GmbH    *
 *   Copyright (C) 2019-2024 Members of R3B Collaboration                     *
 *                                                                            *
 *             This software is distributed under the terms of the            *
 *                 GNU General Public Licence (GPL) version 3,                *
 *                    copied verbatim in the file "LICENSE".                  *
 *                                                                            *
 * In applying this license GSI does not waive the privileges and immunities  *
 * granted to it by virtue of its status as an Intergovernmental Organization *
 * or submit itself to any jurisdiction.                                      *
 ******************************************************************************/

#include "FairLogger.h"
#include "FairRootManager.h"
#include "FairRunOnline.h"
#include "FairRuntimeDb.h"
#include "TCanvas.h"
#include "TClonesArray.h"
#include "TFolder.h"
#include "TH1.h"
#include "TH2.h"
#include "THttpServer.h"
#include "TMath.h"

#include "R3BEventHeader.h"
#include "R3BFrsSciMappedData.h"
#include "R3BFrsSciPosCalData.h"
#include "R3BFrsSciTcalData.h"
#include "R3BFrsSciTofCalData.h"
#include "R3BOnlineSpectraFrsSci.h"

R3BOnlineSpectraFrsSci::R3BOnlineSpectraFrsSci()
    : R3BOnlineSpectraFrsSci("FrsSciOnlineSpectra", 1)
{
}

R3BOnlineSpectraFrsSci::R3BOnlineSpectraFrsSci(const char* name, Int_t iVerbose)
    : FairTask(name, iVerbose)
    , fMapped(NULL)
    , fTcal(NULL)
    , fPosCal(NULL)
    , fTofCal(NULL)
    , fNEvents(0)
    , fNbDets(3)
    , fNbPmts(3)
    , fNbTofs(3)
    , fpos_range_min(-5.)
    , fpos_range_max(5.)
    , ftof_range_min(0.)
    , ftof_range_max(2500)
{
}

R3BOnlineSpectraFrsSci::~R3BOnlineSpectraFrsSci()
{
    LOG(debug) << "R3BOnlineSpectraFrsSci::Destructor";
    if (fMapped)
        delete fMapped;
    if (fTcal)
        delete fTcal;
    if (fPosCal)
        delete fPosCal;
    if (fTofCal)
        delete fTofCal;
}

InitStatus R3BOnlineSpectraFrsSci::Init()
{
    LOG(info) << "R3BOnlineSpectraFrsSci::Init()";

    char Name1[255];
    char Name2[255];

    // try to get a handle on the EventHeader. EventHeader may not be
    // present though and hence may be null. Take care when using.

    FairRootManager* mgr = FairRootManager::Instance();
    if (NULL == mgr)
        LOG(fatal) << "R3BOnlineSpectraFrsSci::Init FairRootManager not found";

    FairRunOnline* run = FairRunOnline::Instance();
    run->GetHttpServer()->Register("", this);

    // --- ------------ --- //
    // --- MAPPED LEVEL --- //
    // --- ------------ --- //

    // === get access to mapped data ===//
    fMapped = dynamic_cast<TClonesArray*>(mgr->GetObject("FrsSciMappedData"));
    if (!fMapped)
    {
        LOG(fatal) << "FrsSciMapped not found";
        return (kFATAL);
    }

    // === declare and create TCanvas and Histogram === //

    sprintf(Name1, "MapFinetime");
    cMap_FT = new TCanvas(Name1, Name1, 10, 10, 800, 700);
    cMap_FT->Divide(fNbPmts, fNbDets);
    fh1_Map_finetime = new TH1I*[fNbDets * fNbPmts];

    sprintf(Name1, "MapMult1D");
    cMap_Mult1D = new TCanvas(Name1, Name1, 10, 10, 800, 700);
    cMap_Mult1D->Divide(fNbPmts, fNbDets);
    fh1_Map_mult = new TH1I*[fNbDets * fNbPmts];

    sprintf(Name1, "MapMult2D");
    cMap_Mult2D = new TCanvas(Name1, Name1, 10, 10, 800, 700);
    cMap_Mult2D->Divide(2, fNbDets);
    fh2_Map_mult = new TH2I*[fNbDets];
    fh2_Map_multRvsL = new TH2I*[fNbDets];

    for (UShort_t i = 0; i < fNbDets; i++)
    {
        for (UShort_t j = 0; j < fNbPmts; j++)
        {
            // === TH1F: fine time === //
            sprintf(Name1, "FrsSci%i_FineTimePmt%i", i + 1, j + 1);
            fh1_Map_finetime[i * fNbPmts + j] = new TH1I(Name1, Name1, 1000, 0, 1000);
            fh1_Map_finetime[i * fNbPmts + j]->GetXaxis()->SetTitle("Fine time");
            fh1_Map_finetime[i * fNbPmts + j]->GetYaxis()->SetTitle("Counts per bin");
            fh1_Map_finetime[i * fNbPmts + j]->GetXaxis()->CenterTitle(true);
            fh1_Map_finetime[i * fNbPmts + j]->GetYaxis()->CenterTitle(true);
            fh1_Map_finetime[i * fNbPmts + j]->GetXaxis()->SetLabelSize(0.05);
            fh1_Map_finetime[i * fNbPmts + j]->GetXaxis()->SetTitleSize(0.05);
            fh1_Map_finetime[i * fNbPmts + j]->GetYaxis()->SetLabelSize(0.05);
            fh1_Map_finetime[i * fNbPmts + j]->GetYaxis()->SetTitleSize(0.05);
            cMap_FT->cd(i * fNbPmts + j + 1);
            fh1_Map_finetime[i * fNbPmts + j]->Draw("");

            // === TH1I: 1D-mult at map level === //
            sprintf(Name1, "FrsSci%i_Pmt%i_MultPerEvent_Mapped", i + 1, j + 1);
            sprintf(Name2,
                    "FrsSci%i_Pmt%i_MultPerEvent_Mapped (blue no condition on TPAT, red condition on TPAT = 1 or 2)",
                    i + 1,
                    j + 1);
            fh1_Map_mult[i * fNbPmts + j] = new TH1I(Name1, Name1, 70, -0.5, 69.5);
            fh1_Map_mult[i * fNbPmts + j]->GetXaxis()->SetTitle("Multiplicity per event");
            fh1_Map_mult[i * fNbPmts + j]->GetYaxis()->SetTitle("Counts");
            fh1_Map_mult[i * fNbPmts + j]->GetXaxis()->CenterTitle(true);
            fh1_Map_mult[i * fNbPmts + j]->GetYaxis()->CenterTitle(true);
            fh1_Map_mult[i * fNbPmts + j]->GetXaxis()->SetLabelSize(0.05);
            fh1_Map_mult[i * fNbPmts + j]->GetXaxis()->SetTitleSize(0.05);
            fh1_Map_mult[i * fNbPmts + j]->GetYaxis()->SetLabelSize(0.05);
            fh1_Map_mult[i * fNbPmts + j]->GetYaxis()->SetTitleSize(0.05);
            fh1_Map_mult[i * fNbPmts + j]->SetLineWidth(2);
            fh1_Map_mult[i * fNbPmts + j]->SetLineColor(kBlue);
            cMap_Mult1D->cd(i * fNbPmts + j + 1);
            gPad->SetLogy();
            fh1_Map_mult[i * fNbPmts + j]->Draw("");
        }

        // === TH1F: multiplicity per event and channel at mapped level === //
        sprintf(Name1, "FrsSci%i_MultPerEvent", i + 1);
        fh2_Map_mult[i] = new TH2I(Name1, Name1, fNbPmts + 1, -0.5, fNbPmts + 0.5, 20, -0.5, 19.5);
        fh2_Map_mult[i]->GetXaxis()->SetTitle("channel: 1=PMT R,    2=PMT L,    3=COMMON REF");
        fh2_Map_mult[i]->GetYaxis()->SetTitle("multiplicity per event");
        fh2_Map_mult[i]->GetXaxis()->CenterTitle(true);
        fh2_Map_mult[i]->GetYaxis()->CenterTitle(true);
        fh2_Map_mult[i]->GetXaxis()->SetLabelSize(0.05);
        fh2_Map_mult[i]->GetXaxis()->SetTitleSize(0.05);
        fh2_Map_mult[i]->GetYaxis()->SetLabelSize(0.05);
        fh2_Map_mult[i]->GetYaxis()->SetTitleSize(0.05);
        cMap_Mult2D->cd(i * 2 + 1);
        fh2_Map_mult[i]->Draw("COL");

        // === TH1F: multiplicity per event and channel at mapped level === //
        sprintf(Name1, "FrsSci%i_MultPerEvent_RvsL", i + 1);
        fh2_Map_multRvsL[i] = new TH2I(Name1, Name1, 40, -1.5, 38.5, 40, -1.5, 38.5);
        fh2_Map_multRvsL[i]->GetXaxis()->SetTitle("Multiplicity per event on the Left Pmt");
        fh2_Map_multRvsL[i]->GetYaxis()->SetTitle("Multiplicity per event on the Right Pmt");
        fh2_Map_multRvsL[i]->GetXaxis()->CenterTitle(true);
        fh2_Map_multRvsL[i]->GetYaxis()->CenterTitle(true);
        fh2_Map_multRvsL[i]->GetXaxis()->SetLabelSize(0.05);
        fh2_Map_multRvsL[i]->GetXaxis()->SetTitleSize(0.05);
        fh2_Map_multRvsL[i]->GetYaxis()->SetLabelSize(0.05);
        fh2_Map_multRvsL[i]->GetYaxis()->SetTitleSize(0.05);
        cMap_Mult2D->cd(i * 2 + 2);
        fh2_Map_multRvsL[i]->Draw("COL");

    } // end of loop over fNbDets

    // --- ------------ --- //
    // --- TCAL LEVEL --- //
    // --- ------------ --- //

    // === get access to tcal data ===//
    fTcal = dynamic_cast<TClonesArray*>(mgr->GetObject("FrsSciTcalData"));
    if (!fTcal)
    {
        LOG(info) << "R3BOnlineSpectraFrsSci::Init() :: FrsSciTcal not found";
        LOG(info) << "    No FrsSciTcal spectra ---> ok";
    }
    else
    {
        // === declare and create TCanvas and Histogram === //

        sprintf(Name1, "Tcal1Hit_PosRaw");
        cTcal_Pos = new TCanvas(Name1, Name1, 10, 10, 800, 700);
        cTcal_Pos->Divide(1, fNbDets);
        fh1_Tcal1Hit_PosRaw = new TH1D*[fNbDets];

        sprintf(Name1, "Tcal1Hit_TofRaw");
        cTcal_Tof = new TCanvas(Name1, Name1, 10, 10, 800, 700);
        cTcal_Tof->Divide(1, fNbTofs);
        fh1_Tcal1Hit_TofRaw = new TH1D*[fNbTofs];

        for (UShort_t i = 0; i < fNbDets; i++)
        {
            // === TH1F: Raw Position in Ns if mult1 RIGHT and LEFT === //
            sprintf(Name1, "FrsSci%i_PosRaw_MULT1", i + 1);
            fh1_Tcal1Hit_PosRaw[i] = new TH1D(Name1,
                                              Name1,
                                              static_cast<int>(20. * (fpos_range_max - fpos_range_min)),
                                              fpos_range_min,
                                              fpos_range_max);
            fh1_Tcal1Hit_PosRaw[i]->GetXaxis()->SetTitle("Raw Positon [ns] if mult1 at L and R. 50 ps per bin");
            fh1_Tcal1Hit_PosRaw[i]->GetYaxis()->SetTitle("number of counts with mult1");
            fh1_Tcal1Hit_PosRaw[i]->GetXaxis()->CenterTitle(true);
            fh1_Tcal1Hit_PosRaw[i]->GetYaxis()->CenterTitle(true);
            fh1_Tcal1Hit_PosRaw[i]->GetXaxis()->SetLabelSize(0.05);
            fh1_Tcal1Hit_PosRaw[i]->GetXaxis()->SetTitleSize(0.05);
            fh1_Tcal1Hit_PosRaw[i]->GetYaxis()->SetLabelSize(0.05);
            fh1_Tcal1Hit_PosRaw[i]->GetYaxis()->SetTitleSize(0.05);
            cTcal_Pos->cd(i + 1);
            fh1_Tcal1Hit_PosRaw[i]->Draw();
        }

        if (fNbTofs > 0)
        {
            UShort_t cpt = 0;
            for (UShort_t sta = 1; sta < fNbDets; sta++)
            {
                for (UShort_t sto = sta + 1; sto <= fNbDets; sto++)
                {
                    sprintf(Name1, "TofRaw_FrsSci%i_to_FrsSci%i_MULT1", sta, sto);
                    fh1_Tcal1Hit_TofRaw[cpt] = new TH1D(Name1,
                                                        Name1,
                                                        static_cast<int>(10. * (ftof_range_max - ftof_range_min)),
                                                        ftof_range_min,
                                                        fpos_range_max);
                    fh1_Tcal1Hit_TofRaw[cpt]->GetXaxis()->SetTitle(
                        "Raw Tof [ns] if mult1 at L, R and Tref. 100 ps per bin.");
                    fh1_Tcal1Hit_TofRaw[cpt]->GetYaxis()->SetTitle("number of counts with mult1");
                    fh1_Tcal1Hit_TofRaw[cpt]->GetXaxis()->CenterTitle(true);
                    fh1_Tcal1Hit_TofRaw[cpt]->GetYaxis()->CenterTitle(true);
                    fh1_Tcal1Hit_TofRaw[cpt]->GetXaxis()->SetLabelSize(0.05);
                    fh1_Tcal1Hit_TofRaw[cpt]->GetXaxis()->SetTitleSize(0.05);
                    fh1_Tcal1Hit_TofRaw[cpt]->GetYaxis()->SetLabelSize(0.05);
                    fh1_Tcal1Hit_TofRaw[cpt]->GetYaxis()->SetTitleSize(0.05);
                    cTcal_Tof->cd(cpt + 1);
                    fh1_Tcal1Hit_TofRaw[cpt]->Draw();
                    cpt++;
                }
            }
        }
    }
    // --- ------------ --- //
    // --- CAL LEVELS --- //
    // --- ------------ --- //

    // === get access to pos cal data ===//
    fPosCal = dynamic_cast<TClonesArray*>(mgr->GetObject("FrsSciPosCalData"));
    if (!fPosCal)
    {
        LOG(info) << "R3BOnlineSpectraFrsSci::Init() :: FrsSciPosCal not found";
        LOG(info) << "    No FrsSciPosCal spectra ---> ok";
    }
    else
    {
        // === declare and create TCanvas and Histogram === //

        sprintf(Name1, "Cal_PosRaw");
        cCal_PosRaw = new TCanvas(Name1, Name1, 10, 10, 800, 700);
        cCal_PosRaw->Divide(1, fNbDets);
        fh1_Cal_PosRaw = new TH1D*[fNbDets];

        sprintf(Name1, "Cal_PosCal");
        cCal_PosCal = new TCanvas(Name1, Name1, 10, 10, 800, 700);
        cCal_PosCal->Divide(1, fNbDets);
        fh1_Cal_PosCal = new TH1D*[fNbDets];

        for (UShort_t i = 0; i < fNbDets; i++)
        {
            // === TH1F: Raw Position in Ns at Cal level=== //
            sprintf(Name1, "FrsSci%i_PosRaw", i + 1);
            if (fh1_Tcal1Hit_PosRaw[i])
            {
                Int_t nbins = fh1_Tcal1Hit_PosRaw[i]->GetNbinsX();
                fh1_Cal_PosRaw[i] =
                    new TH1D(Name1,
                             Name1,
                             nbins,
                             fh1_Tcal1Hit_PosRaw[i]->GetBinLowEdge(1),
                             fh1_Tcal1Hit_PosRaw[i]->GetBinLowEdge(nbins) + fh1_Tcal1Hit_PosRaw[i]->GetBinWidth(nbins));
            }
            else
                fh1_Cal_PosRaw[i] = new TH1D(Name1,
                                             Name1,
                                             static_cast<int>(2. * (fpos_range_max - fpos_range_min)),
                                             fpos_range_min,
                                             fpos_range_max);
            fh1_Cal_PosRaw[i]->GetXaxis()->SetTitle("Raw Positon [ns] in red CAL level, in blue TCAL-MULT1");
            fh1_Cal_PosRaw[i]->GetYaxis()->SetTitle("number of counts with mult1");
            fh1_Cal_PosRaw[i]->GetXaxis()->CenterTitle(true);
            fh1_Cal_PosRaw[i]->GetYaxis()->CenterTitle(true);
            fh1_Cal_PosRaw[i]->GetXaxis()->SetLabelSize(0.05);
            fh1_Cal_PosRaw[i]->GetXaxis()->SetTitleSize(0.05);
            fh1_Cal_PosRaw[i]->GetYaxis()->SetLabelSize(0.05);
            fh1_Cal_PosRaw[i]->GetYaxis()->SetTitleSize(0.05);
            cCal_PosRaw->cd(i + 1);
            fh1_Cal_PosRaw[i]->SetLineColor(kRed);
            fh1_Cal_PosRaw[i]->Draw();
            fh1_Tcal1Hit_PosRaw[i]->Draw("same");

            // === TH1F: calibrated Position in Ns at Cal level=== //
            sprintf(Name1, "FrsSci%i_PosCal", i + 1);
            fh1_Cal_PosCal[i] = new TH1D(Name1, Name1, 30000, -150, 150);
            fh1_Cal_PosCal[i]->GetXaxis()->SetTitle("Calibrated Positon [mm]");
            fh1_Cal_PosCal[i]->GetYaxis()->SetTitle("number of counts");
            fh1_Cal_PosCal[i]->GetXaxis()->CenterTitle(true);
            fh1_Cal_PosCal[i]->GetYaxis()->CenterTitle(true);
            fh1_Cal_PosCal[i]->GetXaxis()->SetLabelSize(0.05);
            fh1_Cal_PosCal[i]->GetXaxis()->SetTitleSize(0.05);
            fh1_Cal_PosCal[i]->GetYaxis()->SetLabelSize(0.05);
            fh1_Cal_PosCal[i]->GetYaxis()->SetTitleSize(0.05);
            cCal_PosCal->cd(i + 1);
            fh1_Cal_PosCal[i]->Draw();
        }
    }
    // === get access to tof cal data ===//
    fTofCal = dynamic_cast<TClonesArray*>(mgr->GetObject("FrsSciTofCalData"));
    if (!fTofCal)
    {
        LOG(info) << "R3BOnlineSpectraFrsSci::Init() :: FrsSciTofCal not found";
        LOG(info) << "    No FrsSciTofCal spectra ---> ok";
    }

    else
    {
        sprintf(Name1, "Cal_TofRaw");
        cCal_TofRaw = new TCanvas(Name1, Name1, 10, 10, 800, 700);
        cCal_TofRaw->Divide(1, fNbTofs);
        fh1_Cal_TofRaw = new TH1D*[fNbTofs];

        sprintf(Name1, "Cal_TofCal");
        cCal_TofCal = new TCanvas(Name1, Name1, 10, 10, 800, 700);
        cCal_TofCal->Divide(1, fNbTofs);
        fh1_Cal_TofCal = new TH1D*[fNbTofs];

        sprintf(Name1, "Cal_Beta");
        cCal_Beta = new TCanvas(Name1, Name1, 10, 10, 800, 700);
        cCal_Beta->Divide(1, fNbTofs);
        fh1_Cal_Beta = new TH1D*[fNbTofs];

        sprintf(Name1, "Cal_BRho_AoQ");
        cCal_AoQ = new TCanvas(Name1, Name1, 10, 10, 800, 700);
        cCal_AoQ->Divide(2, fNbTofs);
        fh1_Cal_BRho = new TH1D*[fNbTofs];
        fh1_Cal_AoQ = new TH1D*[fNbTofs];

        sprintf(Name1, "Cal_BRhoVsPos");
        cCal_BRhoVsPos = new TCanvas(Name1, Name1, 10, 10, 800, 700);
        cCal_BRhoVsPos->Divide(2, fNbTofs);
        fh2_Cal_BRhoVsPosSta = new TH2D*[fNbTofs];
        fh2_Cal_BRhoVsPosSto = new TH2D*[fNbTofs];

        sprintf(Name1, "Cal_AoQVsPos");
        cCal_AoQVsPos = new TCanvas(Name1, Name1, 10, 10, 800, 700);
        cCal_AoQVsPos->Divide(2, fNbTofs);
        fh2_Cal_AoQVsPosSta = new TH2D*[fNbTofs];
        fh2_Cal_AoQVsPosSto = new TH2D*[fNbTofs];

        if (fNbTofs > 0)
        {
            UShort_t cpt = 0;
            for (UShort_t sta = 1; sta < fNbDets; sta++)
            {
                for (UShort_t sto = sta + 1; sto <= fNbDets; sto++)
                {
                    sprintf(Name1, "TofRaw_FrsSci%i_to_FrsSci%i", sta, sto);
                    if (fh1_Tcal1Hit_TofRaw[cpt])
                    {
                        Int_t nbins = fh1_Tcal1Hit_TofRaw[cpt]->GetNbinsX();
                        fh1_Cal_TofRaw[cpt] = new TH1D(Name1,
                                                       Name1,
                                                       nbins,
                                                       fh1_Tcal1Hit_TofRaw[cpt]->GetBinLowEdge(1),
                                                       fh1_Tcal1Hit_TofRaw[cpt]->GetBinLowEdge(nbins) +
                                                           fh1_Tcal1Hit_TofRaw[cpt]->GetBinWidth(nbins));
                    }
                    else
                        fh1_Cal_TofRaw[cpt] = new TH1D(Name1, Name1, 150000, 200, 2500);
                    fh1_Cal_TofRaw[cpt]->GetXaxis()->SetTitle("Raw Tof [ns] in red CAL level, in blue TCAL-MULT1");
                    fh1_Cal_TofRaw[cpt]->GetYaxis()->SetTitle("number of counts with mult1");
                    fh1_Cal_TofRaw[cpt]->GetXaxis()->CenterTitle(true);
                    fh1_Cal_TofRaw[cpt]->GetYaxis()->CenterTitle(true);
                    fh1_Cal_TofRaw[cpt]->GetXaxis()->SetLabelSize(0.05);
                    fh1_Cal_TofRaw[cpt]->GetXaxis()->SetTitleSize(0.05);
                    fh1_Cal_TofRaw[cpt]->GetYaxis()->SetLabelSize(0.05);
                    fh1_Cal_TofRaw[cpt]->GetYaxis()->SetTitleSize(0.05);
                    cCal_TofRaw->cd(cpt + 1);
                    fh1_Cal_TofRaw[cpt]->SetLineColor(kRed);
                    fh1_Cal_TofRaw[cpt]->Draw();
                    fh1_Tcal1Hit_TofRaw[cpt]->Draw("same");

                    sprintf(Name1, "TofCal_FrsSci%i_to_FrsSci%i", sta, sto);
                    fh1_Cal_TofCal[cpt] = new TH1D(Name1, Name1, 16000, 200, 1000.);
                    fh1_Cal_TofCal[cpt]->GetXaxis()->SetTitle("calibrated Tof [ns]");
                    fh1_Cal_TofCal[cpt]->GetYaxis()->SetTitle("number of counts");
                    fh1_Cal_TofCal[cpt]->GetXaxis()->CenterTitle(true);
                    fh1_Cal_TofCal[cpt]->GetYaxis()->CenterTitle(true);
                    fh1_Cal_TofCal[cpt]->GetXaxis()->SetLabelSize(0.05);
                    fh1_Cal_TofCal[cpt]->GetXaxis()->SetTitleSize(0.05);
                    fh1_Cal_TofCal[cpt]->GetYaxis()->SetLabelSize(0.05);
                    fh1_Cal_TofCal[cpt]->GetYaxis()->SetTitleSize(0.05);
                    cCal_TofCal->cd(cpt + 1);
                    fh1_Cal_TofCal[cpt]->Draw();

                    sprintf(Name1, "Beta_FrsSci%i_to_FrsSci%i", sta, sto);
                    fh1_Cal_Beta[cpt] = new TH1D(Name1, Name1, 2000, 0.7, 0.72);
                    fh1_Cal_Beta[cpt]->GetXaxis()->SetTitle("Beta");
                    fh1_Cal_Beta[cpt]->GetYaxis()->SetTitle("number of counts");
                    fh1_Cal_Beta[cpt]->GetXaxis()->CenterTitle(true);
                    fh1_Cal_Beta[cpt]->GetYaxis()->CenterTitle(true);
                    fh1_Cal_Beta[cpt]->GetXaxis()->SetLabelSize(0.05);
                    fh1_Cal_Beta[cpt]->GetXaxis()->SetTitleSize(0.05);
                    fh1_Cal_Beta[cpt]->GetYaxis()->SetLabelSize(0.05);
                    fh1_Cal_Beta[cpt]->GetYaxis()->SetTitleSize(0.05);
                    cCal_Beta->cd(cpt + 1);
                    fh1_Cal_Beta[cpt]->Draw();

                    sprintf(Name1, "BRho_FrsSci%i_to_FrsSci%i", sta, sto);
                    fh1_Cal_BRho[cpt] = new TH1D(Name1, Name1, 10000, 9.5, 9.7);
                    fh1_Cal_BRho[cpt]->GetXaxis()->SetTitle("BRho [T.m]");
                    fh1_Cal_BRho[cpt]->GetYaxis()->SetTitle("number of counts");
                    fh1_Cal_BRho[cpt]->GetXaxis()->CenterTitle(true);
                    fh1_Cal_BRho[cpt]->GetYaxis()->CenterTitle(true);
                    fh1_Cal_BRho[cpt]->GetXaxis()->SetLabelSize(0.05);
                    fh1_Cal_BRho[cpt]->GetXaxis()->SetTitleSize(0.05);
                    fh1_Cal_BRho[cpt]->GetYaxis()->SetLabelSize(0.05);
                    fh1_Cal_BRho[cpt]->GetYaxis()->SetTitleSize(0.05);
                    cCal_AoQ->cd(2 * cpt + 1);
                    fh1_Cal_BRho[cpt]->Draw();

                    sprintf(Name1, "AoQ_FrsSci%i_to_FrsSci%i", sta, sto);
                    fh1_Cal_AoQ[cpt] = new TH1D(Name1, Name1, 30000, 0.5, 3.5);
                    fh1_Cal_AoQ[cpt]->GetXaxis()->SetTitle("A/Q");
                    fh1_Cal_AoQ[cpt]->GetYaxis()->SetTitle("number of counts");
                    fh1_Cal_AoQ[cpt]->GetXaxis()->CenterTitle(true);
                    fh1_Cal_AoQ[cpt]->GetYaxis()->CenterTitle(true);
                    fh1_Cal_AoQ[cpt]->GetXaxis()->SetLabelSize(0.05);
                    fh1_Cal_AoQ[cpt]->GetXaxis()->SetTitleSize(0.05);
                    fh1_Cal_AoQ[cpt]->GetYaxis()->SetLabelSize(0.05);
                    fh1_Cal_AoQ[cpt]->GetYaxis()->SetTitleSize(0.05);
                    cCal_AoQ->cd(2 * (cpt + 1));
                    fh1_Cal_AoQ[cpt]->Draw();

                    sprintf(Name1, "BRho_Vs_PosFrsSci%i", sta);
                    fh2_Cal_BRhoVsPosSta[cpt] = new TH2D(Name1, Name1, 150, -2, 0.5, 500, 9.5, 9.7);
                    fh2_Cal_BRhoVsPosSta[cpt]->GetXaxis()->SetTitle("Pos [mm] if calibrated");
                    fh2_Cal_BRhoVsPosSta[cpt]->GetYaxis()->SetTitle("BRho [T.m]");
                    fh2_Cal_BRhoVsPosSta[cpt]->GetXaxis()->CenterTitle(true);
                    fh2_Cal_BRhoVsPosSta[cpt]->GetYaxis()->CenterTitle(true);
                    fh2_Cal_BRhoVsPosSta[cpt]->GetXaxis()->SetLabelSize(0.05);
                    fh2_Cal_BRhoVsPosSta[cpt]->GetXaxis()->SetTitleSize(0.05);
                    fh2_Cal_BRhoVsPosSta[cpt]->GetYaxis()->SetLabelSize(0.05);
                    fh2_Cal_BRhoVsPosSta[cpt]->GetYaxis()->SetTitleSize(0.05);
                    cCal_BRhoVsPos->cd(2 * cpt + 1);
                    fh2_Cal_BRhoVsPosSta[cpt]->Draw("col");

                    sprintf(Name1, "BRho_Vs_PosRawFrsSci%i", sto);
                    fh2_Cal_BRhoVsPosSto[cpt] = new TH2D(Name1, Name1, 150, 0.5, 2., 500, 9.5, 9.7);
                    fh2_Cal_BRhoVsPosSto[cpt]->GetXaxis()->SetTitle("Pos [mm] if calibrated");
                    fh2_Cal_BRhoVsPosSto[cpt]->GetYaxis()->SetTitle("BRho [T.m]");
                    fh2_Cal_BRhoVsPosSto[cpt]->GetXaxis()->CenterTitle(true);
                    fh2_Cal_BRhoVsPosSto[cpt]->GetYaxis()->CenterTitle(true);
                    fh2_Cal_BRhoVsPosSto[cpt]->GetXaxis()->SetLabelSize(0.05);
                    fh2_Cal_BRhoVsPosSto[cpt]->GetXaxis()->SetTitleSize(0.05);
                    fh2_Cal_BRhoVsPosSto[cpt]->GetYaxis()->SetLabelSize(0.05);
                    fh2_Cal_BRhoVsPosSto[cpt]->GetYaxis()->SetTitleSize(0.05);
                    cCal_BRhoVsPos->cd(2 * (cpt + 1));
                    fh2_Cal_BRhoVsPosSto[cpt]->Draw("col");

                    sprintf(Name1, "AoQ_Vs_PosFrsSci%i", sta);
                    fh2_Cal_AoQVsPosSta[cpt] = new TH2D(Name1, Name1, 150, -2, 0.5, 600, 2.99, 3.15);
                    fh2_Cal_AoQVsPosSta[cpt]->GetXaxis()->SetTitle("Pos [mm] if calibrated");
                    fh2_Cal_AoQVsPosSta[cpt]->GetYaxis()->SetTitle("AoQ [T.m]");
                    fh2_Cal_AoQVsPosSta[cpt]->GetXaxis()->CenterTitle(true);
                    fh2_Cal_AoQVsPosSta[cpt]->GetYaxis()->CenterTitle(true);
                    fh2_Cal_AoQVsPosSta[cpt]->GetXaxis()->SetLabelSize(0.05);
                    fh2_Cal_AoQVsPosSta[cpt]->GetXaxis()->SetTitleSize(0.05);
                    fh2_Cal_AoQVsPosSta[cpt]->GetYaxis()->SetLabelSize(0.05);
                    fh2_Cal_AoQVsPosSta[cpt]->GetYaxis()->SetTitleSize(0.05);
                    cCal_AoQVsPos->cd(2 * cpt + 1);
                    fh2_Cal_AoQVsPosSta[cpt]->Draw("col");

                    sprintf(Name1, "AoQ_Vs_PosRawFrsSci%i", sto);
                    fh2_Cal_AoQVsPosSto[cpt] = new TH2D(Name1, Name1, 150, 0.5, 2., 600, 2.99, 3.15);
                    fh2_Cal_AoQVsPosSto[cpt]->GetXaxis()->SetTitle("Pos [mm] if calibrated");
                    fh2_Cal_AoQVsPosSto[cpt]->GetYaxis()->SetTitle("AoQ [T.m]");
                    fh2_Cal_AoQVsPosSto[cpt]->GetXaxis()->CenterTitle(true);
                    fh2_Cal_AoQVsPosSto[cpt]->GetYaxis()->CenterTitle(true);
                    fh2_Cal_AoQVsPosSto[cpt]->GetXaxis()->SetLabelSize(0.05);
                    fh2_Cal_AoQVsPosSto[cpt]->GetXaxis()->SetTitleSize(0.05);
                    fh2_Cal_AoQVsPosSto[cpt]->GetYaxis()->SetLabelSize(0.05);
                    fh2_Cal_AoQVsPosSto[cpt]->GetYaxis()->SetTitleSize(0.05);
                    cCal_AoQVsPos->cd(2 * (cpt + 1));
                    fh2_Cal_AoQVsPosSto[cpt]->Draw("col");

                    cpt++;
                }
            }
        }
    }

    // --- --------------- --- //
    // --- MAIN FOLDER-Sci --- //
    // --- --------------- --- //
    TFolder* mainfolMap = new TFolder("FrsSciMap", "FrsSci Map info");
    mainfolMap->Add(cMap_FT);
    mainfolMap->Add(cMap_Mult1D);
    mainfolMap->Add(cMap_Mult2D);
    run->AddObject(mainfolMap);

    if (fTcal)
    {
        TFolder* mainfolTcal = new TFolder("FrsSciTcal", "FrsSci Tcal info");
        mainfolTcal->Add(cTcal_Pos);
        mainfolTcal->Add(cTcal_Tof);
        run->AddObject(mainfolTcal);
    }
    if (fPosCal || fTofCal)
    {
        TFolder* mainfolCal = new TFolder("FrsSciCal", "FrsSci Cal info");
        if (fPosCal)
        {
            mainfolCal->Add(cCal_PosRaw);
            mainfolCal->Add(cCal_PosCal);
        }
        if (fTofCal)
        {
            mainfolCal->Add(cCal_TofRaw);
            mainfolCal->Add(cCal_TofCal);
            mainfolCal->Add(cCal_Beta);
            mainfolCal->Add(cCal_AoQ);
            mainfolCal->Add(cCal_BRhoVsPos);
            mainfolCal->Add(cCal_AoQVsPos);
        }
        run->AddObject(mainfolCal);
    }

    // Register command to reset histograms
    run->GetHttpServer()->RegisterCommand("Reset_FRSSCI_HIST", Form("/Objects/%s/->Reset_Histo()", GetName()));

    return kSUCCESS;
}

void R3BOnlineSpectraFrsSci::Reset_Histo()
{
    LOG(info) << "R3BOnlineSpectraFrsSci::Reset_Histo";
    for (UShort_t i = 0; i < fNbDets; i++)
    {
        fh2_Map_mult[i]->Reset();
        fh2_Map_multRvsL[i]->Reset();

        // === FINE TIME === //
        for (UShort_t j = 0; j < fNbPmts; j++)
        {
            fh1_Map_finetime[i * fNbPmts + j]->Reset();
            fh1_Map_mult[i * fNbPmts + j]->Reset();
        }

        if (fTcal)
            fh1_Tcal1Hit_PosRaw[i]->Reset();
        if (fPosCal)
        {
            fh1_Cal_PosRaw[i]->Reset();
            fh1_Cal_PosCal[i]->Reset();
        }
    }

    for (UShort_t i = 0; i < fNbTofs; i++)
    {
        if (fTcal)
            fh1_Tcal1Hit_TofRaw[i]->Reset();
        if (fTofCal)
        {
            fh1_Cal_TofRaw[i]->Reset();
            fh1_Cal_TofCal[i]->Reset();
            fh1_Cal_Beta[i]->Reset();
            fh1_Cal_BRho[i]->Reset();
            fh1_Cal_AoQ[i]->Reset();
            fh2_Cal_BRhoVsPosSta[i]->Reset();
            fh2_Cal_BRhoVsPosSto[i]->Reset();
            fh2_Cal_AoQVsPosSta[i]->Reset();
            fh2_Cal_AoQVsPosSto[i]->Reset();
        }
    }
}

void R3BOnlineSpectraFrsSci::Exec(Option_t* option)
{

    FairRootManager* mgr = FairRootManager::Instance();
    if (NULL == mgr)
    {
        LOG(fatal) << "R3BOnlineSpectraFrsSci::Exec FairRootManager not found";
    }
    UInt_t nHits;
    UShort_t iDet;  // 0-based
    UShort_t iPmt;  // 0-based
    UShort_t iRank; // 0-based
    UShort_t cpt;

    UInt_t multMap[fNbDets * fNbPmts];
    UInt_t multTcal[fNbDets * fNbPmts];
    Double_t Traw[fNbDets * fNbPmts];

    for (UShort_t i = 0; i < fNbDets; i++)
    {
        for (UShort_t j = 0; j < fNbPmts; j++)
        {
            multMap[i * fNbPmts + j] = 0;
            multTcal[i * fNbPmts + j] = 0;
        }
    }

    if (fMapped)
    {
        if (fMapped->GetEntriesFast() > 0)
        {
            // --- ---------------- --- //
            // --- read mapped data --- //
            // --- ---------------- --- //
            nHits = fMapped->GetEntriesFast();
            for (UInt_t ihit = 0; ihit < nHits; ihit++)
            {
                R3BFrsSciMappedData* hitmapped = dynamic_cast<R3BFrsSciMappedData*>(fMapped->At(ihit));
                if (!hitmapped)
                    continue;
                iDet = hitmapped->GetDetector() - 1;
                iPmt = hitmapped->GetPmt() - 1;
                multMap[iDet * fNbPmts + iPmt]++;
                fh1_Map_finetime[iDet * fNbPmts + iPmt]->Fill(hitmapped->GetTimeFine());
            } // end of loop over mapped data

            // --- ----------------------------------------- --- //
            // --- filling some histogramms outside the loop --- //
            // --- ----------------------------------------- --- //
            for (UShort_t i = 0; i < fNbDets; i++)
            {
                if (fNbPmts > 1)
                    fh2_Map_multRvsL[i]->Fill(multMap[i * fNbPmts + 1], multMap[i * fNbPmts]);

                for (UShort_t j = 0; j < fNbPmts; j++)
                {
                    fh1_Map_mult[i * fNbPmts + j]->Fill(multMap[i * fNbPmts + j]);
                    fh2_Map_mult[i]->Fill(j + 1, multMap[i * fNbPmts + j]);
                }
            }

        } // --- end of if (fMapped->GetEntries() > 0) --- //

        if (fTcal)
        {
            if (fTcal->GetEntries() > 0)
            {
                // --- -------------- --- //
                // --- read tcal data --- //
                // --- -------------- --- //
                nHits = fTcal->GetEntriesFast();
                for (UInt_t ihit = 0; ihit < nHits; ihit++)
                {
                    R3BFrsSciTcalData* hittcal = dynamic_cast<R3BFrsSciTcalData*>(fTcal->At(ihit));
                    if (!hittcal)
                        continue;
                    iDet = hittcal->GetDetector() - 1;
                    iPmt = hittcal->GetPmt() - 1;
                    multTcal[iDet * fNbPmts + iPmt]++;
                    Traw[iDet * fNbPmts + iPmt] = hittcal->GetRawTimeNs();
                } // end of loop over tcal data

                // --- ------------------- --- //
                // --- filling histogramms --- //
                // --- ------------------- --- //
                if (fNbPmts > 1 && fNbDets > 0)
                {
                    // PosRaw [ns] = Tright - Tleft -> x increasing from right to left
                    for (UShort_t i = 0; i < fNbDets; i++)
                    {
                        if (multTcal[i * fNbPmts] == 1 && multTcal[i * fNbPmts + 1] == 1)
                        {
                            fh1_Tcal1Hit_PosRaw[i]->Fill((Float_t)(Traw[i * fNbPmts] - Traw[i * fNbPmts + 1]));
                        }
                    }
                    // TofRaw [ns] relatively to Tref
                    if (fNbPmts > 2 && fNbDets > 1 && fNbTofs > 0)
                    {
                        cpt = 0;
                        for (UShort_t sta = 0; sta < fNbDets - 1; sta++)
                        {
                            for (UShort_t sto = sta + 1; sto < fNbDets; sto++)
                            {
                                if (multTcal[sta * fNbPmts] == 1 && multTcal[sta * fNbPmts + 1] == 1 &&
                                    multTcal[sta * fNbPmts + 2] == 1 && multTcal[sto * fNbPmts] == 1 &&
                                    multTcal[sto * fNbPmts + 1] == 1 && multTcal[sto * fNbPmts + 2] == 1)
                                {
                                    fh1_Tcal1Hit_TofRaw[cpt]->Fill(
                                        //    0.5 * (Traw[sto * fNbPmts] + Traw[sto * fNbPmts + 1]) -
                                        //    0.5 * (Traw[sta * fNbPmts] + Traw[sta * fNbPmts + 1])) ;
                                        0.5 * (Traw[sto * fNbPmts] + Traw[sto * fNbPmts + 1]) -
                                        Traw[sto * fNbPmts + 2] -
                                        0.5 * (Traw[sta * fNbPmts] + Traw[sta * fNbPmts + 1]) +
                                        Traw[sta * fNbPmts + 2]);
                                    cpt++;
                                } // end of if mult = 1
                            }     // end of loop over stop FrsSci
                        }         // end of loop over start FrsSco
                    }
                }
            } // --- end of if Tcal->GetEntries() >0 --- //

            // --- ---------------------------- --- //
            // --- read cal data and fill histo --- //
            // --- ---------------------------- --- //
            if (fPosCal)
            {
                if (fPosCal->GetEntriesFast() > 0)
                {
                    // --- ---------------- --- //
                    // --- read poscal data --- //
                    // --- ---------------- --- //
                    nHits = fPosCal->GetEntriesFast();
                    for (UInt_t ihit = 0; ihit < nHits; ihit++)
                    {
                        R3BFrsSciPosCalData* hitposcal = dynamic_cast<R3BFrsSciPosCalData*>(fPosCal->At(ihit));
                        if (!hitposcal)
                            continue;
                        iDet = hitposcal->GetDetector() - 1;
                        fh1_Cal_PosRaw[iDet]->Fill(hitposcal->GetRawPosNs());
                        fh1_Cal_PosCal[iDet]->Fill(hitposcal->GetCalPosMm());
                    } // end of loop over pos cal data
                }     // --- end of if fPosCal->GetEntries() >0 --- //
            }         // end of if fPosCal
            if (fTofCal)
            {
                if (fTofCal->GetEntriesFast() > 0)
                {
                    // --- ---------------- --- //
                    // --- read tofcal data --- //
                    // --- ---------------- --- //
                    nHits = fTofCal->GetEntriesFast();
                    for (UInt_t ihit = 0; ihit < nHits; ihit++)
                    {
                        R3BFrsSciTofCalData* hittofcal = dynamic_cast<R3BFrsSciTofCalData*>(fTofCal->At(ihit));
                        if (!hittofcal)
                            continue;
                        iRank = hittofcal->GetRank() - 1;
                        fh1_Cal_TofRaw[iRank]->Fill(hittofcal->GetRawTofNs());
                        fh1_Cal_TofCal[iRank]->Fill(hittofcal->GetCalTofNs());
                        fh1_Cal_Beta[iRank]->Fill(hittofcal->GetBeta());
                        fh1_Cal_BRho[iRank]->Fill(hittofcal->GetBRho());
                        fh1_Cal_AoQ[iRank]->Fill(hittofcal->GetAoQ());
                        fh2_Cal_BRhoVsPosSta[iRank]->Fill(hittofcal->GetCalPosStaMm(), hittofcal->GetBRho());
                        fh2_Cal_BRhoVsPosSto[iRank]->Fill(hittofcal->GetCalPosStoMm(), hittofcal->GetBRho());
                        fh2_Cal_AoQVsPosSta[iRank]->Fill(hittofcal->GetCalPosStaMm(), hittofcal->GetAoQ());
                        fh2_Cal_AoQVsPosSto[iRank]->Fill(hittofcal->GetCalPosStoMm(), hittofcal->GetAoQ());
                    } // end of loop over tof cal data
                }     // --- end of if fTofCal->GetEntries() >0 --- //
            }         // end of if fTofCal

        } // end of if Tcal

        fNEvents++;

    } // end of if fMapped
}

// -----   Public method Finish   -----------------------------------------------
void R3BOnlineSpectraFrsSci::FinishEvent()
{
    if (fMapped)
    {
        fMapped->Clear();
    }
    if (fTcal)
    {
        fTcal->Clear();
    }
    if (fPosCal)
    {
        fPosCal->Clear();
    }
    if (fTofCal)
    {
        fTofCal->Clear();
    }
}

void R3BOnlineSpectraFrsSci::FinishTask()
{
    if (fMapped)
    {
        for (UShort_t i = 0; i < fNbDets; i++)
        {
            fh2_Map_mult[i]->Write();
            fh2_Map_multRvsL[i]->Write();
            for (UShort_t j = 0; j < fNbPmts; j++)
            {
                fh1_Map_finetime[i * fNbPmts + j]->Write();
                fh1_Map_mult[i * fNbPmts + j]->Write();
            }
        }
        if (fTcal)
        {
            for (UShort_t i = 0; i < fNbDets; i++)
            {
                fh1_Tcal1Hit_PosRaw[i]->Write();
            }
            for (UShort_t i = 0; i < fNbTofs; i++)
            {
                fh1_Tcal1Hit_TofRaw[i]->Write();
            }
        }
        if (fPosCal)
        {
            for (UShort_t i = 0; i < fNbDets; i++)
            {
                fh1_Cal_PosRaw[i]->Write();
                fh1_Cal_PosRaw[i]->Write();
            }
        }
        if (fTofCal)
        {
            for (UShort_t i = 0; i < fNbTofs; i++)
            {
                fh1_Cal_TofRaw[i]->Write();
                fh1_Cal_TofCal[i]->Write();
                fh1_Cal_Beta[i]->Write();
                fh1_Cal_BRho[i]->Write();
                fh1_Cal_AoQ[i]->Write();
                fh2_Cal_BRhoVsPosSta[i]->Write();
                fh2_Cal_BRhoVsPosSto[i]->Write();
                fh2_Cal_AoQVsPosSta[i]->Write();
                fh2_Cal_AoQVsPosSto[i]->Write();
            }
        }
    }
}
ClassImp(R3BOnlineSpectraFrsSci);
