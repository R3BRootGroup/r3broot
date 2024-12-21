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
#include "R3BOnlineSpectraSci2.h"
#include "R3BSci2MappedData.h"
#include "R3BSci2TcalData.h"

R3BOnlineSpectraSci2::R3BOnlineSpectraSci2()
    : R3BOnlineSpectraSci2("Sci2OnlineSpectra", 1)
{
}

R3BOnlineSpectraSci2::R3BOnlineSpectraSci2(const char* name, Int_t iVerbose)
    : FairTask(name, iVerbose)
    , fMapped(NULL)
    , fTcal(NULL)
    , fNEvents(0)
    , fNbDetectors(1)
    , fNbChannels(3)
    , fTpat1(-1)
    , fTpat2(-1)
{
}

R3BOnlineSpectraSci2::~R3BOnlineSpectraSci2()
{
    LOG(debug) << "R3BOnlineSpectraSci2::Destructor";
    if (fMapped)
        delete fMapped;
    if (fTcal)
        delete fTcal;
}

InitStatus R3BOnlineSpectraSci2::Init()
{
    LOG(info) << "R3BOnlineSpectraSci2::Init()";

    // try to get a handle on the EventHeader. EventHeader may not be
    // present though and hence may be null. Take care when using.

    FairRootManager* mgr = FairRootManager::Instance();
    if (NULL == mgr)
        LOG(fatal) << "R3BOnlineSpectraSci2::Init FairRootManager not found";

    fEventHeader = dynamic_cast<R3BEventHeader*>(mgr->GetObject("EventHeader."));
    if (!fEventHeader)
        fEventHeader = dynamic_cast<R3BEventHeader*>(mgr->GetObject("R3BEventHeader"));

    FairRunOnline* run = FairRunOnline::Instance();
    run->GetHttpServer()->Register("", this);

    // --- ------------------------- --- //
    // --- get access to mapped data --- //
    // --- ------------------------- --- //
    fMapped = dynamic_cast<TClonesArray*>(mgr->GetObject("Sci2Mapped"));
    if (!fMapped)
    {
        LOG(fatal) << "Sci2Mapped not found";
        return (kFATAL);
    }

    // --- ------------------------ --- //
    // --- get access to tcal data  --- //
    // --- ------------------------ --- //
    fTcal = dynamic_cast<TClonesArray*>(mgr->GetObject("Sci2Tcal"));
    if (!fTcal)
    {
        LOG(warn) << "Sci2Tcal not found: is OK";
    }

    // --- ------------------------------ --- //
    // --- declare TCanvas and Histograms --- //
    // --- ------------------------------ --- //

    // === 1D-MULTIPLICITY without condition on TPAT === //
    fh1_multMap = new TH1I*[fNbDetectors * fNbChannels];
    fh1_multTcal = new TH1I*[fNbDetectors * fNbChannels];

    // === 1D-MULTIPLICITY with condition on TPAT === //
    fh1_multMap_condTpat = new TH1I*[fNbDetectors * fNbChannels];
    fh1_multTcal_condTpat = new TH1I*[fNbDetectors * fNbChannels];

    // === MAPPED LEVEL: FINE TIME AND MULT PER EVENT === //
    cMapped = new TCanvas*[fNbDetectors];
    fh1_finetime = new TH1I*[fNbDetectors * fNbChannels];
    fh2_mult = new TH2I*[fNbDetectors];
    fh2_mult_RvsL = new TH2I*[fNbDetectors];
    fh2_mult_RvsL_condTpat = new TH2I*[fNbDetectors];

    // === POSITION AT TCAL, SINGLE TCAL AND CAL LEVELS === //
    cPos = new TCanvas*[fNbDetectors];
    fh1_RawPos_TcalMult1 = new TH1F*[fNbDetectors];

    // === 2D multiplicity at Mapped level: Tref vs PMT === //
    fh2_mult_TrefVsPmt = new TH2I*[fNbDetectors * (fNbChannels - 1)];
    fh2_mult_TrefVsPmt_condTpat = new TH2I*[fNbDetectors * (fNbChannels - 1)];

    // --- ----------------------------- --- //
    // --- create TCanvas and Histograms --- //
    // --- ------------------------------ --- //
    char Name1[255];
    char Name2[255];

    // === TCanvas: 1D-MultMap === //
    sprintf(Name1, "MultMap");
    cMultMap = new TCanvas(Name1, Name1, 10, 10, 800, 700);
    cMultMap->Divide(fNbChannels, fNbDetectors);

    // === TCanvas: 1D-MultTcal === //
    sprintf(Name1, "MultTcal");
    cMultTcal = new TCanvas(Name1, Name1, 10, 10, 800, 700);
    cMultTcal->Divide(fNbChannels, fNbDetectors);

    // === TCanvas: 2D-MultMap === //
    sprintf(Name1, "MultMap2D");
    cMultMap2D = new TCanvas(Name1, Name1, 10, 10, 800, 700);
    cMultMap2D->Divide(fNbDetectors, 4);

    // === TCanvas: 2D-MultMap === //
    sprintf(Name1, "MultMap2D_RvsL");
    cMultMap2D_RvsL = new TCanvas(Name1, Name1, 10, 10, 800, 700);
    cMultMap2D_RvsL->Divide(2, fNbDetectors);

    for (Int_t i = 0; i < fNbDetectors; i++)
    {
        // === TCanvas: mapped level === //
        sprintf(Name1, "Sci2%i_MultAndFt", i + 1);
        cMapped[i] = new TCanvas(Name1, Name1, 10, 10, 800, 700);
        cMapped[i]->Divide(2, 2);

        // === TCanvas: position === //
        sprintf(Name1, "Sci2%i_Pos", i + 1);
        cPos[i] = new TCanvas(Name1, Name1, 10, 10, 800, 700);

        for (Int_t j = 0; j < fNbChannels; j++)
        {
            // === TH1I: 1D-mult at map level === //
            sprintf(Name1, "Sci2%i_Pmt%i_MultPerEvent_Mapped", i + 1, j + 1);
            sprintf(Name2,
                    "Sci2%i_Pmt%i_MultPerEvent_Mapped (blue no condition on TPAT, red condition on TPAT = 1 or 2)",
                    i + 1,
                    j + 1);
            fh1_multMap[i * fNbChannels + j] = new TH1I(Name1, Name1, 70, -0.5, 69.5);
            fh1_multMap[i * fNbChannels + j]->GetXaxis()->SetTitle("Multiplicity per event");
            fh1_multMap[i * fNbChannels + j]->GetYaxis()->SetTitle("Counts");
            fh1_multMap[i * fNbChannels + j]->GetXaxis()->CenterTitle(true);
            fh1_multMap[i * fNbChannels + j]->GetYaxis()->CenterTitle(true);
            fh1_multMap[i * fNbChannels + j]->GetXaxis()->SetLabelSize(0.05);
            fh1_multMap[i * fNbChannels + j]->GetXaxis()->SetTitleSize(0.05);
            fh1_multMap[i * fNbChannels + j]->GetYaxis()->SetLabelSize(0.05);
            fh1_multMap[i * fNbChannels + j]->GetYaxis()->SetTitleSize(0.05);
            fh1_multMap[i * fNbChannels + j]->SetLineWidth(2);
            fh1_multMap[i * fNbChannels + j]->SetLineColor(kBlue);
            cMultMap->cd(i * fNbChannels + j + 1);
            gPad->SetLogy();
            fh1_multMap[i * fNbChannels + j]->Draw("");

            sprintf(Name1, "Sci2%i_Pmt%i_MultPerEvent_Mapped_condTpat", i + 1, j + 1);
            fh1_multMap_condTpat[i * fNbChannels + j] = new TH1I(Name1, Name1, 70, -0.5, 69.5);
            fh1_multMap_condTpat[i * fNbChannels + j]->GetXaxis()->SetTitle("Multiplicity per event");
            fh1_multMap_condTpat[i * fNbChannels + j]->GetYaxis()->SetTitle("Counts");
            fh1_multMap_condTpat[i * fNbChannels + j]->GetXaxis()->CenterTitle(true);
            fh1_multMap_condTpat[i * fNbChannels + j]->GetYaxis()->CenterTitle(true);
            fh1_multMap_condTpat[i * fNbChannels + j]->GetXaxis()->SetLabelSize(0.05);
            fh1_multMap_condTpat[i * fNbChannels + j]->GetXaxis()->SetTitleSize(0.05);
            fh1_multMap_condTpat[i * fNbChannels + j]->GetYaxis()->SetLabelSize(0.05);
            fh1_multMap_condTpat[i * fNbChannels + j]->GetYaxis()->SetTitleSize(0.05);
            fh1_multMap_condTpat[i * fNbChannels + j]->SetLineWidth(1);
            fh1_multMap_condTpat[i * fNbChannels + j]->SetLineColor(kRed);
            cMultMap->cd(i * fNbChannels + j + 1);
            gPad->SetLogy();
            fh1_multMap_condTpat[i * fNbChannels + j]->Draw("sames");

            // === TH1I: 1D-mult at tcal level === //
            sprintf(Name1, "Sci2%i_Pmt%i_MultPerEvent_Tcal", i + 1, j + 1);
            sprintf(Name2,
                    "Sci2%i_Pmt%i_MultPerEvent_Tcal (blue no condition on Tpat, red condition on TPAT = 1 or 2)",
                    i + 1,
                    j + 1);
            fh1_multTcal[i * fNbChannels + j] = new TH1I(Name1, Name1, 70, -0.5, 69.5);
            fh1_multTcal[i * fNbChannels + j]->GetXaxis()->SetTitle("Multiplicity per event");
            fh1_multTcal[i * fNbChannels + j]->GetYaxis()->SetTitle("Counts");
            fh1_multTcal[i * fNbChannels + j]->GetXaxis()->CenterTitle(true);
            fh1_multTcal[i * fNbChannels + j]->GetYaxis()->CenterTitle(true);
            fh1_multTcal[i * fNbChannels + j]->GetXaxis()->SetLabelSize(0.05);
            fh1_multTcal[i * fNbChannels + j]->GetXaxis()->SetTitleSize(0.05);
            fh1_multTcal[i * fNbChannels + j]->GetYaxis()->SetLabelSize(0.05);
            fh1_multTcal[i * fNbChannels + j]->GetYaxis()->SetTitleSize(0.05);
            fh1_multTcal[i * fNbChannels + j]->SetLineWidth(2);
            fh1_multTcal[i * fNbChannels + j]->SetLineColor(kBlue);
            cMultTcal->cd(i * fNbChannels + j + 1);
            gPad->SetLogy();
            fh1_multTcal[i * fNbChannels + j]->Draw("");

            sprintf(Name1, "Sci2%i_Pmt%i_MultPerEvent_Tcal_condTpat", i + 1, j + 1);
            fh1_multTcal_condTpat[i * fNbChannels + j] = new TH1I(Name1, Name1, 70, -0.5, 69.5);
            fh1_multTcal_condTpat[i * fNbChannels + j]->GetXaxis()->SetTitle("Multiplicity per event");
            fh1_multTcal_condTpat[i * fNbChannels + j]->GetYaxis()->SetTitle("Counts");
            fh1_multTcal_condTpat[i * fNbChannels + j]->GetXaxis()->CenterTitle(true);
            fh1_multTcal_condTpat[i * fNbChannels + j]->GetYaxis()->CenterTitle(true);
            fh1_multTcal_condTpat[i * fNbChannels + j]->GetXaxis()->SetLabelSize(0.05);
            fh1_multTcal_condTpat[i * fNbChannels + j]->GetXaxis()->SetTitleSize(0.05);
            fh1_multTcal_condTpat[i * fNbChannels + j]->GetYaxis()->SetLabelSize(0.05);
            fh1_multTcal_condTpat[i * fNbChannels + j]->GetYaxis()->SetTitleSize(0.05);
            fh1_multTcal_condTpat[i * fNbChannels + j]->SetLineWidth(1);
            fh1_multTcal_condTpat[i * fNbChannels + j]->SetLineColor(kRed);
            cMultTcal->cd(i * fNbChannels + j + 1);
            gPad->SetLogy();
            fh1_multTcal_condTpat[i * fNbChannels + j]->Draw("sames");

            // === TH1F: fine time === //
            sprintf(Name1, "Sci2%i_FineTimeCh%i", i + 1, j + 1);
            fh1_finetime[i * fNbChannels + j] = new TH1I(Name1, Name1, 1000, 0, 1000);
            fh1_finetime[i * fNbChannels + j]->GetXaxis()->SetTitle("Fine time");
            fh1_finetime[i * fNbChannels + j]->GetYaxis()->SetTitle("Counts per bin");
            fh1_finetime[i * fNbChannels + j]->GetXaxis()->CenterTitle(true);
            fh1_finetime[i * fNbChannels + j]->GetYaxis()->CenterTitle(true);
            fh1_finetime[i * fNbChannels + j]->GetXaxis()->SetLabelSize(0.05);
            fh1_finetime[i * fNbChannels + j]->GetXaxis()->SetTitleSize(0.05);
            fh1_finetime[i * fNbChannels + j]->GetYaxis()->SetLabelSize(0.05);
            fh1_finetime[i * fNbChannels + j]->GetYaxis()->SetTitleSize(0.05);
            cMapped[i]->cd(j + 1);
            fh1_finetime[i * fNbChannels + j]->Draw("");
        }

        for (Int_t pmt = 0; pmt < fNbChannels - 1; pmt++)
        {
            sprintf(Name1, "Sci2%i_MultMap_TrefVsPmt%i", i + 1, pmt + 1);
            fh2_mult_TrefVsPmt[i * (fNbChannels - 1) + pmt] = new TH2I(Name1, Name1, 70, -0.5, 69.5, 5, -0.5, 4.5);
            fh2_mult_TrefVsPmt[i * (fNbChannels - 1) + pmt]->GetXaxis()->SetTitle("Mult Pmt");
            fh2_mult_TrefVsPmt[i * (fNbChannels - 1) + pmt]->GetYaxis()->SetTitle("Mult Tref");
            fh2_mult_TrefVsPmt[i * (fNbChannels - 1) + pmt]->GetXaxis()->CenterTitle(true);
            fh2_mult_TrefVsPmt[i * (fNbChannels - 1) + pmt]->GetYaxis()->CenterTitle(true);
            fh2_mult_TrefVsPmt[i * (fNbChannels - 1) + pmt]->GetXaxis()->SetLabelSize(0.05);
            fh2_mult_TrefVsPmt[i * (fNbChannels - 1) + pmt]->GetXaxis()->SetTitleSize(0.05);
            fh2_mult_TrefVsPmt[i * (fNbChannels - 1) + pmt]->GetYaxis()->SetLabelSize(0.05);
            fh2_mult_TrefVsPmt[i * (fNbChannels - 1) + pmt]->GetYaxis()->SetTitleSize(0.05);
            cMultMap2D->cd(i * (fNbChannels - 1) + pmt + 1);
            fh2_mult_TrefVsPmt[i * (fNbChannels - 1) + pmt]->Draw("col");

            sprintf(Name1, "Sci2%i_MultMap_TrefVsPmt%ii_condTpat", i + 1, pmt + 1);
            sprintf(Name2, "Sci2%i_MultMap_TrefVsPmt%i for Tpat = 1 or 2", i + 1, pmt + 1);
            fh2_mult_TrefVsPmt_condTpat[i * (fNbChannels - 1) + pmt] =
                new TH2I(Name1, Name2, 70, -0.5, 69.5, 5, -0.5, 4.5);
            fh2_mult_TrefVsPmt_condTpat[i * (fNbChannels - 1) + pmt]->GetXaxis()->SetTitle("Mult Pmt");
            fh2_mult_TrefVsPmt_condTpat[i * (fNbChannels - 1) + pmt]->GetYaxis()->SetTitle("Mult Tref");
            fh2_mult_TrefVsPmt_condTpat[i * (fNbChannels - 1) + pmt]->GetXaxis()->CenterTitle(true);
            fh2_mult_TrefVsPmt_condTpat[i * (fNbChannels - 1) + pmt]->GetYaxis()->CenterTitle(true);
            fh2_mult_TrefVsPmt_condTpat[i * (fNbChannels - 1) + pmt]->GetXaxis()->SetLabelSize(0.05);
            fh2_mult_TrefVsPmt_condTpat[i * (fNbChannels - 1) + pmt]->GetXaxis()->SetTitleSize(0.05);
            fh2_mult_TrefVsPmt_condTpat[i * (fNbChannels - 1) + pmt]->GetYaxis()->SetLabelSize(0.05);
            fh2_mult_TrefVsPmt_condTpat[i * (fNbChannels - 1) + pmt]->GetYaxis()->SetTitleSize(0.05);
            cMultMap2D->cd(i * (fNbChannels - 1) + pmt + 3);
            fh2_mult_TrefVsPmt_condTpat[i * (fNbChannels - 1) + pmt]->Draw("col");
        }

        // === TH1F: multiplicity per event and channel at mapped level === //
        sprintf(Name1, "Sci2%i_MultPerEvent", i + 1);
        fh2_mult[i] = new TH2I(Name1, Name1, fNbChannels, 0.5, fNbChannels + 0.5, 20, -0.5, 19.5);
        fh2_mult[i]->GetXaxis()->SetTitle("channel: 1=PMT R,    2=PMT L,    3=COMMON REF");
        fh2_mult[i]->GetYaxis()->SetTitle("multiplicity per event");
        fh2_mult[i]->GetXaxis()->CenterTitle(true);
        fh2_mult[i]->GetYaxis()->CenterTitle(true);
        fh2_mult[i]->GetXaxis()->SetLabelSize(0.05);
        fh2_mult[i]->GetXaxis()->SetTitleSize(0.05);
        fh2_mult[i]->GetYaxis()->SetLabelSize(0.05);
        fh2_mult[i]->GetYaxis()->SetTitleSize(0.05);
        cMapped[i]->cd(4);
        fh2_mult[i]->Draw("COL");

        // === TH1F: multiplicity per event and channel at mapped level === //
        sprintf(Name1, "Sci2%i_MultPerEvent_RvsL", i + 1);
        fh2_mult_RvsL[i] = new TH2I(Name1, Name1, 40, -1.5, 38.5, 40, -1.5, 38.5);
        fh2_mult_RvsL[i]->GetXaxis()->SetTitle("Multiplicity per event on the Left Pmt");
        fh2_mult_RvsL[i]->GetYaxis()->SetTitle("Multiplicity per event on the Right Pmt");
        fh2_mult_RvsL[i]->GetXaxis()->CenterTitle(true);
        fh2_mult_RvsL[i]->GetYaxis()->CenterTitle(true);
        fh2_mult_RvsL[i]->GetXaxis()->SetLabelSize(0.05);
        fh2_mult_RvsL[i]->GetXaxis()->SetTitleSize(0.05);
        fh2_mult_RvsL[i]->GetYaxis()->SetLabelSize(0.05);
        fh2_mult_RvsL[i]->GetYaxis()->SetTitleSize(0.05);
        cMultMap2D_RvsL->cd(i + 1);
        fh2_mult_RvsL[i]->Draw("COL");

        // === TH1F: multiplicity per event and channel at mapped level === //
        sprintf(Name1, "Sci2%i_MultPerEvent_RvsL_condTpat", i + 1);
        fh2_mult_RvsL_condTpat[i] = new TH2I(Name1, Name1, 40, -1.5, 38.5, 40, -1.5, 38.5);
        fh2_mult_RvsL_condTpat[i]->GetXaxis()->SetTitle("Multiplicity per event on the Left Pmt");
        fh2_mult_RvsL_condTpat[i]->GetYaxis()->SetTitle("Multiplicity per event on the Right Pmt");
        fh2_mult_RvsL_condTpat[i]->GetXaxis()->CenterTitle(true);
        fh2_mult_RvsL_condTpat[i]->GetYaxis()->CenterTitle(true);
        fh2_mult_RvsL_condTpat[i]->GetXaxis()->SetLabelSize(0.05);
        fh2_mult_RvsL_condTpat[i]->GetXaxis()->SetTitleSize(0.05);
        fh2_mult_RvsL_condTpat[i]->GetYaxis()->SetLabelSize(0.05);
        fh2_mult_RvsL_condTpat[i]->GetYaxis()->SetTitleSize(0.05);
        cMultMap2D_RvsL->cd(fNbDetectors + i + 1);
        fh2_mult_RvsL_condTpat[i]->Draw("COL");

        // === TH1F: raw position at tcal level if mult=1 === //
        sprintf(Name1, "Sci2%i_RawPos_Tcal_Mult1", i + 1);
        fh1_RawPos_TcalMult1[i] = new TH1F(Name1, Name1, 200000, -100, 100);
        fh1_RawPos_TcalMult1[i]->GetXaxis()->SetTitle(
            "(RIGHT,Wix. side) -->  raw position [ns, 1ps/bin] --> (LEFT,Mes. side) -->");
        fh1_RawPos_TcalMult1[i]->GetYaxis()->SetTitle("Counts per bin");
        fh1_RawPos_TcalMult1[i]->GetXaxis()->CenterTitle(true);
        fh1_RawPos_TcalMult1[i]->GetYaxis()->CenterTitle(true);
        fh1_RawPos_TcalMult1[i]->GetXaxis()->SetLabelSize(0.05);
        fh1_RawPos_TcalMult1[i]->GetXaxis()->SetTitleSize(0.05);
        fh1_RawPos_TcalMult1[i]->GetYaxis()->SetLabelSize(0.05);
        fh1_RawPos_TcalMult1[i]->GetYaxis()->SetTitleSize(0.05);
        cPos[i]->cd();
        fh1_RawPos_TcalMult1[i]->Draw("");

    } // end of loop over fNbDetectors

    // --- --------------- --- //
    // --- MAIN FOLDER-Sci --- //
    // --- --------------- --- //
    TFolder* mainfolSciMult = new TFolder("Sci2Mult", "Sci2 Mult info");
    mainfolSciMult->Add(cMultMap);
    mainfolSciMult->Add(cMultTcal);
    mainfolSciMult->Add(cMultMap2D);
    mainfolSciMult->Add(cMultMap2D_RvsL);

    TFolder* mainfolSci = new TFolder("Sci2_SA", "Sci2 info");
    for (Int_t i = 0; i < fNbDetectors; i++)
    {
        mainfolSci->Add(cMapped[i]);
        mainfolSci->Add(cPos[i]);
    }
    run->AddObject(mainfolSciMult);
    run->AddObject(mainfolSci);

    // Register command to reset histograms
    run->GetHttpServer()->RegisterCommand("Reset_SCI2_HIST", Form("/Objects/%s/->Reset_Histo()", GetName()));

    return kSUCCESS;
}

void R3BOnlineSpectraSci2::Reset_Histo()
{
    LOG(info) << "R3BOnlineSpectraSci2::Reset_Histo";
    for (Int_t i = 0; i < fNbDetectors; i++)
    {
        // === MULTIPLICITY === //
        fh1_multMap[i]->Reset();
        fh1_multTcal[i]->Reset();
        fh1_multMap_condTpat[i]->Reset();
        fh1_multTcal_condTpat[i]->Reset();
        fh2_mult[i]->Reset();
        fh2_mult_RvsL[i]->Reset();
        fh2_mult_RvsL_condTpat[i]->Reset();
        for (Int_t j = 0; j < fNbChannels - 1; j++)
        {
            fh2_mult_TrefVsPmt[i * (fNbChannels - 1) + j]->Reset();
            fh2_mult_TrefVsPmt_condTpat[i * (fNbChannels - 1) + j]->Reset();
        }

        // === FINE TIME === //
        for (Int_t j = 0; j < fNbChannels; j++)
            fh1_finetime[i * fNbChannels + j]->Reset();

        // === POSITION === //
        fh1_RawPos_TcalMult1[i]->Reset();
    }
}

void R3BOnlineSpectraSci2::Exec(Option_t* option)
{

    FairRootManager* mgr = FairRootManager::Instance();
    if (NULL == mgr)
        LOG(fatal) << "R3BOnlineSpectraSci2::Exec FairRootManager not found";

    // --- -------------- --- //
    // --- TPAT CONDITION --- //
    // --- -------------- --- //
    Bool_t BeamOrFission = kTRUE;
    /*if (fEventHeader->GetTpat() > 0)
    {
        if ((fEventHeader->GetTpat() & 0x1) == 1 || // beam
            (fEventHeader->GetTpat() & 0x2) == 2)   // fission

    }*/

    // fTpat = 1-16; fTpat_bit = 0-15
    Int_t fTpat_bit1 = fTpat1 - 1;
    Int_t fTpat_bit2 = fTpat2 - 1;
    Int_t tpatbin;
    if (fTpat1 >= 0 && fTpat2 >= 0)
    {
        for (int i = 0; i < 16; i++)
        {
            tpatbin = (fEventHeader->GetTpat() & (1 << i));
            if (tpatbin != 0 && (i < fTpat_bit1 || i > fTpat_bit2))
            {
                BeamOrFission = kFALSE;
            }
        }
    }

    // --- --------------- --- //
    // --- local variables --- //
    // --- --------------- --- //

    Int_t nHits;
    Int_t iDet; // 0-based
    Int_t iCh;  // 0-based
    Float_t iRawPos;
    Float_t RawPos[fNbDetectors];

    Int_t multMap[fNbDetectors * fNbChannels];
    Int_t multTcal[fNbDetectors * fNbChannels];
    Double_t iRawTimeNs[fNbDetectors * fNbChannels];

    for (Int_t i = 0; i < fNbDetectors; i++)
    {
        for (Int_t j = 0; j < fNbChannels; j++)
        {
            multMap[i * fNbChannels + j] = 0;
            multTcal[i * fNbChannels + j] = 0;
            iRawTimeNs[i * fNbChannels + j] = 0;
        }
        RawPos[i] = -100000.;
    }

    if (fMapped && fMapped->GetEntriesFast() > 0)
    {
        // --- ---------------- --- //
        // --- read mapped data --- //
        // --- ---------------- --- //
        nHits = fMapped->GetEntriesFast();
        for (Int_t ihit = 0; ihit < nHits; ihit++)
        {
            R3BSci2MappedData* hitmapped = dynamic_cast<R3BSci2MappedData*>(fMapped->At(ihit));
            if (!hitmapped)
                continue;
            iDet = hitmapped->GetDetector() - 1;
            iCh = hitmapped->GetChannel() - 1;
            multMap[iDet * fNbChannels + iCh]++;
            fh1_finetime[iDet * fNbChannels + iCh]->Fill(hitmapped->GetTimeFine());
        } // end of loop over mapped data

        if (fTcal && fTcal->GetEntriesFast())
        {
            // --- -------------- --- //
            // --- read tcal data --- //
            // --- -------------- --- //
            nHits = fTcal->GetEntriesFast();
            for (Int_t ihit = 0; ihit < nHits; ihit++)
            {
                R3BSci2TcalData* hittcal = dynamic_cast<R3BSci2TcalData*>(fTcal->At(ihit));
                if (!hittcal)
                    continue;
                iDet = hittcal->GetDetector() - 1;
                iCh = hittcal->GetChannel() - 1;
                multTcal[iDet * fNbChannels + iCh]++;
                iRawTimeNs[iDet * fNbChannels + iCh] = hittcal->GetRawTimeNs();
            } // --- end of loop over Tcal data --- //
        }     //--- end of if Tcal data --- //

        // --- ----------------------------------------- --- //
        // --- filling some histogramms outside the loop --- //
        // --- ----------------------------------------- --- //
        Float_t delta = 0;
        for (Int_t i = 0; i < fNbDetectors; i++)
        {
            fh2_mult_RvsL[i]->Fill(multMap[i * fNbChannels + 1], multMap[i * fNbChannels]);
            for (Int_t j = 0; j < (fNbChannels - 1); j++)
                fh2_mult_TrefVsPmt[i * (fNbChannels - 1) + j]->Fill(multMap[i * fNbChannels + j],
                                                                    multMap[i * fNbChannels + 2]);
            if (BeamOrFission == kTRUE)
            {
                fh2_mult_RvsL_condTpat[i]->Fill(multMap[i * fNbChannels + 1], multMap[i * fNbChannels]);
                for (Int_t j = 0; j < (fNbChannels - 1); j++)
                    fh2_mult_TrefVsPmt_condTpat[i * (fNbChannels - 1) + j]->Fill(multMap[i * fNbChannels + j],
                                                                                 multMap[i * fNbChannels + 2]);
            }

            for (Int_t j = 0; j < fNbChannels; j++)
            {
                fh2_mult[i]->Fill(j + 1, multMap[i * fNbChannels + j]);
                fh1_multMap[i * fNbChannels + j]->Fill(multMap[i * fNbChannels + j]);
                fh1_multTcal[i * fNbChannels + j]->Fill(multTcal[i * fNbChannels + j]);
                if (BeamOrFission == kTRUE)
                {
                    fh1_multMap_condTpat[i * fNbChannels + j]->Fill(multMap[i * fNbChannels + j]);
                    fh1_multTcal_condTpat[i * fNbChannels + j]->Fill(multTcal[i * fNbChannels + j]);
                }
            }
            if ((multTcal[i * fNbChannels] == 1) && (multTcal[i * fNbChannels + 1] == 1))
            {
                // TrawRIGHT-TrawLEFT = 5*(CCr-CCl)+(FTl-FTr) : x is increasing from RIGHT to LEFT
                iRawPos = iRawTimeNs[i * fNbChannels] - iRawTimeNs[i * fNbChannels + 1];
                fh1_RawPos_TcalMult1[i]->Fill(iRawPos);
            }
        }

    } // --- end of if Mapped data --- //

    fNEvents++;
}

// -----   Public method Reset   ------------------------------------------------
void R3BOnlineSpectraSci2::Reset() {}

// -----   Public method Finish   -----------------------------------------------
void R3BOnlineSpectraSci2::FinishEvent()
{
    if (fMapped)
    {
        fMapped->Clear();
    }
    if (fTcal)
    {
        fTcal->Clear();
    }
}

void R3BOnlineSpectraSci2::FinishTask()
{
    for (Int_t i = 0; i < fNbDetectors; i++)
    {
        if (fMapped)
        {
            cMapped[i]->Write();
            fh2_mult[i]->Write();
            fh2_mult_RvsL[i]->Write();
            fh2_mult_RvsL_condTpat[i]->Write();
            for (Int_t j = 0; j < (fNbChannels - 1); j++)
            {
                fh2_mult_TrefVsPmt[i * (fNbChannels - 1) + j]->Write();
                fh2_mult_TrefVsPmt_condTpat[i * (fNbChannels - 1) + j]->Write();
            }
        }
        if (fTcal)
        {
            cPos[i]->Write();
            fh1_RawPos_TcalMult1[i]->Write();
        }
        for (Int_t j = 0; j < fNbChannels; j++)
        {
            if (fMapped)
            {
                fh1_finetime[i * fNbChannels + j]->Write();
                fh1_multMap[i * fNbChannels + j]->Write();
            }
            if (fTcal)
            {
                fh1_multTcal[i * fNbChannels + j]->Write();
            }
        }
    } // end of loop over fNbDetectors
}
ClassImp(R3BOnlineSpectraSci2);
