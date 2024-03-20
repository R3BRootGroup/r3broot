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

// ------------------------------------------------------------
// -----           R3BIncomingIDOnlineSpectra             -----
// -----    Created 05/05/22 by J.L. Rodriguez-Sanchez    -----
// -----           Fill PID online histograms             -----
// ------------------------------------------------------------

#include "R3BIncomingIDOnlineSpectra.h"
#include "R3BEventHeader.h"
#include "R3BFrsData.h"
#include "R3BLogger.h"
#include "R3BLosHitData.h"
#include "R3BMwpcHitData.h"
#include "R3BTGeoPar.h"

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

R3BIncomingIDOnlineSpectra::R3BIncomingIDOnlineSpectra()
    : R3BIncomingIDOnlineSpectra("IncomingIDOnlineSpectra", 1)
{
}

R3BIncomingIDOnlineSpectra::R3BIncomingIDOnlineSpectra(const TString& name, Int_t iVerbose)
    : FairTask(name, iVerbose)
    , fHitFrs(NULL)
    , fHitLos(NULL)
    , fMwpc0HitDataCA(NULL)
    , fMwpc1HitDataCA(NULL)
    , header(nullptr)
    , fNEvents(0)
    , fTpat(-1)
    , fStaId(1)
    , fMin_Z(0.)
    , fMax_Z(20.)
    , fMin_Aq(1.6)
    , fMax_Aq(3.9)
    , fMin_Brho(6)
    , fMax_Brho(10)
    , fMin_RawTof(-50000)
    , fMax_RawTof(50000)
    , fMin_Beta(0.)
    , fMax_Beta(1.)
    , fMin_Z_gate(0.)
    , fMax_Z_gate(20.)
    , fMin_Aq_gate(1.6)
    , fMax_Aq_gate(3.9)
{
}

R3BIncomingIDOnlineSpectra::~R3BIncomingIDOnlineSpectra()
{
    R3BLOG(info, "");
    if (fHitFrs)
        delete fHitFrs;
    if (fHitLos)
        delete fHitLos;
    if (fMwpc0HitDataCA)
        delete fMwpc0HitDataCA;
    if (fMwpc1HitDataCA)
        delete fMwpc1HitDataCA;
}

// -----   Public method SetParContainers   --------------------------------
void R3BIncomingIDOnlineSpectra::SetParContainers()
{
    // Parameter Container
    FairRuntimeDb* rtdb = FairRuntimeDb::instance();
    R3BLOG_IF(fatal, NULL == rtdb, "FairRuntimeDb not found");

    fMw0GeoPar = dynamic_cast<R3BTGeoPar*>(rtdb->getContainer("Mwpc0GeoPar"));
    R3BLOG_IF(error, !fMw0GeoPar, "Could not get access to Mwpc0GeoPar container.");

    fMw1GeoPar = dynamic_cast<R3BTGeoPar*>(rtdb->getContainer("Mwpc1GeoPar"));
    R3BLOG_IF(error, !fMw1GeoPar, "Could not get access to Mwpc1GeoPar container.");
    return;
}

InitStatus R3BIncomingIDOnlineSpectra::Init()
{
    R3BLOG(info, "");
    FairRootManager* mgr = FairRootManager::Instance();
    R3BLOG_IF(fatal, NULL == mgr, "FairRootManager not found");

    FairRunOnline* run = FairRunOnline::Instance();
    run->GetHttpServer()->Register("", this);

    header = dynamic_cast<R3BEventHeader*>(mgr->GetObject("EventHeader."));
    R3BLOG_IF(error, !header, "Branch EventHeader. not found");

    // get access to mapped data of FRS
    fHitFrs = dynamic_cast<TClonesArray*>(mgr->GetObject("FrsData"));
    R3BLOG_IF(fatal, !fHitFrs, "Branch FrsData not found");
    fHitLos = dynamic_cast<TClonesArray*>(mgr->GetObject("LosHit"));
    R3BLOG_IF(fatal, !fHitLos, "Branch LosHitData not found");
    fMwpc0HitDataCA = dynamic_cast<TClonesArray*>(mgr->GetObject("Mwpc0HitData"));
    R3BLOG_IF(fatal, !fMwpc0HitDataCA, "Branch fMwpc0HitDataCA not found");
    fMwpc1HitDataCA = dynamic_cast<TClonesArray*>(mgr->GetObject("Mwpc1HitData"));
    R3BLOG_IF(fatal, !fMwpc1HitDataCA, "Branch fMwpc1HitDataCA not found");

    // Create histograms for detectors
    TString Name1;
    TString Name2;

    cBeta = new TCanvas("Beta_frs", "Beta info", 10, 10, 800, 700);

    // Hit data, beta
    Name1 = "fh1_beta_frs";
    Name2 = "FRS: Beta";
    fh1_beta = new TH1F(Name1, Name2, 1000, fMin_Beta, fMax_Beta);
    fh1_beta->GetXaxis()->SetTitle("Beta");
    fh1_beta->GetYaxis()->SetTitle("Counts");
    fh1_beta->GetYaxis()->SetTitleOffset(1.15);
    fh1_beta->GetXaxis()->CenterTitle(true);
    fh1_beta->GetYaxis()->CenterTitle(true);
    fh1_beta->GetXaxis()->SetLabelSize(0.045);
    fh1_beta->GetXaxis()->SetTitleSize(0.045);
    fh1_beta->GetYaxis()->SetLabelSize(0.045);
    fh1_beta->GetYaxis()->SetTitleSize(0.045);
    fh1_beta->SetFillColor(2);
    fh1_beta->SetLineColor(1);
    fh1_beta->Draw("");

    cTof = new TCanvas("TofRaw_frs", "TofRaw info", 10, 10, 800, 700);

    Name1 = "fh1_tofraw_frs";
    Name2 = "FRS: Raw Tof";
    fh1_tof = new TH1F(Name1, Name2, 20000, fMin_RawTof, fMax_RawTof);
    fh1_tof->GetXaxis()->SetTitle("ToF [ns]");
    fh1_tof->GetYaxis()->SetTitle("Counts");
    fh1_tof->GetYaxis()->SetTitleOffset(1.15);
    fh1_tof->GetXaxis()->CenterTitle(true);
    fh1_tof->GetYaxis()->CenterTitle(true);
    fh1_tof->GetXaxis()->SetLabelSize(0.045);
    fh1_tof->GetXaxis()->SetTitleSize(0.045);
    fh1_tof->GetYaxis()->SetLabelSize(0.045);
    fh1_tof->GetYaxis()->SetTitleSize(0.045);
    fh1_tof->SetFillColor(2);
    fh1_tof->SetLineColor(1);
    fh1_tof->Draw("");

    // Hit data, brho
    cBrho = new TCanvas("Brho_frs", "Brho info", 10, 10, 800, 700);

    Name1 = "fh1_brho_frs";
    Name2 = "FRS: Brho S2-Cave";
    fh1_brho = new TH1F(Name1, Name2, 1500, fMin_Brho, fMax_Brho);
    fh1_brho->GetXaxis()->SetTitle("Brho [Tm]");
    fh1_brho->GetYaxis()->SetTitle("Counts");
    fh1_brho->GetYaxis()->SetTitleOffset(1.15);
    fh1_brho->GetXaxis()->CenterTitle(true);
    fh1_brho->GetYaxis()->CenterTitle(true);
    fh1_brho->GetXaxis()->SetLabelSize(0.045);
    fh1_brho->GetXaxis()->SetTitleSize(0.045);
    fh1_brho->GetYaxis()->SetLabelSize(0.045);
    fh1_brho->GetYaxis()->SetTitleSize(0.045);
    fh1_brho->SetFillColor(2);
    fh1_brho->SetLineColor(1);
    fh1_brho->Draw("");

    // Hit data, Xs2_vs_Beta
    cXs2vsBeta = new TCanvas("Xs2_vs_beta_frs", "Xs2_vs_Beta 2D info", 10, 10, 800, 700);

    Name1 = "fh2_Xs2_vs_beta_frs";
    Name2 = "FRS: Xs2 vs #beta";
    fh2_Xs2vsbeta = new TH2F(Name1, Name2, 800, -100., 100., 1000, 0.45, 0.99);
    fh2_Xs2vsbeta->GetXaxis()->SetTitle("X at S2 [mm]");
    fh2_Xs2vsbeta->GetYaxis()->SetTitle("FRS-#beta");
    fh2_Xs2vsbeta->GetYaxis()->SetTitleOffset(1.1);
    fh2_Xs2vsbeta->GetXaxis()->CenterTitle(true);
    fh2_Xs2vsbeta->GetYaxis()->CenterTitle(true);
    fh2_Xs2vsbeta->GetXaxis()->SetLabelSize(0.045);
    fh2_Xs2vsbeta->GetXaxis()->SetTitleSize(0.045);
    fh2_Xs2vsbeta->GetYaxis()->SetLabelSize(0.045);
    fh2_Xs2vsbeta->GetYaxis()->SetTitleSize(0.045);
    fh2_Xs2vsbeta->Draw("col");

    cAoQvsPosS2 = new TCanvas("AoQvsPosS2", "Pos-S2 vs AoQ", 10, 10, 800, 700);

    Name1 = "fh2_PosS2vsAoQ_frs";
    Name2 = "FRS: Pos-S2 vs AoQ with mult==1";
    fh2_Pos2vsAoQ_m1 = new TH2F(Name1, Name2, 2000, -100, 100, 1000, fMin_Aq, fMax_Aq);
    fh2_Pos2vsAoQ_m1->GetXaxis()->SetTitle("Pos-S2 [mm]");
    fh2_Pos2vsAoQ_m1->GetYaxis()->SetTitle("AoQ");
    fh2_Pos2vsAoQ_m1->GetYaxis()->SetTitleOffset(1.1);
    fh2_Pos2vsAoQ_m1->GetXaxis()->CenterTitle(true);
    fh2_Pos2vsAoQ_m1->GetYaxis()->CenterTitle(true);
    fh2_Pos2vsAoQ_m1->GetXaxis()->SetLabelSize(0.045);
    fh2_Pos2vsAoQ_m1->GetXaxis()->SetTitleSize(0.045);
    fh2_Pos2vsAoQ_m1->GetYaxis()->SetLabelSize(0.045);
    fh2_Pos2vsAoQ_m1->GetYaxis()->SetTitleSize(0.045);
    fh2_Pos2vsAoQ_m1->Draw("colz");

    // Hit data, Aq_vs_q
    cAqvsq = new TCanvas("Aq_vs_q_frs", "A/q_vs_q 2D info", 10, 10, 800, 700);

    Name1 = "fh2_Aq_vs_q_frs";
    Name2 = "FRS: A/q vs q";
    fh2_Aqvsq = new TH2F(Name1, Name2, 1000, fMin_Aq, fMax_Aq, 1900, fMin_Z, fMax_Z);
    fh2_Aqvsq->GetXaxis()->SetTitle("A/q");
    fh2_Aqvsq->GetYaxis()->SetTitle("Z [Charge units]");
    fh2_Aqvsq->GetYaxis()->SetTitleOffset(1.1);
    fh2_Aqvsq->GetXaxis()->CenterTitle(true);
    fh2_Aqvsq->GetYaxis()->CenterTitle(true);
    fh2_Aqvsq->GetXaxis()->SetLabelSize(0.045);
    fh2_Aqvsq->GetXaxis()->SetTitleSize(0.045);
    fh2_Aqvsq->GetYaxis()->SetLabelSize(0.045);
    fh2_Aqvsq->GetYaxis()->SetTitleSize(0.045);
    fh2_Aqvsq->Draw("colz");

    // Isotope gated ion images
    TString RangeGates =
        Form("IsoGated_Z_%1.1f_%1.1f_Aq_%1.2f_%1.2f", fMin_Z_gate, fMax_Z_gate, fMin_Aq_gate, fMax_Aq_gate);
    cIsoGated = new TCanvas("Isotope_gated_images", "Isotope gated ion images", 10, 10, 800, 700);
    cIsoGated->Divide(2, 2);

    cIsoGated->cd(1);
    Name1 = "Z_xcave";
    Name2 = "Z vs Xcave";
    fh2_Z_xc = new TH2F(Name1, Name2, 1000, -100, 100, 1000, fMin_Z, fMax_Z);
    fh2_Z_xc->GetXaxis()->SetTitle("XCave [mm]");
    fh2_Z_xc->GetYaxis()->SetTitle("Z [charge units]");
    fh2_Z_xc->GetYaxis()->SetTitleOffset(1.1);
    fh2_Z_xc->GetXaxis()->CenterTitle(true);
    fh2_Z_xc->GetYaxis()->CenterTitle(true);
    fh2_Z_xc->GetXaxis()->SetLabelSize(0.045);
    fh2_Z_xc->GetXaxis()->SetTitleSize(0.045);
    fh2_Z_xc->GetYaxis()->SetLabelSize(0.045);
    fh2_Z_xc->GetYaxis()->SetTitleSize(0.045);
    fh2_Z_xc->Draw("colz");

    cIsoGated->cd(2);
    Name1 = RangeGates + "_Z_xcave";
    Name2 = RangeGates + " Z vs Xcave";
    fh2_IsoGated_Z_xc = new TH2F(Name1, Name2, 1000, -100, 100, 1000, fMin_Z, fMax_Z);
    fh2_IsoGated_Z_xc->GetXaxis()->SetTitle("XCave [mm]");
    fh2_IsoGated_Z_xc->GetYaxis()->SetTitle("Z [charge units]");
    fh2_IsoGated_Z_xc->GetYaxis()->SetTitleOffset(1.1);
    fh2_IsoGated_Z_xc->GetXaxis()->CenterTitle(true);
    fh2_IsoGated_Z_xc->GetYaxis()->CenterTitle(true);
    fh2_IsoGated_Z_xc->GetXaxis()->SetLabelSize(0.045);
    fh2_IsoGated_Z_xc->GetXaxis()->SetTitleSize(0.045);
    fh2_IsoGated_Z_xc->GetYaxis()->SetLabelSize(0.045);
    fh2_IsoGated_Z_xc->GetYaxis()->SetTitleSize(0.045);
    fh2_IsoGated_Z_xc->Draw("colz");

    cIsoGated->cd(3);
    Name1 = RangeGates + "_xs2_xcave";
    Name2 = RangeGates + " Xs2 vs Xcave";
    fh2_IsoGated_xs2_xc = new TH2F(Name1, Name2, 1000, -100, 100, 1000, -100, 100);
    fh2_IsoGated_xs2_xc->GetXaxis()->SetTitle("Xs2 [mm]");
    fh2_IsoGated_xs2_xc->GetYaxis()->SetTitle("Xcave [mm]");
    fh2_IsoGated_xs2_xc->GetYaxis()->SetTitleOffset(1.1);
    fh2_IsoGated_xs2_xc->GetXaxis()->CenterTitle(true);
    fh2_IsoGated_xs2_xc->GetYaxis()->CenterTitle(true);
    fh2_IsoGated_xs2_xc->GetXaxis()->SetLabelSize(0.045);
    fh2_IsoGated_xs2_xc->GetXaxis()->SetTitleSize(0.045);
    fh2_IsoGated_xs2_xc->GetYaxis()->SetLabelSize(0.045);
    fh2_IsoGated_xs2_xc->GetYaxis()->SetTitleSize(0.045);
    fh2_IsoGated_xs2_xc->Draw("colz");

    cIsoGated->cd(4);
    Name1 = RangeGates + "_xcave_acave";
    Name2 = RangeGates + " Xcave vs Anglecave";
    fh2_IsoGated_xc_anglec = new TH2F(Name1, Name2, 1000, -100, 100, 1000, -100, 100);
    fh2_IsoGated_xc_anglec->GetXaxis()->SetTitle("Xcave [mm]");
    fh2_IsoGated_xc_anglec->GetYaxis()->SetTitle("Angle Cave [mrad]");
    fh2_IsoGated_xc_anglec->GetYaxis()->SetTitleOffset(1.1);
    fh2_IsoGated_xc_anglec->GetXaxis()->CenterTitle(true);
    fh2_IsoGated_xc_anglec->GetYaxis()->CenterTitle(true);
    fh2_IsoGated_xc_anglec->GetXaxis()->SetLabelSize(0.045);
    fh2_IsoGated_xc_anglec->GetXaxis()->SetTitleSize(0.045);
    fh2_IsoGated_xc_anglec->GetYaxis()->SetLabelSize(0.045);
    fh2_IsoGated_xc_anglec->GetYaxis()->SetTitleSize(0.045);
    fh2_IsoGated_xc_anglec->Draw("colz");

    cLosE_Tof = new TCanvas("LosE_Tof", "DeltaE in LOS and Tof", 10, 10, 800, 700);

    Name1 = "LOS-E_vs_ToF";
    Name2 = "LOS Energy vs Raw ToF (S2-LOS);Raw ToF (S2-LOS) / ns;LOS Z";
    fh2_LosE_Tof = new TH2F(Name1, Name2, 3000, 0, 3000, 900, 0, 30);
    fh2_LosE_Tof->GetYaxis()->SetTitleOffset(1.1);
    fh2_LosE_Tof->GetXaxis()->CenterTitle(true);
    fh2_LosE_Tof->GetYaxis()->CenterTitle(true);
    fh2_LosE_Tof->GetXaxis()->SetLabelSize(0.045);
    fh2_LosE_Tof->GetXaxis()->SetTitleSize(0.045);
    fh2_LosE_Tof->GetYaxis()->SetLabelSize(0.045);
    fh2_LosE_Tof->GetYaxis()->SetTitleSize(0.045);
    fh2_LosE_Tof->Draw("colz");

    cLosE_Tof2 = new TCanvas(
        "LosE_with_ToF_measurement", "LOS Energy distribution with (black) and without (red) S2", 10, 10, 800, 700);

    Name1 = "LOS-Z_without_S2";
    Name2 = "LOS Energy with (black) without (red) good S2 hit;LOS Z;Counts";
    fh1_LosE_withoutTof = new TH1F(Name1, Name2, 900, 0, 30);
    fh1_LosE_withoutTof->GetYaxis()->SetTitleOffset(1.1);
    fh1_LosE_withoutTof->GetXaxis()->CenterTitle(true);
    fh1_LosE_withoutTof->GetYaxis()->CenterTitle(true);
    fh1_LosE_withoutTof->GetXaxis()->SetLabelSize(0.045);
    fh1_LosE_withoutTof->GetXaxis()->SetTitleSize(0.045);
    fh1_LosE_withoutTof->GetYaxis()->SetLabelSize(0.045);
    fh1_LosE_withoutTof->GetYaxis()->SetTitleSize(0.045);
    fh1_LosE_withoutTof->SetLineColor(kRed);
    fh1_LosE_withoutTof->SetLineWidth(2);
    fh1_LosE_withoutTof->Draw("");

    Name1 = "LOS-Z_with_S2";
    Name2 = "LOS Energy with (black) without (red) good S2 hit;LOS Z;Counts";
    fh1_LosE_withTof = new TH1F(Name1, Name2, 900, 0, 30);
    fh1_LosE_withTof->GetYaxis()->SetTitleOffset(1.1);
    fh1_LosE_withTof->GetXaxis()->CenterTitle(true);
    fh1_LosE_withTof->GetYaxis()->CenterTitle(true);
    fh1_LosE_withTof->GetXaxis()->SetLabelSize(0.045);
    fh1_LosE_withTof->GetXaxis()->SetTitleSize(0.045);
    fh1_LosE_withTof->GetYaxis()->SetLabelSize(0.045);
    fh1_LosE_withTof->GetYaxis()->SetTitleSize(0.045);
    fh1_LosE_withTof->Draw("same");

    // MAIN FOLDER-INCOMINGID
    TFolder* mainfol = new TFolder("FRS-IncomingID", "FRS incomingID info");
    mainfol->Add(cAoQvsPosS2);
    mainfol->Add(cBeta);
    mainfol->Add(cTof);
    mainfol->Add(cBrho);
    mainfol->Add(cXs2vsBeta);
    mainfol->Add(cAqvsq);
    mainfol->Add(cIsoGated);
    mainfol->Add(cLosE_Tof);
    mainfol->Add(cLosE_Tof2);
    run->AddObject(mainfol);

    // Register command to reset histograms
    run->GetHttpServer()->RegisterCommand("Reset_IncomingID_HIST", Form("/Objects/%s/->Reset_Histo()", GetName()));

    return kSUCCESS;
}

void R3BIncomingIDOnlineSpectra::Reset_Histo()
{
    R3BLOG(info, "");
    fh2_Pos2vsAoQ_m1->Reset();
    fh1_beta->Reset();
    fh1_tof->Reset();
    fh1_brho->Reset();
    fh2_Aqvsq->Reset();
    fh2_Xs2vsbeta->Reset();
    fh2_Z_xc->Reset();
    fh2_IsoGated_Z_xc->Reset();
    fh2_IsoGated_xs2_xc->Reset();
    fh2_IsoGated_xc_anglec->Reset();
    fh2_LosE_Tof->Reset();
    fh1_LosE_withTof->Reset();
    fh1_LosE_withoutTof->Reset();
}

void R3BIncomingIDOnlineSpectra::Exec(Option_t* option)
{
    if ((fTpat >= 0) && (header) && ((header->GetTpat() & fTpat) != fTpat))
        return;

    // Fill Hit data
    if (fHitFrs && fHitFrs->GetEntriesFast() > 0)
    {
        Int_t nHits = fHitFrs->GetEntriesFast();
        for (Int_t ihit = 0; ihit < nHits; ihit++)
        {
            R3BFrsData* hit = dynamic_cast<R3BFrsData*>(fHitFrs->At(ihit));
            if (!hit)
                continue;
            if (hit->GetStaId() != fStaId)
                continue;
            fh2_Pos2vsAoQ_m1->Fill(hit->GetXS2(), hit->GetAq());
            fh1_tof->Fill(hit->GetTof());
            fh1_beta->Fill(hit->GetBeta());
            fh1_brho->Fill(hit->GetBrho());
            fh2_Aqvsq->Fill(hit->GetAq(), hit->GetZ());
            fh2_Xs2vsbeta->Fill(hit->GetXS2(), hit->GetBeta());

            auto nHits_Mw0 = fMwpc0HitDataCA->GetEntriesFast();
            auto nHits_Mw1 = fMwpc1HitDataCA->GetEntriesFast();
            for (Int_t iMw0 = 0; iMw0 < nHits_Mw0; iMw0++)
            {
                auto hit_mw0 = dynamic_cast<R3BMwpcHitData*>(fMwpc0HitDataCA->At(iMw0));
                if (!hit_mw0)
                    continue;
                auto mwpc0x = hit_mw0->GetX() + fMw0GeoPar->GetPosX() * 10.; // mm
                for (Int_t iMw1 = 0; iMw1 < nHits_Mw1; iMw1++)
                {
                    auto hit_mw1 = dynamic_cast<R3BMwpcHitData*>(fMwpc1HitDataCA->At(iMw1));
                    if (!hit_mw1)
                        continue;
                    auto mwpc1x = hit_mw1->GetX() + fMw1GeoPar->GetPosX() * 10.; // mm
                    auto XCave = mwpc0x;
                    auto AngleCave =
                        (mwpc0x - mwpc1x) / (fMw0GeoPar->GetPosZ() - fMw1GeoPar->GetPosZ()) / 10. * 1000.; // mrad

                    fh2_Z_xc->Fill(XCave, hit->GetZ());
                    // Plot PID gated histograms below
                    if (hit->GetAq() < fMin_Aq_gate || hit->GetAq() > fMax_Aq_gate || hit->GetZ() < fMin_Z_gate ||
                        hit->GetZ() > fMax_Z_gate)
                        continue;
                    fh2_IsoGated_Z_xc->Fill(XCave, hit->GetZ());
                    fh2_IsoGated_xs2_xc->Fill(hit->GetXS2(), XCave);
                    fh2_IsoGated_xc_anglec->Fill(XCave, AngleCave);
                }
            }
            // making los pid
            if (fHitLos && fHitLos->GetEntriesFast() > 0)
            {
                Int_t nHitsLos = fHitLos->GetEntriesFast();
                for (Int_t ihitLos = 0; ihitLos < nHitsLos; ihitLos++)
                {
                    R3BLosHitData* hitLos = dynamic_cast<R3BLosHitData*>(fHitLos->At(ihitLos));
                    if (!hitLos)
                        continue;
                    fh2_LosE_Tof->Fill(hit->GetTof(), hitLos->GetZ());
                    fh1_LosE_withTof->Fill(hitLos->GetZ());
                }
            }
        } // for hit frs
    }
    else // for the case of no ToF in FRS
    {
        if (fHitLos && fHitLos->GetEntriesFast() > 0)
        {
            Int_t nHitsLos = fHitLos->GetEntriesFast();
            for (Int_t ihitLos = 0; ihitLos < nHitsLos; ihitLos++)
            {
                R3BLosHitData* hitLos = dynamic_cast<R3BLosHitData*>(fHitLos->At(ihitLos));
                if (!hitLos)
                    continue;
                fh1_LosE_withoutTof->Fill(hitLos->GetZ());
            }
        }
    }
    fNEvents += 1;
}

void R3BIncomingIDOnlineSpectra::FinishEvent()
{
    if (fHitFrs)
    {
        fHitFrs->Clear();
    }

    if (fHitLos)
    {
        fHitLos->Clear();
    }

    if (fMwpc0HitDataCA)
    {
        fMwpc0HitDataCA->Clear();
    }

    if (fMwpc1HitDataCA)
    {
        fMwpc1HitDataCA->Clear();
    }
}

void R3BIncomingIDOnlineSpectra::FinishTask()
{
    if (fHitFrs && fMwpc0HitDataCA && fMwpc1HitDataCA)
    {
        cAoQvsPosS2->Write();
        cBeta->Write();
        cTof->Write();
        cBrho->Write();
        cXs2vsBeta->Write();
        cAqvsq->Write();
        cIsoGated->Write();
        cLosE_Tof->Write();
        cLosE_Tof2->Write();
    }
}

ClassImp(R3BIncomingIDOnlineSpectra)
