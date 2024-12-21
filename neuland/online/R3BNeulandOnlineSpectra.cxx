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

#include "R3BNeulandOnlineSpectra.h"
#include "R3BEventHeader.h"
#include "R3BShared.h"

#include <FairRunOnline.h>

#include <FairRootManager.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <THttpServer.h>
#include <iostream>
#include <limits>

R3BNeulandOnlineSpectra::R3BNeulandOnlineSpectra()
    : FairTask("R3BNeulandOnlineSpectra", 0)
    , fNeulandMappedData("NeulandMappedData")
    , fNeulandCalData("NeulandCalData")
    , fNeulandHits("NeulandHits")
{
}

InitStatus R3BNeulandOnlineSpectra::Init()
{
    // Initialize random number:
    std::srand(std::time(0)); // use current time as seed for random generator

    auto run = FairRunOnline::Instance();
    if (run != nullptr && run->GetHttpServer() != nullptr)
    {
        fIsOnline = true;
    }

    if (fIsOnline)
    {
        run->GetHttpServer()->Register("/Tasks", this);
        run->GetHttpServer()->RegisterCommand("Reset_Neuland", Form("/Tasks/%s/->ResetHistos()", GetName()));
        run->GetHttpServer()->RegisterCommand("Reset_Neuland_Mapped",
                                              Form("/Tasks/%s/->ResetHistosMapped()", GetName()));
    }

    auto ioman = FairRootManager::Instance();
    if (ioman == nullptr)
    {
        throw std::runtime_error("R3BNeulandOnlineSpectra: No FairRootManager");
    }

    fEventHeader = dynamic_cast<R3BEventHeader*>(ioman->GetObject("EventHeader."));
    if (fEventHeader == nullptr)
    {
        throw std::runtime_error("R3BNeulandOnlineSpectra: No R3BEventHeader");
    }

    fNeulandMappedData.Init();
    fNeulandCalData.Init();
    fNeulandHits.Init();

    auto canvasMapped = new TCanvas("NeulandMapped", "NeulandMapped", 10, 10, 850, 850);
    canvasMapped->Divide(1, 2);

    canvasMapped->cd(1);
    ahMappedBar1[0] = R3B::root_owned<TH1D>("hMappedBar1fLE", "Mapped: Bars fine 1LE", fNBars, 0.5, fNBars + 0.5);
    ahMappedBar1[0]->Draw();
    ahMappedBar1[1] = R3B::root_owned<TH1D>("hMappedBar1fTE", "Mapped: Bars fine 1TE", fNBars, 0.5, fNBars + 0.5);
    ahMappedBar1[1]->SetLineColor(1);
    ahMappedBar1[1]->Draw("same");
    ahMappedBar1[2] = R3B::root_owned<TH1D>("hMappedBar1cLE", "Mapped: Bars coarse 1LE", fNBars, 0.5, fNBars + 0.5);
    ahMappedBar1[2]->SetLineColor(2);
    ahMappedBar1[2]->Draw("same");
    ahMappedBar1[3] = R3B::root_owned<TH1D>("hMappedBar1cTE", "Mapped: Bars coarse 1TE", fNBars, 0.5, fNBars + 0.5);
    ahMappedBar1[3]->SetLineColor(6);
    ahMappedBar1[3]->Draw("same");

    canvasMapped->cd(2);
    ahMappedBar2[0] = R3B::root_owned<TH1D>("hMappedBar2fLE", "Mapped: Bars fine 2LE", fNBars, 0.5, fNBars + 0.5);
    ahMappedBar2[0]->Draw();
    ahMappedBar2[1] = R3B::root_owned<TH1D>("hMappedBar2fTE", "Mapped: Bars fine 2TE", fNBars, 0.5, fNBars + 0.5);
    ahMappedBar2[1]->SetLineColor(1);
    ahMappedBar2[1]->Draw("same");
    ahMappedBar2[2] = R3B::root_owned<TH1D>("hMappedBar2cLE", "Mapped: Bars coarse 2LE", fNBars, 0.5, fNBars + 0.5);
    ahMappedBar2[2]->SetLineColor(2);
    ahMappedBar2[2]->Draw("same");
    ahMappedBar2[3] = R3B::root_owned<TH1D>("hMappedBar2cTE", "Mapped: Bars coarse 2TE", fNBars, 0.5, fNBars + 0.5);
    ahMappedBar2[3]->SetLineColor(6);
    ahMappedBar2[3]->Draw("same");

    canvasMapped->cd(0);
    if (fIsOnline)
    {
        run->AddObject(canvasMapped);
    }

    auto canvasJump = new TCanvas("NeulandJumps", "NeulandJumps", 10, 10, 850, 850);
    canvasJump->Divide(1, 2);

    hJumpsvsEvnt = R3B::root_owned<TH2D>("hJumpsvsEvnt", "Jumps vs Evnt", 1000, 0, 10000000, 1000, -11000, 11000);
    canvasJump->cd(1);
    gPad->SetLogz();
    hJumpsvsEvnt->Draw("colz");

    hJumpsvsEvntzoom =
        R3B::root_owned<TH2D>("hJumpsvsEvntzoom", "Jumps vs Evnt zoomed", 1000, 0, 10000000, 1000, -200, 200);
    canvasJump->cd(2);
    gPad->SetLogz();
    hJumpsvsEvntzoom->Draw("colz");

    canvasJump->cd(0);
    if (fIsOnline)
    {
        run->AddObject(canvasJump);
    }

    auto canvasCal = new TCanvas("NeulandCal", "NeulandCal", 10, 10, 850, 850);
    canvasCal->Divide(2, 3);

    hTstart = R3B::root_owned<TH1D>("hTstart", "Tstart", 1000, -10000., 50000.);
    canvasCal->cd(1);
    hTstart->Draw();

    hNstart = R3B::root_owned<TH1D>("hNstart", "Nstart", 20, 0.5, 20.5);
    hTestJump = R3B::root_owned<TH2D>("hTestJump", "Test Jump", fNBars, 0.5, fNBars + 0.5, 1000, -11000., 11000.);
    canvasCal->cd(2);
    // hNstart->Draw();
    gPad->SetLogz();
    hTestJump->Draw("colz");

    ahCalTvsBar[0] = R3B::root_owned<TH2D>(
        "hCalT1vsBar", "CalLevel: Time1 vs Bars  ", fNBars, 0.5, fNBars + 0.5, 2000, -5000, 15000);
    canvasCal->cd(3);
    ahCalTvsBar[0]->Draw("colz");

    ahCalTvsBar[1] = R3B::root_owned<TH2D>(
        "hCalT2vsBar", "CalLevel: Time2 vs Bars  ", fNBars, 0.5, fNBars + 0.5, 2000, -5000, 15000);
    canvasCal->cd(4);
    ahCalTvsBar[1]->Draw("colz");

    ahCalEvsBar[0] =
        R3B::root_owned<TH2D>("hCalE1vsBar", "CalLevel: Energy1 vs Bars", fNBars, 0.5, fNBars + 0.5, 600, 0, 600);
    canvasCal->cd(5);
    ahCalEvsBar[0]->Draw("colz");

    ahCalEvsBar[1] =
        R3B::root_owned<TH2D>("hCalE2vsBar", "CalLevel: Energy2 vs Bars", fNBars, 0.5, fNBars + 0.5, 600, 0, 600);
    canvasCal->cd(6);
    ahCalEvsBar[1]->Draw("colz");

    canvasCal->cd(0);
    if (fIsOnline)
    {
        run->AddObject(canvasCal);
    }

    auto canvasHit = new TCanvas("NeulandHit", "NeulandHit", 10, 10, 850, 850);
    canvasHit->Divide(2, 2);

    hHitEvsBar =
        R3B::root_owned<TH2D>("hHitEvsBar", "HitLevel: Energy vs Bars ", fNBars, 0.5, fNBars + 0.5, 1000, 0, 60);
    canvasHit->cd(1);
    hHitEvsBar->Draw("colz");

    hTdiffvsBar = R3B::root_owned<TH2D>("hTdiffvsBar", "Tdiff vs Bars", fNBars, 0.5, fNBars + 0.5, 1000, -60, 60);
    canvasHit->cd(2);
    hTdiffvsBar->Draw("colz");

    hToFvsBar = R3B::root_owned<TH2D>("hTofvsBar", "Tof vs Bars", fNBars, 0.5, fNBars + 0.5, 6000, -100, 400);
    hToFcvsBar = R3B::root_owned<TH2D>("hTofcvsBar", "Tofc vs Bars", fNBars, 0.5, fNBars + 0.5, 6000, -100, 400);
    canvasHit->cd(3);
    // hToFcvsBar->Draw("colz");
    hToFvsBar->Draw("colz");

    hTofvsEhit = R3B::root_owned<TH2D>("hTofvsEhit", "Tof vs Ehit", 2000, 0, 120, 6000, -100, 400);
    hTofcvsEhit = R3B::root_owned<TH2D>("hTofcvsEhit", "Tofc vs Ehit", 2000, 0, 120, 6000, -100, 400);
    canvasHit->cd(4);
    // hTofcvsEhit->Draw("colz");
    hTofvsEhit->Draw("colz");

    canvasHit->cd(0);
    if (fIsOnline)
    {
        run->AddObject(canvasHit);
    }

    auto canvasHitCosmics = new TCanvas("NeulandHitCosmics", "NeulandHitCosmics", 10, 10, 850, 850);
    canvasHitCosmics->Divide(2, 2);

    hHitEvsBarCosmics = R3B::root_owned<TH2D>(
        "hHitEvsBarCosmics", "HitLevel: Energy vs Bars cosmics", fNBars, 0.5, fNBars + 0.5, 1000, 0, 60);
    canvasHitCosmics->cd(1);
    hHitEvsBarCosmics->Draw("colz");

    hTdiffvsBarCosmics =
        R3B::root_owned<TH2D>("hTdiffvsBarCosmics", "Tdiff vs Bars cosmics", fNBars, 0.5, fNBars + 0.5, 1000, -60, 60);
    canvasHitCosmics->cd(2);
    hTdiffvsBarCosmics->Draw("colz");

    hDT675 =
        R3B::root_owned<TH2D>("hDT675", "Thit - Thit675 vs Bars cosmics", fNBars, 0.5, fNBars + 0.5, 1000, -20, 20);
    canvasHitCosmics->cd(3);
    hDT675->Draw("colz");

    hDT675c = R3B::root_owned<TH2D>(
        "hDT675c", "Thit - Thit675 vs Bars cosmics corrected", fNBars, 0.5, fNBars + 0.5, 1000, -20, 20);

    hDT625 =
        R3B::root_owned<TH2D>("hDT625", "Thit - Thit625 vs Bars cosmics", fNBars, 0.5, fNBars + 0.5, 1000, -20, 20);
    canvasHitCosmics->cd(4);
    hDT625->Draw("colz");

    hDT625c = R3B::root_owned<TH2D>(
        "hDT625c", "Thit - Thit625 vs Bars cosmics corrected", fNBars, 0.5, fNBars + 0.5, 1000, -20, 20);

    canvasHitCosmics->cd(0);
    if (fIsOnline)
    {
        run->AddObject(canvasHitCosmics);
    }

    auto canvasPlaneXY = new TCanvas("NeulandPlaneXY", "NeulandPlaneXY", 10, 10, 850, 850);
    canvasPlaneXY->Divide(5, 6);
    for (unsigned int i = 0; i < fNPlanes; i++)
    {
        ahXYperPlane[i] = R3B::root_owned<TH2D>("hHitXYPlane" + TString::Itoa(i, 10),
                                                "Hit XY Plane" + TString::Itoa(i, 10),
                                                300,
                                                -150,
                                                150,
                                                300,
                                                -150,
                                                150);
        canvasPlaneXY->cd(i + 1);
        ahXYperPlane[i]->Draw("colz");
    }
    canvasPlaneXY->cd(0);
    if (fIsOnline)
    {
        run->AddObject(canvasPlaneXY);
    }

    auto canvasPlaneSofia = new TCanvas("Timing", "Timing", 10, 10, 850, 850);
    canvasPlaneSofia->Divide(2, 2);
    hTofvsX = R3B::root_owned<TH2D>("hTofvsX", "Tof vs X", 3000, -200., 200., 6000, 0, 300);
    hTofcvsX = R3B::root_owned<TH2D>("hTofcvsX", "Tofc vs X", 1000, -200., 200., 3000, 0, 300);
    hTofvsY = R3B::root_owned<TH2D>("hTofvsY", "Tof vs Y", 3000, -200., 200., 6000, 0, 300);
    hTofcvsY = R3B::root_owned<TH2D>("hTofcvsY", "Tofc vs Y", 3000, -200., 200., 6000, 0, 300);
    hTofvsZ = R3B::root_owned<TH2D>("hTofvsZ", "Tof vs Z", 26, 0, 26, 3000, 0, 300);
    hTofcvsZ = R3B::root_owned<TH2D>("hTofcvsZ", "Tofc vs Z", 26, 0, 26, 3000, 0, 300);

    hSofiaTime = R3B::root_owned<TH1D>("hSofiaTime", "hSofiaTime", 50000, 0, 50000);
    // canvasPlaneSofia->cd(1);
    // hSofiaTime->Draw();

    hNeuLANDvsStart =
        R3B::root_owned<TH2D>("hNeuLANDvsStart", "hNeuLANDvsStart", 3000, -10000, 40000, 1000, -10000, 10000);
    canvasPlaneSofia->cd(1);
    hNeuLANDvsStart->Draw("colz");

    hTOF = R3B::root_owned<TH1D>("hTOF", "hTOF", 6000, -11000, 11000);
    hTOFc = R3B::root_owned<TH1D>("hTOFc", "hTOFc", 6000, -100, 500);
    canvasPlaneSofia->cd(2);
    gPad->SetLogy();
    // hTOF->Draw();
    // hTOFc->SetLineColor(2);
    hTOFc->Draw("");

    canvasPlaneSofia->cd(3);
    gPad->SetLogz();
    hTofcvsX->Draw("colz");
    // hTofvsX->Draw("colz");

    canvasPlaneSofia->cd(4);
    gPad->SetLogz();
    // hTofcvsZ->Draw("colz");
    hTofvsZ->Draw("colz");

    canvasPlaneSofia->cd(0);
    if (fIsOnline)
    {
        run->AddObject(canvasPlaneSofia);
    }

    return kSUCCESS;
}

void R3BNeulandOnlineSpectra::Exec(Option_t*)
{
    const double clight = 29.9792458;

    const double start = fEventHeader->GetTStart();

    if (((UInt_t)fEventHeader->GetEventno() % 10000000) < 10000)
    {
        hJumpsvsEvnt->Reset();
        hJumpsvsEvntzoom->Reset();
    }

    const auto mappedData = fNeulandMappedData.Retrieve();
    const auto calData = fNeulandCalData.Retrieve();
    const auto hits = fNeulandHits.Retrieve();

    for (const auto& mapped : mappedData)
    {
        const auto plane = mapped->GetPlaneId();
        const auto barp = mapped->GetBarId();
        const auto bar = (plane - 1) * 50 + barp;

        if (mapped->GetFineTime1LE() > 0)
            ahMappedBar1[0]->Fill(bar);
        if (mapped->GetFineTime1TE() > 0)
            ahMappedBar1[1]->Fill(bar);
        if (mapped->GetCoarseTime1LE() > 0)
            ahMappedBar1[2]->Fill(bar);
        if (mapped->GetCoarseTime1TE() > 0)
            ahMappedBar1[3]->Fill(bar);
        if (mapped->GetFineTime2LE() > 0)
            ahMappedBar2[0]->Fill(bar);
        if (mapped->GetFineTime2TE() > 0)
            ahMappedBar2[1]->Fill(bar);
        if (mapped->GetCoarseTime2LE() > 0)
            ahMappedBar2[2]->Fill(bar);
        if (mapped->GetCoarseTime2TE() > 0)
            ahMappedBar2[3]->Fill(bar);
    }

    for (const auto& data : calData)
    {
        const auto side = data->GetSide() - 1; // [1,2] -> [0,1]
        const auto bar = data->GetBarId();
        ahCalTvsBar[side]->Fill(bar, data->GetTime());
        ahCalEvsBar[side]->Fill(bar, data->GetQdc());
        if (std::isnan(data->GetTriggerTime()))
        {
            hNeuLANDvsStart->Fill(start, data->GetTime());
        }
        else
        {
            hNeuLANDvsStart->Fill(start, data->GetTime() - data->GetTriggerTime());
        }

        for (const auto& datax : calData)
        {
            const auto sidex = datax->GetSide() - 1; // [1,2] -> [0,1]
            const auto barx = datax->GetBarId();

            if (barx != bar)
            {
                hTestJump->Fill(barx, data->GetTime() - datax->GetTime());
                hJumpsvsEvnt->Fill((UInt_t)fEventHeader->GetEventno() % 10000000, data->GetTime() - datax->GetTime());
                hJumpsvsEvntzoom->Fill((UInt_t)fEventHeader->GetEventno() % 10000000,
                                       data->GetTime() - datax->GetTime());
            }
        }
    }

    Double_t randx;

    for (const auto& hit : hits)
    {
        const auto bar = hit->GetPaddle();

        if (IsBeam())
        {
            hTstart->Fill(start);
            if (std::isnan(hit->GetT()))
                continue;
            hHitEvsBar->Fill(bar, hit->GetE());
            hTdiffvsBar->Fill(bar, hit->GetTdcL() - hit->GetTdcR());

            const Double_t tcorr = hit->GetT() - (hit->GetPosition().Mag() - fDistanceToTarget) / clight;
            const Double_t tadj = hit->GetT();

            if (hit->GetE() > 0.)
            { // 7.
                hToFvsBar->Fill(bar, tadj);
                hToFcvsBar->Fill(bar, tcorr);
                hTOFc->Fill(tcorr);
                hTOF->Fill(tadj);
            }

            hTofvsEhit->Fill(hit->GetE(), tadj);
            hTofcvsEhit->Fill(hit->GetE(), tcorr);

            randx = (std::rand() / (float)RAND_MAX);
            const int plane = static_cast<const int>(std::floor((hit->GetPaddle() - 1) / 50)); // ig -1
            ahXYperPlane[plane]->Fill(hit->GetPosition().X() + (plane % 2) * 5. * (randx - 0.5),
                                      hit->GetPosition().Y() + ((plane + 1) % 2) * 5. * (randx - 0.5));

            hTofvsX->Fill(hit->GetPosition().X() + (plane % 2) * 5. * (randx - 0.5), tadj);
            hTofcvsX->Fill(hit->GetPosition().X() + (plane % 2) * 5. * (randx - 0.5), tcorr);
            hTofvsY->Fill(hit->GetPosition().Y() + ((plane + 1) % 2) * 5. * (randx - 0.5), tadj);
            hTofcvsY->Fill(hit->GetPosition().Y() + ((plane + 1) % 2) * 5. * (randx - 0.5), tcorr);
            hTofvsZ->Fill(plane, tadj);
            hTofcvsZ->Fill(plane, tcorr);
        }
        else
        {
            if ((fEventHeader->GetTpat() & fCosmicTpat) == fCosmicTpat)
            { // 0x2000 before,  0x100 fission 2021
                hHitEvsBarCosmics->Fill(bar, hit->GetE());
                hTdiffvsBarCosmics->Fill(bar, hit->GetTdcL() - hit->GetTdcR());

                for (const auto& hitref : hits)
                {
                    if ((hitref->GetPaddle() == 675) && (bar != 675))
                    {
                        hDT675->Fill(
                            bar, (hit->GetTdcL() + hit->GetTdcR()) / 2. - (hitref->GetTdcL() + hitref->GetTdcR()) / 2.);
                        hDT675c->Fill(bar,
                                      (hit->GetTdcL() + hit->GetTdcR()) / 2. -
                                          (hitref->GetTdcL() + hitref->GetTdcR()) / 2. +
                                          copysign(1., (hit->GetPosition() - hitref->GetPosition()).Y()) *
                                              (hit->GetPosition() - hitref->GetPosition()).Mag() / clight);
                    }
                    if ((hitref->GetPaddle() == 625) && (bar != 625))
                    {
                        hDT625->Fill(
                            bar, (hit->GetTdcL() + hit->GetTdcR()) / 2. - (hitref->GetTdcL() + hitref->GetTdcR()) / 2.);
                        hDT625c->Fill(bar,
                                      (hit->GetTdcL() + hit->GetTdcR()) / 2. -
                                          (hitref->GetTdcL() + hitref->GetTdcR()) / 2. +
                                          copysign(1., (hit->GetPosition() - hitref->GetPosition()).Y()) *
                                              (hit->GetPosition() - hitref->GetPosition()).Mag() / clight);
                    }
                }
            }
        }
    }
}

void R3BNeulandOnlineSpectra::FinishTask()
{
    TDirectory* tmp = gDirectory;
    FairRootManager::Instance()->GetOutFile()->cd();

    gDirectory->mkdir("R3BNeulandOnlineSpectra");
    gDirectory->cd("R3BNeulandOnlineSpectra");

    ahMappedBar1[0]->Write();
    ahMappedBar1[1]->Write();
    ahMappedBar1[2]->Write();
    ahMappedBar1[3]->Write();
    ahMappedBar2[0]->Write();
    ahMappedBar2[1]->Write();
    ahMappedBar2[2]->Write();
    ahMappedBar2[3]->Write();

    hTstart->Write();
    hNstart->Write();

    hTestJump->Write();
    hJumpsvsEvnt->Write();
    hJumpsvsEvntzoom->Write();

    ahCalTvsBar[0]->Write();
    ahCalTvsBar[1]->Write();
    ahCalEvsBar[0]->Write();
    ahCalEvsBar[1]->Write();

    hHitEvsBar->Write();

    hHitEvsBarCosmics->Write();

    hTdiffvsBar->Write();
    hToFvsBar->Write();
    hTofvsEhit->Write();
    hToFcvsBar->Write();
    hTofcvsEhit->Write();

    hTofvsX->Write();
    hTofcvsX->Write();
    hTofvsY->Write();
    hTofcvsY->Write();
    hTofvsZ->Write();
    hTofcvsZ->Write();

    hTdiffvsBarCosmics->Write();

    hDT675->Write();
    hDT625->Write();

    hDT675c->Write();
    hDT625c->Write();

    for (unsigned int i = 0; i < fNPlanes; i++)
    {
        ahXYperPlane[i]->Write();
    }

    hSofiaTime->Write();
    hNeuLANDvsStart->Write();
    hTOF->Write();
    hTOFc->Write();

    gDirectory = tmp;
}

void R3BNeulandOnlineSpectra::ResetHistos()
{
    ahMappedBar1[0]->Reset();
    ahMappedBar1[1]->Reset();
    ahMappedBar1[2]->Reset();
    ahMappedBar1[3]->Reset();
    ahMappedBar2[0]->Reset();
    ahMappedBar2[1]->Reset();
    ahMappedBar2[2]->Reset();
    ahMappedBar2[3]->Reset();

    hTstart->Reset();
    hNstart->Reset();

    hTestJump->Reset();
    hJumpsvsEvnt->Reset();
    hJumpsvsEvntzoom->Reset();

    ahCalTvsBar[0]->Reset();
    ahCalTvsBar[1]->Reset();
    ahCalEvsBar[0]->Reset();
    ahCalEvsBar[1]->Reset();

    hHitEvsBar->Reset();

    hHitEvsBarCosmics->Reset();

    hTdiffvsBar->Reset();
    hToFvsBar->Reset();
    hTofvsEhit->Reset();
    hToFcvsBar->Reset();
    hTofcvsEhit->Reset();

    hTofvsX->Reset();
    hTofcvsX->Reset();
    hTofvsY->Reset();
    hTofcvsY->Reset();
    hTofvsZ->Reset();
    hTofcvsZ->Reset();

    hTdiffvsBarCosmics->Reset();

    hDT675->Reset();
    hDT625->Reset();

    for (unsigned int i = 0; i < fNPlanes; i++)
    {
        ahXYperPlane[i]->Reset();
    }

    hSofiaTime->Reset();
    hNeuLANDvsStart->Reset();
    hTOF->Reset();
    hTOFc->Reset();
}

void R3BNeulandOnlineSpectra::ResetHistosMapped()
{
    ahMappedBar1[0]->Reset();
    ahMappedBar1[1]->Reset();
    ahMappedBar1[2]->Reset();
    ahMappedBar1[3]->Reset();
    ahMappedBar2[0]->Reset();
    ahMappedBar2[1]->Reset();
    ahMappedBar2[2]->Reset();
    ahMappedBar2[3]->Reset();
}

bool R3BNeulandOnlineSpectra::IsBeam() const { return !std::isnan(fEventHeader->GetTStart()); }

ClassImp(R3BNeulandOnlineSpectra)
