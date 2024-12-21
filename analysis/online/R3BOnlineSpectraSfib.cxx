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

/*
 * This task should fill histograms with detector variables which allow
 * to test the detectors online
 *
 */

#include "R3BOnlineSpectraSfib.h"
#include "R3BLosCalData.h"
#include "R3BLosMappedData.h"

#include "R3BSci8CalData.h"
#include "R3BSci8MappedData.h"

#include "R3BBeamMonitorMappedData.h"

#include "R3BTofdCalData.h"
#include "R3BTofdMappedData.h"

#include "R3BPaddleCalData.h"
#include "R3BPaddleTamexMappedData.h"

#include "R3BPspxCalData.h"
#include "R3BPspxHitData.h"
#include "R3BPspxMappedData.h"
#include "R3BPspxPrecalData.h"

#include "R3BEventHeader.h"
#include "R3BTCalEngine.h"

#include "R3BBunchedFiberCalData.h"
#include "R3BBunchedFiberHitData.h"
#include "R3BBunchedFiberMappedData.h"

#include "R3BSfibMappedData.h"

#include "FairLogger.h"
#include "FairRootManager.h"
#include "FairRunAna.h"
#include "FairRunOnline.h"
#include "FairRuntimeDb.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TGaxis.h"
#include "TH1F.h"
#include "TH2F.h"
#include "THttpServer.h"

#include "TClonesArray.h"
#include "TMath.h"
#include <TRandom3.h>
#include <TRandomGen.h>
#include <algorithm>
#include <array>
#include <cmath>
#include <cstdlib>
#include <ctime>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#define IS_NAN(x) TMath::IsNaN(x)
using namespace std;

R3BOnlineSpectraSfib::R3BOnlineSpectraSfib()
    : R3BOnlineSpectraSfib("OnlineSpectraSfib", 1)
{
}

R3BOnlineSpectraSfib::R3BOnlineSpectraSfib(const char* name, Int_t iVerbose)
    : FairTask(name, iVerbose)
    , fTrigger(-1)
    , fTpat(-1)
    , fNEvents(0)
{
}

R3BOnlineSpectraSfib::~R3BOnlineSpectraSfib()
{
    //	delete fhTpat ;
    //	delete fhTrigger;
}

InitStatus R3BOnlineSpectraSfib::Init()
{

    // Initialize random number:
    std::srand(std::time(0)); // use current time as seed for random generator

    LOG(info) << "R3BOnlineSpectraSfib::Init()";

    // try to get a handle on the EventHeader. EventHeader may not be
    // present though and hence may be null. Take care when using.

    FairRootManager* mgr = FairRootManager::Instance();
    if (NULL == mgr)
        LOG(fatal) << "FairRootManager not found";

    header = dynamic_cast<R3BEventHeader*>(mgr->GetObject("R3BEventHeader"));
    if (!header)
        header = dynamic_cast<R3BEventHeader*>(mgr->GetObject("EventHeader."));

    FairRunOnline* run = FairRunOnline::Instance();

    run->GetHttpServer()->Register("/Tasks", this);

    // Get objects for detectors on all levels
    assert(DET_MAX + 1 == sizeof(fDetectorNames) / sizeof(fDetectorNames[0]));
    for (int det = 0; det < DET_MAX; det++)
    {
        fMappedItems.push_back(dynamic_cast<TClonesArray*>(mgr->GetObject(Form("%sMapped", fDetectorNames[det]))));
        if (NULL == fMappedItems.at(det))
        {
            printf("Could not find mapped data for '%s'.\n", fDetectorNames[det]);
        }
        fCalItems.push_back(dynamic_cast<TClonesArray*>(mgr->GetObject(Form("%sCal", fDetectorNames[det]))));
        fHitItems.push_back(dynamic_cast<TClonesArray*>(mgr->GetObject(Form("%sHit", fDetectorNames[det]))));
    }

    //------------------------------------------------------------------------
    // create histograms of all detectors
    //------------------------------------------------------------------------

    if (fMappedItems.at(DET_SFIB))
    {

        TCanvas* cSfib_det = new TCanvas("Sfib_det", "Sfib det", 10, 10, 510, 510);
        cSfib_det->Divide(2, 3);

        TCanvas* cSfib_det2 = new TCanvas("Sfib_det2", "Sfib det2", 10, 10, 510, 510);
        cSfib_det2->Divide(2, 3);

        TCanvas* cSfib_det3 = new TCanvas("Sfib_det3", "Sfib det3", 10, 10, 510, 510);
        cSfib_det3->Divide(2, 3);

        TCanvas* cSfib_det4 = new TCanvas("Sfib_det4", "Sfib det4", 10, 10, 510, 510);
        cSfib_det4->Divide(2, 3);

        TCanvas* cSfib_det5 = new TCanvas("Sfib_det5", "Sfib det5", 10, 10, 510, 510);
        cSfib_det5->Divide(2, 3);

        TCanvas* cSfib_det6 = new TCanvas("Sfib_det6", "Sfib det6", 10, 10, 510, 510);
        cSfib_det6->Divide(2, 3);

        TCanvas* cSfib_det7 = new TCanvas("Sfib_det7", "Sfib det7", 10, 10, 510, 510);
        cSfib_det6->Divide(2, 2);

        fh_sfib_Tot_top = new TH2F("sfib_Tot_top", "SFIB ToT@Top", 256, 0, 256, 100, 0., 600.);
        fh_sfib_Tot_top->GetXaxis()->SetTitle("Channel");
        fh_sfib_Tot_top->GetYaxis()->SetTitle("ToT / ns");

        fh_sfib_Tot_top1d = new TH1F("sfib_Tot_top1d", "ToT@Top Dist", 500, 0., 500.);
        fh_sfib_Tot_top1d->GetXaxis()->SetTitle("ToT / ns");
        fh_sfib_Tot_top1d->GetYaxis()->SetTitle("Counts");

        fh_sfib_Tot_top_g =
            new TH2F("sfib_Tot_top_g", "SFIB ToT@Top gate with Tot@Bot>100 ", 256, 0, 256, 100, 0., 600.);
        fh_sfib_Tot_top_g->GetXaxis()->SetTitle("Channel");
        fh_sfib_Tot_top_g->GetYaxis()->SetTitle("ToT / ns");

        fh_sfib_Tot_top1d_g = new TH1F("sfib_Tot_top1d_g", "ToT@Top Dist gate with Tot@Bot>100", 500, 0., 500.);
        fh_sfib_Tot_top1d_g->GetXaxis()->SetTitle("ToT / ns");
        fh_sfib_Tot_top1d_g->GetYaxis()->SetTitle("Counts");

        fh_sfib_Tot_bot = new TH2F("sfib_Tot_bot", "SFIB ToT@Bot", 256, 0, 256, 100, 0., 600.);
        fh_sfib_Tot_bot->GetXaxis()->SetTitle("Channel");
        fh_sfib_Tot_bot->GetYaxis()->SetTitle("ToT / ns");

        fh_sfib_Tot_bot1d = new TH1F("sfib_Tot_bot1d", "ToT@Bot Dist", 500, 0., 500.);
        fh_sfib_Tot_bot1d->GetXaxis()->SetTitle("ToT / ns");
        fh_sfib_Tot_bot1d->GetYaxis()->SetTitle("Counts");

        fh_sfib_Tot_bot_g =
            new TH2F("sfib_Tot_bot_g", "SFIB ToT@Bot gate with ToT@Top>350", 256, 0, 256, 100, 0., 600.);
        fh_sfib_Tot_bot_g->GetXaxis()->SetTitle("Channel");
        fh_sfib_Tot_bot_g->GetYaxis()->SetTitle("ToT / ns");

        fh_sfib_Tot_bot1d_g = new TH1F("sfib_Tot_bot1d_g", "ToT@Bot Dist gate with Tot@Top>350", 500, 0., 500.);
        fh_sfib_Tot_bot1d_g->GetXaxis()->SetTitle("ToT / ns");
        fh_sfib_Tot_bot1d_g->GetYaxis()->SetTitle("Counts");

        fh_sfib_channels = new TH1F("sfib_channels", "SFIB Channels", 256, 0, 256);
        fh_sfib_channels->GetXaxis()->SetTitle("Channel");
        fh_sfib_channels->GetYaxis()->SetTitle("Counts");

        fh_sfib_channels_top = new TH1F("sfib_channels_top", "SFIB Channels Top", 256, 0, 256);
        fh_sfib_channels_top->GetXaxis()->SetTitle("Channel");
        fh_sfib_channels_top->GetYaxis()->SetTitle("Counts");

        fh_sfib_channels_bot = new TH1F("sfib_channels_bot", "SFIB Channels Bottom", 256, 0, 256);
        fh_sfib_channels_bot->GetXaxis()->SetTitle("Channel");
        fh_sfib_channels_bot->GetYaxis()->SetTitle("Counts");

        fh_sfib_channels_topvsbot =
            new TH2F("sfib_channels_topvsbot", "SFIB Channels Top vs Bottom", 256, 0, 256, 256, 0, 256);
        fh_sfib_channels_topvsbot->GetXaxis()->SetTitle("Top Channels");
        fh_sfib_channels_topvsbot->GetYaxis()->SetTitle("Bottom Channels");

        fh_sfib_channels_g = new TH1F("sfib_channels_g", "SFIB Channels gate", 256, 0, 256);
        fh_sfib_channels_g->GetXaxis()->SetTitle("Channel");
        fh_sfib_channels_g->GetYaxis()->SetTitle("Counts");

        // fh_sfib_tofd = new TH2F("sfib_fib_tofd","SFIB Channel vs TOFD paddle",256,0,256,180,0,180);
        // fh_sfib_tofd->GetXaxis()->SetTitle("SFIB Channel");
        // fh_sfib_tofd->GetYaxis()->SetTitle("TOFD paddle");

        // eng run Dec 2019, no tofd but psp in front of sfib
        // fh_sfib_pspx = new TH2F("sfib_fib_pspx","SFIB Channel vs PSP strip",256,0,256,64,0,64);
        // fh_sfib_pspx->GetXaxis()->SetTitle("SFIB Channel");
        // fh_sfib_pspx->GetYaxis()->SetTitle("PSP strip");

        fh_sfib_tot_tvb = new TH2F("sfib_totb_tvb", "SFIB TOT@Top vs TOT@Bot", 500, 100, 600, 500, 100, 600);
        fh_sfib_tot_tvb->GetXaxis()->SetTitle("ToT@Bot / ns");
        fh_sfib_tot_tvb->GetYaxis()->SetTitle("ToT@Top / ns");

        fh_sfib_fmult = new TH1F("sfib_fmult", "SFIB fiber multiplicity", 32, 0, 32);
        fh_sfib_fmult->GetXaxis()->SetTitle("fiber multiplicity");
        fh_sfib_fmult->GetYaxis()->SetTitle("Counts");

        fh_sfib_cmult = new TH1F("sfib_cmult", "SFIB cluster multiplicity", 32, 0, 32);
        fh_sfib_cmult->GetXaxis()->SetTitle("cluster multiplicity");
        fh_sfib_cmult->GetYaxis()->SetTitle("Counts");

        fh_sfib_fmult_g = new TH1F("sfib_fmult_g", "SFIB fiber multiplicity with gate", 32, 0, 32);
        fh_sfib_fmult_g->GetXaxis()->SetTitle("fiber multiplicity");
        fh_sfib_fmult_g->GetYaxis()->SetTitle("Counts");

        fh_sfib_cmult_g = new TH1F("sfib_cmult_g", "SFIB cluster multiplicity with gate", 32, 0, 32);
        fh_sfib_cmult_g->GetXaxis()->SetTitle("cluster multiplicity");
        fh_sfib_cmult_g->GetYaxis()->SetTitle("Counts");

        fh_sfib_tot_ch_top1 = new TH1F("sfib_tot_ch_top1", "SFIB TOT Top1 Channels", 256, 1, 256);
        fh_sfib_tot_ch_top1->GetXaxis()->SetTitle("ToT@Top1 / ns");
        fh_sfib_tot_ch_top1->GetYaxis()->SetTitle("Counts");
        fh_sfib_tot_ch_top1->SetStats(0);

        fh_sfib_tot_ch_top2 = new TH1F("sfib_tot_ch_top2", "SFIB TOT Top2 Channels", 256, 1, 256);
        fh_sfib_tot_ch_top2->GetXaxis()->SetTitle("ToT@Top2 / ns");
        fh_sfib_tot_ch_top2->GetYaxis()->SetTitle("Counts");
        fh_sfib_tot_ch_top2->SetStats(0);

        fh_sfib_tot_ch_bot1 = new TH1F("sfib_tot_ch_bot1", "SFIB TOT Bottom1 Channels", 256, 1, 256);
        fh_sfib_tot_ch_bot1->GetXaxis()->SetTitle("ToT@Bottom1 / ns");
        fh_sfib_tot_ch_bot1->GetYaxis()->SetTitle("Counts");
        fh_sfib_tot_ch_bot1->SetStats(0);

        fh_sfib_tot_ch_bot2 = new TH1F("sfib_tot_ch_bot2", "SFIB TOT Bottom2 Channels", 256, 1, 256);
        fh_sfib_tot_ch_bot2->GetXaxis()->SetTitle("ToT@Bottom2 / ns");
        fh_sfib_tot_ch_bot2->GetYaxis()->SetTitle("Counts");
        fh_sfib_tot_ch_bot2->SetStats(0);

        fh_sfib_tot_max_top1 = new TH1F("sfib_tot_max_top1", "SFIB TOT Top1 Maximum", 119, 5, 600);
        fh_sfib_tot_max_top1->GetXaxis()->SetTitle("Maximum");
        fh_sfib_tot_max_top1->GetYaxis()->SetTitle("Counts");
        fh_sfib_tot_max_top1->SetStats(0);

        fh_sfib_tot_max_top2 = new TH1F("sfib_tot_max_top2", "SFIB TOT Top2 Maximum", 119, 5, 600);
        fh_sfib_tot_max_top2->GetXaxis()->SetTitle("Maximum");
        fh_sfib_tot_max_top2->GetYaxis()->SetTitle("Counts");
        fh_sfib_tot_max_top2->SetStats(0);

        fh_sfib_tot_max_bot1 = new TH1F("sfib_tot_max_bot1", "SFIB TOT Bottom1 Maximum", 119, 5, 600);
        fh_sfib_tot_max_bot1->GetXaxis()->SetTitle("Maximum");
        fh_sfib_tot_max_bot1->GetYaxis()->SetTitle("Counts");
        fh_sfib_tot_max_bot1->SetStats(0);

        fh_sfib_tot_max_bot2 = new TH1F("sfib_tot_max_bot2", "SFIB TOT Bottom2 Maximum", 119, 5, 600);
        fh_sfib_tot_max_bot2->GetXaxis()->SetTitle("Maximum");
        fh_sfib_tot_max_bot2->GetYaxis()->SetTitle("Counts");
        fh_sfib_tot_max_bot2->SetStats(0);

        fh_sfib_tot_ch_max_top1 =
            new TH2F("sfib_tot_ch_max_top1", "SFIB TOT Top1 Channel vs Maximum", 256, 1, 256, 500, 100, 600);
        fh_sfib_tot_ch_max_top1->GetXaxis()->SetTitle("ToT@Top1 / ns");
        fh_sfib_tot_ch_max_top1->GetYaxis()->SetTitle("Maximum");
        fh_sfib_tot_ch_max_top1->SetStats(0);

        fh_sfib_tot_ch_max_top2 =
            new TH2F("sfib_tot_ch_max_top2", "SFIB TOT Top2 Channel vs Maximum", 256, 1, 256, 500, 100, 600);
        fh_sfib_tot_ch_max_top2->GetXaxis()->SetTitle("ToT@Top2 / ns");
        fh_sfib_tot_ch_max_top2->GetYaxis()->SetTitle("Maximum");
        fh_sfib_tot_ch_max_top2->SetStats(0);

        fh_sfib_tot_ch_max_bot1 =
            new TH2F("sfib_tot_ch_max_bot1", "SFIB TOT Bottom1 Channel vs Maximum", 256, 1, 256, 500, 100, 600);
        fh_sfib_tot_ch_max_bot1->GetXaxis()->SetTitle("ToT@Bottom1 / ns");
        fh_sfib_tot_ch_max_bot1->GetYaxis()->SetTitle("Maximum");
        fh_sfib_tot_ch_max_bot1->SetStats(0);

        fh_sfib_tot_ch_max_bot2 =
            new TH2F("sfib_tot_ch_max_bot2", "SFIB TOT Bottom2 Channel vs Maximum", 256, 1, 256, 500, 100, 600);
        fh_sfib_tot_ch_max_bot2->GetXaxis()->SetTitle("ToT@Bottom2 / ns");
        fh_sfib_tot_ch_max_bot2->GetYaxis()->SetTitle("Maximum");
        fh_sfib_tot_ch_max_bot2->SetStats(0);

        fh_sfib_tot_ch_topbot1 =
            new TH2F("sfib_tot_ch_topbot1", "SFIB TOT TopvsBottom1 Channel", 256, 1, 256, 256, 1, 256);
        fh_sfib_tot_ch_topbot1->GetXaxis()->SetTitle("ToT@Top1 / ns");
        fh_sfib_tot_ch_topbot1->GetYaxis()->SetTitle("ToT!Bottom1 / ns");
        fh_sfib_tot_ch_topbot1->SetStats(0);

        fh_sfib_multi = new TH1F("sfib_multi", "SFIB Multiplicity", 5, 1, 6);
        fh_sfib_multi->GetXaxis()->SetTitle("Multiplicity");
        fh_sfib_multi->GetYaxis()->SetTitle("Counts");
        fh_sfib_multi->SetStats(0);

        fh_sfib_clus = new TH1F("sfib_clus", "SFIB Clustering", 5, 1, 6);
        fh_sfib_clus->GetXaxis()->SetTitle("Clustering");
        fh_sfib_clus->GetYaxis()->SetTitle("Counts");
        fh_sfib_clus->SetStats(0);

        fh_sfib_multi_clus = new TH2F("sfib_multi_clus", "SFIB Multiplicity vs Clustering", 4, 1, 5, 4, 1, 5);
        fh_sfib_multi_clus->GetXaxis()->SetTitle("Multiplicity");
        fh_sfib_multi_clus->GetYaxis()->SetTitle("Clustering");
        fh_sfib_multi_clus->SetStats(0);

        fh_sfib_multi_ch_topbot1 =
            new TH2F("sfib_tot_multi_topbot1", "SFIB Multiplicity TopvsBottom1 Channel", 40, 161, 200, 40, 161, 200);
        fh_sfib_multi_ch_topbot1->GetXaxis()->SetTitle("ToT@Top1 / ns");
        fh_sfib_multi_ch_topbot1->GetYaxis()->SetTitle("ToT@Bottom1 / ns");
        fh_sfib_multi_ch_topbot1->SetStats(0);

        fh_sfib_tot_top_multi1 = new TH1F("sfib_tot_top_multi1", "SFIB TOT Top Multiplicity1", 119, 5, 600);
        fh_sfib_tot_top_multi1->GetXaxis()->SetTitle("ToT@Top1");
        fh_sfib_tot_top_multi1->GetYaxis()->SetTitle("Counts");
        fh_sfib_tot_top_multi1->SetStats(0);

        fh_sfib_tot_top_multi2 = new TH1F("sfib_tot_top_multi2", "SFIB TOT Top Multiplicity2", 119, 5, 600);
        fh_sfib_tot_top_multi2->GetXaxis()->SetTitle("ToT@Top2");
        fh_sfib_tot_top_multi2->GetYaxis()->SetTitle("Counts");
        fh_sfib_tot_top_multi2->SetStats(0);

        cSfib_det->cd(1);
        fh_sfib_channels_top->Draw();

        cSfib_det->cd(2);
        fh_sfib_channels_bot->Draw();

        cSfib_det->cd(3);
        xyline = new TF1("xyline", "x", 0, 256);
        fh_sfib_channels_topvsbot->Draw("colz");
        xyline->Draw("same");

        cSfib_det->cd(4);
        fh_sfib_Tot_top1d->Draw();

        cSfib_det->cd(5);
        fh_sfib_Tot_bot1d->Draw();

        cSfib_det->cd(6);
        fh_sfib_Tot_top1d_g->Draw();

        cSfib_det2->cd(1);
        fh_sfib_Tot_top->Draw("colz");

        cSfib_det2->cd(2);
        fh_sfib_Tot_top1d->Draw();

        cSfib_det2->cd(3);
        fh_sfib_Tot_top_g->Draw("colz");

        cSfib_det2->cd(4);
        fh_sfib_Tot_top1d_g->Draw();

        cSfib_det2->cd(5);
        fh_sfib_Tot_bot->Draw("colz");

        cSfib_det2->cd(6);
        fh_sfib_Tot_bot1d->Draw();

        cSfib_det3->cd(1);
        fh_sfib_Tot_bot_g->Draw("colz");

        cSfib_det3->cd(2);
        fh_sfib_Tot_bot1d_g->Draw();

        cSfib_det3->cd(3);
        fh_sfib_channels->Draw();

        cSfib_det3->cd(4);
        fh_sfib_channels_g->Draw();

        // cSfib_det->cd(11);
        // fh_sfib_tofd->Draw("colz");
        // fh_sfib_pspx->Draw("colz");

        cSfib_det3->cd(5);
        fh_sfib_tot_tvb->Draw("colz");

        cSfib_det3->cd(6);
        gPad->SetLogy();
        fh_sfib_fmult->Draw();

        cSfib_det4->cd(1);
        gPad->SetLogy();
        fh_sfib_cmult->Draw();

        cSfib_det4->cd(2);
        gPad->SetLogy();
        fh_sfib_fmult_g->Draw();

        cSfib_det4->cd(3);
        gPad->SetLogy();
        fh_sfib_cmult_g->Draw();

        cSfib_det4->cd(4);
        fh_sfib_tot_ch_top1->Draw();

        cSfib_det4->cd(5);
        fh_sfib_tot_ch_top2->Draw();

        cSfib_det4->cd(6);
        fh_sfib_tot_ch_bot1->Draw();

        cSfib_det5->cd(1);
        fh_sfib_tot_ch_top1->Draw();

        cSfib_det5->cd(2);
        fh_sfib_tot_max_top1->Draw();

        cSfib_det5->cd(3);
        fh_sfib_tot_max_top2->Draw();

        cSfib_det5->cd(4);
        fh_sfib_tot_max_bot1->Draw();

        cSfib_det5->cd(5);
        fh_sfib_tot_max_bot2->Draw();

        cSfib_det5->cd(6);
        fh_sfib_tot_ch_max_top1->Draw();

        cSfib_det6->cd(1);
        fh_sfib_tot_ch_max_top2->Draw();

        cSfib_det6->cd(2);
        fh_sfib_tot_ch_max_bot1->Draw();

        cSfib_det6->cd(3);
        fh_sfib_tot_ch_max_bot2->Draw();

        cSfib_det6->cd(4);
        fh_sfib_tot_ch_topbot1->Draw("colz");

        cSfib_det6->cd(5);
        fh_sfib_multi->Draw();

        cSfib_det6->cd(6);
        fh_sfib_clus->Draw();

        cSfib_det7->cd(1);
        fh_sfib_multi_clus->Draw();

        cSfib_det7->cd(2);
        fh_sfib_multi_ch_topbot1->Draw("colz");

        cSfib_det7->cd(3);
        fh_sfib_tot_top_multi1->Draw();

        cSfib_det7->cd(4);
        fh_sfib_tot_top_multi2->Draw();

        run->AddObject(cSfib_det);
        run->GetHttpServer()->RegisterCommand("Reset_SFIB", Form("/Tasks/%s/->Reset_SFIB_Histo()", GetName()));
    }

    // -------------------------------------------------------------------------

    return kSUCCESS;
}

void R3BOnlineSpectraSfib::Reset_SFIB_Histo()
{
    fh_sfib_Tot_top->Reset();
    fh_sfib_Tot_top1d->Reset();
    fh_sfib_Tot_top_g->Reset();
    fh_sfib_Tot_top1d_g->Reset();
    fh_sfib_Tot_bot->Reset();
    fh_sfib_Tot_bot1d->Reset();
    fh_sfib_Tot_bot_g->Reset();
    fh_sfib_Tot_bot1d_g->Reset();
    fh_sfib_channels->Reset();
    fh_sfib_channels_top->Reset();
    fh_sfib_channels_bot->Reset();
    fh_sfib_channels_topvsbot->Reset();
    fh_sfib_tot_tvb->Reset();
    // fh_sfib_tofd->Write();
    // fh_sfib_pspx->Write();
    fh_sfib_fmult->Reset();
    fh_sfib_cmult->Reset();
    fh_sfib_fmult_g->Reset();
    fh_sfib_cmult_g->Reset();
    fh_sfib_tot_ch_top1->Reset();
    fh_sfib_tot_ch_top2->Reset();
    fh_sfib_tot_ch_bot1->Reset();
    fh_sfib_tot_ch_bot2->Reset();
    fh_sfib_tot_max_top1->Reset();
    fh_sfib_tot_max_top2->Reset();
    fh_sfib_tot_max_bot1->Reset();
    fh_sfib_tot_max_bot2->Reset();
    fh_sfib_tot_ch_max_top1->Reset();
    fh_sfib_tot_ch_max_top2->Reset();
    fh_sfib_tot_ch_max_bot1->Reset();
    fh_sfib_tot_ch_max_bot2->Reset();
    fh_sfib_tot_ch_topbot1->Reset();
    fh_sfib_multi->Reset();
    fh_sfib_clus->Reset();
    fh_sfib_multi_clus->Reset();
    fh_sfib_multi_ch_topbot1->Reset();
    fh_sfib_tot_top_multi1->Reset();
    fh_sfib_tot_top_multi2->Reset();
    fh_sfib_channels_top->Reset();
    fh_sfib_channels_bot->Reset();
    fh_sfib_channels_topvsbot->Reset();
}

void R3BOnlineSpectraSfib::Exec(Option_t* option)
{
    //  cout << "fNEvents " << fNEvents << endl;

    FairRootManager* mgr = FairRootManager::Instance();
    if (NULL == mgr)
    {
        // FairLogger::GetLogger()->Fatal(MESSAGE_ORIGIN, "FairRootManager not found");
        LOG(error) << "FairRootManager not found";
        return;
    }

    //----------------------------------------------------------------------
    // Sipm Fiber
    //----------------------------------------------------------------------

    // Double_t sfib_tot[2]={-100};//sfib_tot[i],i=0(top),1(bot)
    // Double_t sfib_tot_max[2][5]={-100};//sfib_tot_ch[i][j], i=0(top),1(bot),j=1st,2nd,3rd.. largest TOT value
    // Int_t sfib_tot_ch[2][5]={-100}; //sfib_tot_ch[i][j], i=0(top),1(bot),j=1st,2nd,3rd.. largest channel
    if (fMappedItems.at(DET_SFIB))
    {
        uint32_t lt[2][256][4];
        memset(lt, 0, sizeof lt);
        char mask[256];
        memset(mask, 0, sizeof mask);

        auto det = fMappedItems.at(DET_SFIB);
        Int_t nMapped = det->GetEntriesFast();
        for (Int_t imapped = 0; imapped < nMapped; imapped++)
        {
            auto mapped = dynamic_cast<R3BSfibMappedData const*>(det->At(imapped));
            if (!mapped)
            {
                continue; // should not happen
            }
            auto& lt_this = lt[mapped->IsTop() ? 0 : 1][mapped->GetChannel() - 1];
            lt_this[mapped->IsLeading() ? 0 : 1] = mapped->GetCoarse();
            lt_this[mapped->IsLeading() ? 2 : 3] = mapped->GetFine();

            mask[mapped->GetChannel() - 1] |= mapped->IsTop() ? 1 : 2;
        }
        Int_t fib_mult = 0;
        Int_t fib_mult_g = 0;
        Int_t fib_fired[256] = { 0 };
        Int_t fib_fired_g[256] = { 0 };

        Double_t sfib_tot[2];        // sfib_tot[i],i=0(top),1(bot)
        Double_t sfib_tot_max[2][5]; // sfib_tot_ch[i][j], i=0(top),1(bot),j=1st,2nd,3rd.. largest TOT value
        Int_t sfib_tot_ch[2][5];     // sfib_tot_ch[i][j], i=0(top),1(bot),j=1st,2nd,3rd.. largest channel

        for (int i = 0; i < 2; i++)
        {
            sfib_tot[2] = -100;
            for (int j = 0; j < 5; j++)
            {
                sfib_tot_max[i][j] = -100.; // sfib_tot_ch[i][j], i=0(top),1(bot),j=1st,2nd,3rd.. largest TOT value
                sfib_tot_ch[i][j] = -100;
            }
        }

        // std::ofstream run202("run202.txt",std::ios::app);
        // run202 << fNEvents;
        for (int ch = 0; ch < 256; ++ch)
        {
            sfib_tot[0] = (16 * lt[0][ch][1] + lt[0][ch][3]) - (16 * lt[0][ch][0] + lt[0][ch][2]); //[0] for top
            sfib_tot[1] = (16 * lt[1][ch][1] + lt[1][ch][3]) - (16 * lt[1][ch][0] + lt[1][ch][2]); //[1] for bot
            //
            // Top
            //
            if ((sfib_tot[0] > sfib_tot_max[0][1]) && (sfib_tot[0] < 4000))
            {
                // sfib_tot_max[0][6]=sfib_tot_max[0][5];
                sfib_tot_max[0][5] = sfib_tot_max[0][4];
                sfib_tot_max[0][4] = sfib_tot_max[0][3];
                sfib_tot_max[0][3] = sfib_tot_max[0][2];
                sfib_tot_max[0][2] = sfib_tot_max[0][1];
                sfib_tot_max[0][1] = sfib_tot[0];
                // sfib_tot_ch[0][6]=sfib_tot_ch[0][5];
                sfib_tot_ch[0][5] = sfib_tot_ch[0][4];
                sfib_tot_ch[0][4] = sfib_tot_ch[0][3];
                sfib_tot_ch[0][3] = sfib_tot_ch[0][2];
                sfib_tot_ch[0][2] = sfib_tot_ch[0][1];
                sfib_tot_ch[0][1] = ch;
            }
            else if ((sfib_tot[0] > sfib_tot_max[0][2]) && (sfib_tot[0] < 4000))
            {
                // sfib_tot_max[0][6]=sfib_tot_max[0][5];
                sfib_tot_max[0][5] = sfib_tot_max[0][4];
                sfib_tot_max[0][4] = sfib_tot_max[0][3];
                sfib_tot_max[0][3] = sfib_tot_max[0][2];
                sfib_tot_max[0][2] = sfib_tot[0];
                sfib_tot_max[0][1] = sfib_tot_max[0][1];
                // sfib_tot_ch[0][6]=sfib_tot_ch[0][5];
                sfib_tot_ch[0][5] = sfib_tot_ch[0][4];
                sfib_tot_ch[0][4] = sfib_tot_ch[0][3];
                sfib_tot_ch[0][3] = sfib_tot_ch[0][2];
                sfib_tot_ch[0][2] = ch;
                sfib_tot_ch[0][1] = sfib_tot_ch[0][1];
            }
            else if ((sfib_tot[0] > sfib_tot_max[0][3]) && (sfib_tot[0] < 4000))
            {
                // sfib_tot_max[0][6]=sfib_tot_max[0][5];
                sfib_tot_max[0][5] = sfib_tot_max[0][4];
                sfib_tot_max[0][4] = sfib_tot_max[0][3];
                sfib_tot_max[0][3] = sfib_tot[0];
                sfib_tot_max[0][2] = sfib_tot_max[0][2];
                sfib_tot_max[0][1] = sfib_tot_max[0][1];
                // sfib_tot_ch[0][6]=sfib_tot_ch[0][5];
                sfib_tot_ch[0][5] = sfib_tot_ch[0][4];
                sfib_tot_ch[0][4] = sfib_tot_ch[0][3];
                sfib_tot_ch[0][3] = ch;
                sfib_tot_ch[0][2] = sfib_tot_ch[0][2];
                sfib_tot_ch[0][1] = sfib_tot_ch[0][1];
            }
            else if ((sfib_tot[0] > sfib_tot_max[0][4]) && (sfib_tot[0] < 4000))
            {
                // sfib_tot_max[0][6]=sfib_tot_max[0][5];
                sfib_tot_max[0][5] = sfib_tot_max[0][4];
                sfib_tot_max[0][4] = sfib_tot[0];
                sfib_tot_max[0][3] = sfib_tot_max[0][3];
                sfib_tot_max[0][2] = sfib_tot_max[0][2];
                sfib_tot_max[0][1] = sfib_tot_max[0][1];
                // sfib_tot_ch[0][6]=sfib_tot_ch[0][5];
                sfib_tot_ch[0][5] = sfib_tot_ch[0][4];
                sfib_tot_ch[0][4] = ch;
                sfib_tot_ch[0][3] = sfib_tot_ch[0][3];
                sfib_tot_ch[0][2] = sfib_tot_ch[0][2];
                sfib_tot_ch[0][1] = sfib_tot_ch[0][1];
            }
            else if ((sfib_tot[0] > sfib_tot_max[0][5]) && (sfib_tot[0] < 4000))
            {
                // sfib_tot_max[0][6]=sfib_tot_max[0][5];
                sfib_tot_max[0][5] = sfib_tot[0];
                sfib_tot_max[0][4] = sfib_tot_max[0][4];
                sfib_tot_max[0][3] = sfib_tot_max[0][3];
                sfib_tot_max[0][2] = sfib_tot_max[0][2];
                sfib_tot_max[0][1] = sfib_tot_max[0][1];
                // sfib_tot_ch[0][6]=sfib_tot_ch[0][5];
                sfib_tot_ch[0][5] = ch;
                sfib_tot_ch[0][4] = sfib_tot_ch[0][4];
                sfib_tot_ch[0][3] = sfib_tot_ch[0][3];
                sfib_tot_ch[0][2] = sfib_tot_ch[0][2];
                sfib_tot_ch[0][1] = sfib_tot_ch[0][1];
            }
            // else if((sfib_tot[0]>sfib_tot_max[0][6])&&(sfib_tot[0]<4000))
            //{
            //  sfib_tot_max[0][6]=sfib_tot[0];
            //  sfib_tot_max[0][5]=sfib_tot_max[0][5];
            //  sfib_tot_max[0][4]=sfib_tot_max[0][4];
            //  sfib_tot_max[0][3]=sfib_tot_max[0][3];
            //  sfib_tot_max[0][2]=sfib_tot_max[0][2];
            //  sfib_tot_max[0][1]=sfib_tot_max[0][1];
            //  sfib_tot_ch[0][6]=ch;
            //  sfib_tot_ch[0][5]=sfib_tot_ch[0][5];
            //  sfib_tot_ch[0][4]=sfib_tot_ch[0][4];
            //  sfib_tot_ch[0][3]=sfib_tot_ch[0][3];
            //  sfib_tot_ch[0][2]=sfib_tot_ch[0][2];
            //  sfib_tot_ch[0][1]=sfib_tot_ch[0][1];
            //}
            else
            {
                // sfib_tot_max[0][6]=sfib_tot_max[0][6];
                sfib_tot_max[0][5] = sfib_tot_max[0][5];
                sfib_tot_max[0][4] = sfib_tot_max[0][4];
                sfib_tot_max[0][3] = sfib_tot_max[0][3];
                sfib_tot_max[0][2] = sfib_tot_max[0][2];
                sfib_tot_max[0][1] = sfib_tot_max[0][1];
                // sfib_tot_ch[0][6]=sfib_tot_ch[0][6];
                sfib_tot_ch[0][5] = sfib_tot_ch[0][5];
                sfib_tot_ch[0][4] = sfib_tot_ch[0][4];
                sfib_tot_ch[0][3] = sfib_tot_ch[0][3];
                sfib_tot_ch[0][2] = sfib_tot_ch[0][2];
                sfib_tot_ch[0][1] = sfib_tot_ch[0][1];
            }
            //
            // Bottom
            //
            if ((sfib_tot[1] > sfib_tot_max[1][1]) && (sfib_tot[1] < 4000))
            {
                // sfib_tot_max[1][6]=sfib_tot_max[1][5];
                sfib_tot_max[1][5] = sfib_tot_max[1][4];
                sfib_tot_max[1][4] = sfib_tot_max[1][3];
                sfib_tot_max[1][3] = sfib_tot_max[1][2];
                sfib_tot_max[1][2] = sfib_tot_max[1][1];
                sfib_tot_max[1][1] = sfib_tot[1];
                // sfib_tot_ch[1][6]=sfib_tot_ch[1][5];
                sfib_tot_ch[1][5] = sfib_tot_ch[1][4];
                sfib_tot_ch[1][4] = sfib_tot_ch[1][3];
                sfib_tot_ch[1][3] = sfib_tot_ch[1][2];
                sfib_tot_ch[1][2] = sfib_tot_ch[1][1];
                sfib_tot_ch[1][1] = ch;
            }
            else if ((sfib_tot[1] > sfib_tot_max[1][2]) && (sfib_tot[1] < 4000))
            {
                // sfib_tot_max[1][6]=sfib_tot_max[1][5];
                sfib_tot_max[1][5] = sfib_tot_max[1][4];
                sfib_tot_max[1][4] = sfib_tot_max[1][3];
                sfib_tot_max[1][3] = sfib_tot_max[1][2];
                sfib_tot_max[1][2] = sfib_tot[1];
                sfib_tot_max[1][1] = sfib_tot_max[1][1];
                // sfib_tot_ch[1][6]=sfib_tot_ch[1][5];
                sfib_tot_ch[1][5] = sfib_tot_ch[1][4];
                sfib_tot_ch[1][4] = sfib_tot_ch[1][3];
                sfib_tot_ch[1][3] = sfib_tot_ch[1][2];
                sfib_tot_ch[1][2] = ch;
                sfib_tot_ch[1][1] = sfib_tot_ch[1][1];
            }
            else if ((sfib_tot[1] > sfib_tot_max[1][3]) && (sfib_tot[1] < 4000))
            {
                // sfib_tot_max[1][6]=sfib_tot_max[1][5];
                sfib_tot_max[1][5] = sfib_tot_max[1][4];
                sfib_tot_max[1][4] = sfib_tot_max[1][3];
                sfib_tot_max[1][3] = sfib_tot[1];
                sfib_tot_max[1][2] = sfib_tot_max[1][2];
                sfib_tot_max[1][1] = sfib_tot_max[1][1];
                // sfib_tot_ch[1][6]=sfib_tot_ch[1][5];
                sfib_tot_ch[1][5] = sfib_tot_ch[1][4];
                sfib_tot_ch[1][4] = sfib_tot_ch[1][3];
                sfib_tot_ch[1][3] = ch;
                sfib_tot_ch[1][2] = sfib_tot_ch[1][2];
                sfib_tot_ch[1][1] = sfib_tot_ch[1][1];
            }
            else if ((sfib_tot[1] > sfib_tot_max[1][4]) && (sfib_tot[1] < 4000))
            {
                // sfib_tot_max[1][6]=sfib_tot_max[1][5];
                sfib_tot_max[1][5] = sfib_tot_max[1][4];
                sfib_tot_max[1][4] = sfib_tot[1];
                sfib_tot_max[1][3] = sfib_tot_max[1][3];
                sfib_tot_max[1][2] = sfib_tot_max[1][2];
                sfib_tot_max[1][1] = sfib_tot_max[1][1];
                // sfib_tot_ch[1][6]=sfib_tot_ch[1][5];
                sfib_tot_ch[1][5] = sfib_tot_ch[1][4];
                sfib_tot_ch[1][4] = ch;
                sfib_tot_ch[1][3] = sfib_tot_ch[1][3];
                sfib_tot_ch[1][2] = sfib_tot_ch[1][2];
                sfib_tot_ch[1][1] = sfib_tot_ch[1][1];
            }
            else if ((sfib_tot[1] > sfib_tot_max[1][5]) && (sfib_tot[1] < 4000))
            {
                // sfib_tot_max[1][6]=sfib_tot_max[1][5];
                sfib_tot_max[1][5] = sfib_tot[1];
                sfib_tot_max[1][4] = sfib_tot_max[1][4];
                sfib_tot_max[1][3] = sfib_tot_max[1][3];
                sfib_tot_max[1][2] = sfib_tot_max[1][2];
                sfib_tot_max[1][1] = sfib_tot_max[1][1];
                // sfib_tot_ch[1][6]=sfib_tot_ch[1][5];
                sfib_tot_ch[1][5] = ch;
                sfib_tot_ch[1][4] = sfib_tot_ch[1][4];
                sfib_tot_ch[1][3] = sfib_tot_ch[1][3];
                sfib_tot_ch[1][2] = sfib_tot_ch[1][2];
                sfib_tot_ch[1][1] = sfib_tot_ch[1][1];
            }
            // else if((sfib_tot[1]>sfib_tot_max[1][6])&&(sfib_tot[1]<4000))
            //{
            //  sfib_tot_max[1][6]=sfib_tot[1];
            //  sfib_tot_max[1][5]=sfib_tot_max[1][5];
            //  sfib_tot_max[1][4]=sfib_tot_max[1][4];
            //  sfib_tot_max[1][3]=sfib_tot_max[1][3];
            //  sfib_tot_max[1][2]=sfib_tot_max[1][2];
            //  sfib_tot_max[1][1]=sfib_tot_max[1][1];
            //  sfib_tot_ch[1][6]=ch;
            //  sfib_tot_ch[1][5]=sfib_tot_ch[1][5];
            //  sfib_tot_ch[1][4]=sfib_tot_ch[1][4];
            //  sfib_tot_ch[1][3]=sfib_tot_ch[1][3];
            //  sfib_tot_ch[1][2]=sfib_tot_ch[1][2];
            //  sfib_tot_ch[1][1]=sfib_tot_ch[1][1];
            //}
            else
            {
                // sfib_tot_max[1][6]=sfib_tot_max[1][6];
                sfib_tot_max[1][5] = sfib_tot_max[1][5];
                sfib_tot_max[1][4] = sfib_tot_max[1][4];
                sfib_tot_max[1][3] = sfib_tot_max[1][3];
                sfib_tot_max[1][2] = sfib_tot_max[1][2];
                sfib_tot_max[1][1] = sfib_tot_max[1][1];
                // sfib_tot_ch[1][6]=sfib_tot_ch[1][6];
                sfib_tot_ch[1][5] = sfib_tot_ch[1][5];
                sfib_tot_ch[1][4] = sfib_tot_ch[1][4];
                sfib_tot_ch[1][3] = sfib_tot_ch[1][3];
                sfib_tot_ch[1][2] = sfib_tot_ch[1][2];
                sfib_tot_ch[1][1] = sfib_tot_ch[1][1];
            }
            /*if (ch>-1&&ch<256&&(sfib_tot[0]>0||sfib_tot[1]>0))
              {
                cout << fNEvents;
                cout << " ";
                cout << " TOP= ";
                cout << ch;
                cout << " TOT= ";
                cout << sfib_tot[0]; //Top
                cout << " ";
                cout << (12 * lt[0][ch][1] + lt[0][ch][3]); //Leading
                cout << " ";
                cout << (12 * lt[0][ch][0] + lt[0][ch][2]); //Trailing
                cout << "\n";
                }*/
            // run202 << "Writing this to a file.\n"<<endl;
            if (ch > -1 && ch < 256 &&
                ((sfib_tot[0] > 0 && sfib_tot[0] < 10000) || (sfib_tot[1] > 0 && sfib_tot[1] < 10000)))
            // if (ch>-1&&ch<256&&sfib_tot[0]>0);
            // ASCII
            {
                //	string someString("This is a string.");
                //	printf("%s\n", someString);
                //
                // run202 << " TOP ";
                // run202 << ch;
                // run202 << " ";
                // run202 << sfib_tot[0]; //Top
                /* run202 << " ";
                run202 << (12 * lt[0][ch][1] + lt[0][ch][3]); //Leading
                run202 << " ";
          6yes0	  run202 << (12 * lt[0][ch][0] + lt[0][ch][2]); //Trailing0
                run202 << " ";*/
                // run202.close();
                //	  return 0;
                //}
                //	if (ch>-1&&ch<256&&sfib_tot[1]>0);
                // ASCII
                //{
                // std::ofstream run202("run202.txt",std::ios::app);
                // run202 << " Bottom ";
                // run202 << ch+1000;
                // run202 << " ";
                // run202 << sfib_tot[1]; //Bottom
                /*run202 << " ";
                run202 << (12 * lt[1][ch][1] + lt[1][ch][3]); //Leading
                run202 << " ";
                run202 << (12 * lt[1][ch][0] + lt[1][ch][2]); //Trailing*/
                // run202 << " ";
                // run202 << /n<< std::endl;
                // return 0;
                // run202.close();
            }
            if (lt[0][ch][0] && lt[1][ch][0])
            {
                fh_sfib_tot_tvb->Fill(sfib_tot[1], sfib_tot[0]);
            }
            if (lt[0][ch][0])
            {
                fh_sfib_channels_top->Fill(ch);
                fh_sfib_Tot_top->Fill(ch, sfib_tot[0]);
                fh_sfib_Tot_top1d->Fill(sfib_tot[0]);
                //	    cout << " ";
                //	    cout << fh_sfib_Tot_top1d->Fill(sfib_tot[0]); //Top
                for (int ch1 = 0; ch1 < 256; ch1++)
                {
                    if (lt[1][ch1][0])
                    {
                        fh_sfib_channels_topvsbot->Fill(ch, ch1);
                        /*if (ch/8 == ch1/8){
                          fh_sfib_pspx->Fill(ch,pspx0_m1_ch);
                          }*/
                    }
                }
                if (sfib_tot[1] > 100)
                {
                    fh_sfib_Tot_top_g->Fill(ch, sfib_tot[0]);
                    fh_sfib_Tot_top1d_g->Fill(sfib_tot[0]);
                }
            }
            if (lt[1][ch][0])
            {
                fh_sfib_channels_bot->Fill(ch);
                fh_sfib_Tot_bot->Fill(ch, sfib_tot[1]);
                fh_sfib_Tot_bot1d->Fill(sfib_tot[1]);
                if (sfib_tot[0] > 100)
                {
                    fh_sfib_Tot_bot_g->Fill(ch, sfib_tot[1]);
                    fh_sfib_Tot_bot1d_g->Fill(sfib_tot[1]);
                }
            }
            // multiplicity analysis for any connection of SiPM
            if (3 == mask[ch])
            {
                fh_sfib_channels->Fill(ch);
                fib_mult++;
                fib_fired[ch] = 1;
                if ((sfib_tot[1] > 100) && (sfib_tot[0] > 100))
                {
                    fib_mult_g++;
                    fib_fired_g[ch] = 1;
                    fh_sfib_channels_g->Fill(ch);
                    /*for (int ibar=0;ibar<176;ibar++)
                      {
                    if (XPad[ibar]==1) fh_sfib_tofd->Fill(ch,ibar);
                    }*/
                }
            }
            if (fib_mult > 0)
            {
                fh_sfib_fmult->Fill(fib_mult);
                if (fib_mult_g > 0)
                    fh_sfib_fmult_g->Fill(fib_mult_g);
                Int_t fib_cluster_mult = fib_mult;
                Int_t fib_cluster_mult_g = fib_mult_g;
                for (int i = 0; i < 256; i++)
                {
                    if (fib_fired[i] != 1)
                        continue;
                    if (i < 255 && fib_fired[i + 1] == 1)
                    {
                        fib_cluster_mult--;
                        if (i < 255 && fib_fired_g[i] == 1 && fib_fired_g[i + 1] == 1)
                            fib_cluster_mult_g--;
                    }
                }
                fh_sfib_cmult->Fill(fib_cluster_mult);
                if (fib_mult_g > 0)
                    fh_sfib_cmult_g->Fill(fib_cluster_mult_g);
            }
        }
        // run202<<std::endl;
        // run202.close();
        // 2019.03.12 multiplicity upto 3 J.Tanaka
        Int_t sfib_multi = 0;
        // mapping of fired fiber 2019.03.12 J.Tanaka
        Int_t sfib_fired[5] = { -100 }; // sfib_fired[multiplicity]
        Int_t sfib_map[256][256] = {
            { -100 }
        }; // convert from Sipm number on top 256 times bot 256 to fiber number 4096. sfib_map[top][bot]
        Double_t sfib_multi_tot[2][5] = { { -100 } };
        Double_t sfib_tot_fired[2][2048] = { { -100 } };
        Int_t sfib_multi_ch[2][256] = { { -100 } };
        // Double_t sfib_multi_clus_tot[2][5]={-100};
        // Double_t sfib_multi_clus_fired[2][2048]={-100};

        // attaching the fiber mapping file

        ifstream f("/u/afalduto/sfib_map_forward.txt");
        f.ignore(numeric_limits<streamsize>::max(), '\n');
        unsigned top_ch, bot_ch, fiber_no;
        /*
        while (f) {
          f >> top_ch >> bot_ch >> fiber_no;
          sfib_map[top_ch][bot_ch]=fiber_no;
        }
        */
        if ((sfib_tot_ch[0][1] / 8 == sfib_tot_ch[1][1] / 8) && sfib_tot_ch[0][1] > 0 && sfib_tot_ch[0][1] > 0)
        {
            sfib_multi++;
            sfib_fired[sfib_multi] = sfib_map[sfib_tot_ch[0][1]][sfib_tot_ch[1][1]]; // fiber ID
            sfib_tot_fired[0][sfib_fired[sfib_multi]] = sfib_tot_max[0][1];          // TOT for fired fiber top
            sfib_tot_fired[1][sfib_fired[sfib_multi]] = sfib_tot_max[1][1];          // TOT for fired fiber bot
            sfib_multi_ch[0][sfib_multi] = sfib_tot_ch[0][1];                        // SiPM ID for top
            sfib_multi_tot[0][sfib_multi] = sfib_tot_max[0][1];                      // TOT for SiPM top
            sfib_multi_ch[1][sfib_multi] = sfib_tot_ch[1][1];                        // SiPM ID for bot
            sfib_multi_tot[1][sfib_multi] = sfib_tot_max[1][1];                      // TOT for SiPM bot
        }
        // else
        //    sfib_multi == sfib_multi;
        if ((sfib_tot_ch[0][1] / 8 == sfib_tot_ch[1][2] / 8) && sfib_tot_ch[0][1] > 0 && sfib_tot_ch[1][2] > 0)
        {
            sfib_multi++;
            sfib_fired[sfib_multi] = sfib_map[sfib_tot_ch[0][1]][sfib_tot_ch[1][2]];
            sfib_tot_fired[0][sfib_fired[sfib_multi]] = sfib_tot_max[0][1];
            sfib_tot_fired[1][sfib_fired[sfib_multi]] = sfib_tot_max[1][2];
            sfib_multi_ch[0][sfib_multi] = sfib_tot_ch[0][1];
            sfib_multi_tot[0][sfib_multi] = sfib_tot_max[0][1];
            sfib_multi_ch[1][sfib_multi] = sfib_tot_ch[1][2];
            sfib_multi_tot[1][sfib_multi] = sfib_tot_max[1][2];
        }
        // else
        //    sfib_multi == sfib_multi;
        if ((sfib_tot_ch[0][2] / 8 == sfib_tot_ch[1][1] / 8) && sfib_tot_ch[0][2] > 0 && sfib_tot_ch[1][1] > 0)
        {
            sfib_multi++;
            sfib_fired[sfib_multi] = sfib_map[sfib_tot_ch[0][2]][sfib_tot_ch[1][1]];
            sfib_tot_fired[0][sfib_fired[sfib_multi]] = sfib_tot_max[0][2];
            sfib_tot_fired[1][sfib_fired[sfib_multi]] = sfib_tot_max[1][1];
            sfib_multi_ch[0][sfib_multi] = sfib_tot_ch[0][2];
            sfib_multi_tot[0][sfib_multi] = sfib_tot_max[0][2];
            sfib_multi_ch[1][sfib_multi] = sfib_tot_ch[1][1];
            sfib_multi_tot[1][sfib_multi] = sfib_tot_max[1][1];
        }
        // else
        //    sfib_multi == sfib_multi;
        if ((sfib_tot_ch[0][2] / 8 == sfib_tot_ch[1][2] / 8) && sfib_tot_ch[0][2] > 0 && sfib_tot_ch[1][2] > 0)
        {
            sfib_multi++;
            sfib_fired[sfib_multi] = sfib_map[sfib_tot_ch[0][2]][sfib_tot_ch[1][2]];
            sfib_tot_fired[0][sfib_fired[sfib_multi]] = sfib_tot_max[0][2];
            sfib_tot_fired[1][sfib_fired[sfib_multi]] = sfib_tot_max[1][2];
            sfib_multi_ch[0][sfib_multi] = sfib_tot_ch[0][2];
            sfib_multi_tot[0][sfib_multi] = sfib_tot_max[0][2];
            sfib_multi_ch[1][sfib_multi] = sfib_tot_ch[1][2];
            sfib_multi_tot[1][sfib_multi] = sfib_tot_max[1][2];
        }
        // else
        //    sfib_multi == sfib_multi;
        if ((sfib_tot_ch[0][1] / 8 == sfib_tot_ch[1][3] / 8) && sfib_tot_ch[0][1] > 0 && sfib_tot_ch[1][3] > 0)
        {
            sfib_multi++;
            sfib_fired[sfib_multi] = sfib_map[sfib_tot_ch[0][1]][sfib_tot_ch[1][3]];
            sfib_tot_fired[0][sfib_fired[sfib_multi]] = sfib_tot_max[0][1];
            sfib_tot_fired[1][sfib_fired[sfib_multi]] = sfib_tot_max[1][3];
            sfib_multi_ch[0][sfib_multi] = sfib_tot_ch[0][1];
            sfib_multi_tot[0][sfib_multi] = sfib_tot_max[0][1];
            sfib_multi_ch[1][sfib_multi] = sfib_tot_ch[1][3];
            sfib_multi_tot[1][sfib_multi] = sfib_tot_max[1][3];
        }
        // else
        //    sfib_multi == sfib_multi;
        if ((sfib_tot_ch[0][3] / 8 == sfib_tot_ch[1][1] / 8) && sfib_tot_ch[0][3] > 0 && sfib_tot_ch[1][1] > 0)
        {
            sfib_multi++;
            sfib_fired[sfib_multi] = sfib_map[sfib_tot_ch[0][3]][sfib_tot_ch[1][1]];
            sfib_tot_fired[0][sfib_fired[sfib_multi]] = sfib_tot_max[0][3];
            sfib_tot_fired[1][sfib_fired[sfib_multi]] = sfib_tot_max[1][1];
            sfib_multi_ch[0][sfib_multi] = sfib_tot_ch[0][3];
            sfib_multi_tot[0][sfib_multi] = sfib_tot_max[0][3];
            sfib_multi_ch[1][sfib_multi] = sfib_tot_ch[1][1];
            sfib_multi_tot[1][sfib_multi] = sfib_tot_max[1][1];
        }
        // else
        //    sfib_multi == sfib_multi;
        if ((sfib_tot_ch[0][2] / 8 == sfib_tot_ch[1][3] / 8) && sfib_tot_ch[0][2] > 0 && sfib_tot_ch[1][3] > 0)
        {
            sfib_multi++;
            sfib_fired[sfib_multi] = sfib_map[sfib_tot_ch[0][2]][sfib_tot_ch[1][3]];
            sfib_tot_fired[0][sfib_fired[sfib_multi]] = sfib_tot_max[0][2];
            sfib_tot_fired[1][sfib_fired[sfib_multi]] = sfib_tot_max[1][3];
            sfib_multi_ch[0][sfib_multi] = sfib_tot_ch[0][2];
            sfib_multi_tot[0][sfib_multi] = sfib_tot_max[0][2];
            sfib_multi_ch[1][sfib_multi] = sfib_tot_ch[1][3];
            sfib_multi_tot[1][sfib_multi] = sfib_tot_max[1][3];
        }
        // else
        // sfib_multi == sfib_multi;
        if ((sfib_tot_ch[0][3] / 8 == sfib_tot_ch[1][3] / 8) && sfib_tot_ch[0][3] > 0 && sfib_tot_ch[1][3] > 0)
        {
            sfib_multi++;
            sfib_fired[sfib_multi] = sfib_map[sfib_tot_ch[0][3]][sfib_tot_ch[1][3]];
            sfib_tot_fired[0][sfib_fired[sfib_multi]] = sfib_tot_max[0][3];
            sfib_tot_fired[1][sfib_fired[sfib_multi]] = sfib_tot_max[1][3];
            sfib_multi_ch[0][sfib_multi] = sfib_tot_ch[0][3];
            sfib_multi_tot[0][sfib_multi] = sfib_tot_max[0][3];
            sfib_multi_ch[1][sfib_multi] = sfib_tot_ch[1][3];
            sfib_multi_tot[1][sfib_multi] = sfib_tot_max[1][3];
        }
        // else
        //   sfib_multi == sfib_multi;
        if ((sfib_tot_ch[0][3] / 8 == sfib_tot_ch[1][2] / 8) && sfib_tot_ch[0][3] > 0 && sfib_tot_ch[1][2] > 0)
        {
            sfib_multi++;
            sfib_fired[sfib_multi] = sfib_map[sfib_tot_ch[0][3]][sfib_tot_ch[1][2]];
            sfib_tot_fired[0][sfib_fired[sfib_multi]] = sfib_tot_max[0][3];
            sfib_tot_fired[1][sfib_fired[sfib_multi]] = sfib_tot_max[1][2];
            sfib_multi_ch[0][sfib_multi] = sfib_tot_ch[0][3];
            sfib_multi_tot[0][sfib_multi] = sfib_tot_max[0][3];
            sfib_multi_ch[1][sfib_multi] = sfib_tot_ch[1][2];
            sfib_multi_tot[1][sfib_multi] = sfib_tot_max[1][2];
        }
        // else
        //   sfib_multi == sfib_multi;
        if ((sfib_tot_ch[0][1] / 8 == sfib_tot_ch[1][4] / 8) && sfib_tot_ch[0][1] > 0 && sfib_tot_ch[1][4] > 0)
        {
            sfib_multi++;
            sfib_fired[sfib_multi] = sfib_map[sfib_tot_ch[0][3]][sfib_tot_ch[1][4]];
            sfib_tot_fired[0][sfib_fired[sfib_multi]] = sfib_tot_max[0][1];
            sfib_tot_fired[1][sfib_fired[sfib_multi]] = sfib_tot_max[1][4];
            sfib_multi_ch[0][sfib_multi] = sfib_tot_ch[0][1];
            sfib_multi_tot[0][sfib_multi] = sfib_tot_max[0][1];
            sfib_multi_ch[1][sfib_multi] = sfib_tot_ch[1][4];
            sfib_multi_tot[1][sfib_multi] = sfib_tot_max[1][4];
        }
        // else
        //    sfib_multi == sfib_multi;
        if ((sfib_tot_ch[0][2] / 8 == sfib_tot_ch[1][4] / 8) && sfib_tot_ch[0][2] > 0 && sfib_tot_ch[1][4] > 0)
        {
            sfib_multi++;
            sfib_fired[sfib_multi] = sfib_map[sfib_tot_ch[0][2]][sfib_tot_ch[1][4]];
            sfib_tot_fired[0][sfib_fired[sfib_multi]] = sfib_tot_max[0][2];
            sfib_tot_fired[1][sfib_fired[sfib_multi]] = sfib_tot_max[1][4];
            sfib_multi_ch[0][sfib_multi] = sfib_tot_ch[0][2];
            sfib_multi_tot[0][sfib_multi] = sfib_tot_max[0][2];
            sfib_multi_ch[1][sfib_multi] = sfib_tot_ch[1][4];
            sfib_multi_tot[1][sfib_multi] = sfib_tot_max[1][4];
        }
        // else
        //    sfib_multi == sfib_multi;
        if ((sfib_tot_ch[0][3] / 8 == sfib_tot_ch[1][4] / 8) && sfib_tot_ch[0][3] > 0 && sfib_tot_ch[1][4] > 0)
        {
            sfib_multi++;
            sfib_fired[sfib_multi] = sfib_map[sfib_tot_ch[0][3]][sfib_tot_ch[1][4]];
            sfib_tot_fired[0][sfib_fired[sfib_multi]] = sfib_tot_max[0][3];
            sfib_tot_fired[1][sfib_fired[sfib_multi]] = sfib_tot_max[1][4];
            sfib_multi_ch[0][sfib_multi] = sfib_tot_ch[0][3];
            sfib_multi_tot[0][sfib_multi] = sfib_tot_max[0][3];
            sfib_multi_ch[1][sfib_multi] = sfib_tot_ch[1][4];
            sfib_multi_tot[1][sfib_multi] = sfib_tot_max[1][4];
        }
        // else
        //    sfib_multi == sfib_multi;
        if ((sfib_tot_ch[0][4] / 8 == sfib_tot_ch[1][3] / 8) && sfib_tot_ch[0][4] > 0 && sfib_tot_ch[1][3] > 0)
        {
            sfib_multi++;
            sfib_fired[sfib_multi] = sfib_map[sfib_tot_ch[0][4]][sfib_tot_ch[1][3]];
            sfib_tot_fired[0][sfib_fired[sfib_multi]] = sfib_tot_max[0][4];
            sfib_tot_fired[1][sfib_fired[sfib_multi]] = sfib_tot_max[1][3];
            sfib_multi_ch[0][sfib_multi] = sfib_tot_ch[0][4];
            sfib_multi_tot[0][sfib_multi] = sfib_tot_max[0][4];
            sfib_multi_ch[1][sfib_multi] = sfib_tot_ch[1][3];
            sfib_multi_tot[1][sfib_multi] = sfib_tot_max[1][3];
        }
        // else
        //    sfib_multi == sfib_multi;
        if ((sfib_tot_ch[0][4] / 8 == sfib_tot_ch[1][1] / 8) && sfib_tot_ch[0][4] > 0 && sfib_tot_ch[1][1] > 0)
        {
            sfib_multi++;
            sfib_fired[sfib_multi] = sfib_map[sfib_tot_ch[0][4]][sfib_tot_ch[1][1]];
            sfib_tot_fired[0][sfib_fired[sfib_multi]] = sfib_tot_max[0][4];
            sfib_tot_fired[1][sfib_fired[sfib_multi]] = sfib_tot_max[1][1];
            sfib_multi_ch[0][sfib_multi] = sfib_tot_ch[0][4];
            sfib_multi_tot[0][sfib_multi] = sfib_tot_max[0][4];
            sfib_multi_ch[1][sfib_multi] = sfib_tot_ch[1][1];
            sfib_multi_tot[1][sfib_multi] = sfib_tot_max[1][1];
        }
        // else
        //   sfib_multi == sfib_multi;
        if ((sfib_tot_ch[0][4] / 8 == sfib_tot_ch[1][2] / 8) && sfib_tot_ch[0][4] > 0 && sfib_tot_ch[1][2] > 0)
        {
            sfib_multi++;
            sfib_fired[sfib_multi] = sfib_map[sfib_tot_ch[0][4]][sfib_tot_ch[1][2]];
            sfib_tot_fired[0][sfib_fired[sfib_multi]] = sfib_tot_max[0][4];
            sfib_tot_fired[1][sfib_fired[sfib_multi]] = sfib_tot_max[1][2];
            sfib_multi_ch[0][sfib_multi] = sfib_tot_ch[0][4];
            sfib_multi_tot[0][sfib_multi] = sfib_tot_max[0][4];
            sfib_multi_ch[1][sfib_multi] = sfib_tot_ch[1][2];
            sfib_multi_tot[1][sfib_multi] = sfib_tot_max[1][2];
        }
        // else
        //    sfib_multi == sfib_multi;
        if ((sfib_tot_ch[0][4] / 8 == sfib_tot_ch[1][4] / 8) && sfib_tot_ch[0][4] > 0 && sfib_tot_ch[1][4] > 0)
        {
            sfib_multi++;
            sfib_fired[sfib_multi] = sfib_map[sfib_tot_ch[0][4]][sfib_tot_ch[1][4]];
            sfib_tot_fired[0][sfib_fired[sfib_multi]] = sfib_tot_max[0][4];
            sfib_tot_fired[1][sfib_fired[sfib_multi]] = sfib_tot_max[1][4];
            sfib_multi_ch[0][sfib_multi] = sfib_tot_ch[0][4];
            sfib_multi_tot[0][sfib_multi] = sfib_tot_max[0][4];
            sfib_multi_ch[1][sfib_multi] = sfib_tot_ch[1][4];
            sfib_multi_tot[1][sfib_multi] = sfib_tot_max[1][4];
        }
        // else
        //    sfib_multi == sfib_multi;

        // cout << sfib_multi<<endl; //Top
        //
        // clustering analysis start from here
        //
        Int_t sfib_clus = 0;
        //
        // (1) for multiplicity 1
        //
        /*    if(sfib_multi==1)
          {
        sfib_clus=1;
        sfib_multi_clus_ch[0][1][1]=sfib_multi_ch[0][1];
        sfib_multi_clus_tot[0][1][1]=sfib_multi_tot[0][1];
        sfib_multi_clus_ch[1][1][1]=sfib_multi_ch[1][1];
        sfib_multi_clus_tot[1][1][1]=sfib_multi_tot[1][1];
        sfib_multi_clus_fired[1][1]=sfib_fired[1];
        }
        //cout << sfib_clus << endl;
        //
        */
        // (2) for multiplicity 2
        //
        if (sfib_multi == 2)
        {
            if ((sfib_fired[1] - sfib_fired[2] == 1) || (sfib_fired[2] - sfib_fired[1] == 1))
                sfib_clus = 1;
            else
                sfib_clus = 2;
        }
        //
        // (3) for multiplicity 3
        //
        if (sfib_multi == 3)
        {
            if (abs(sfib_fired[1] - sfib_fired[2]) == 1)
            {
                if (abs(sfib_fired[1] - sfib_fired[3]) == 1 || abs(sfib_fired[2] - sfib_fired[3]) == 1)
                    sfib_clus = 1;
                // condition 1-2-3, 3-1-2, 2-1-3, 3-2-1
                else
                    sfib_clus = 2;
            }
            else if (abs(sfib_fired[1] - sfib_fired[3]) == 1)
            {
                if (abs(sfib_fired[2] - sfib_fired[3]) == 1)
                    sfib_clus = 1;
                // condition 1-3-2, 2-3-1
                else
                    sfib_clus = 2;
            }
            else if (abs(sfib_fired[2] - sfib_fired[3]) == 1)
                sfib_clus = 2;
            else
                sfib_clus = 3;
        }
        //
        // (4) for multiplicity 4
        //
        /*   if(sfib_multi==4)
          {
        if(abs(sfib_fired[1]-sfib_fired[3])==1)
          {
            if(abs(sfib_fired[1]-sfib_fired[4])==1||abs(sfib_fired[2]-sfib_fired[4])==1||abs(sfib_fired[3]-sfib_fired[4])==1)
        sfib_clus=1;
        */

        // cout << sfib_clus<<endl; //Top

        fh_sfib_tot_ch_top1->Fill(sfib_tot_ch[0][1]);
        fh_sfib_tot_max_top1->Fill(sfib_tot_max[0][1]);
        fh_sfib_tot_ch_bot1->Fill(sfib_tot_ch[1][1]);
        fh_sfib_tot_max_bot1->Fill(sfib_tot_max[1][1]);
        fh_sfib_tot_ch_top2->Fill(sfib_tot_ch[0][2]);
        fh_sfib_tot_max_top2->Fill(sfib_tot_max[0][2]);
        fh_sfib_tot_ch_bot2->Fill(sfib_tot_ch[1][2]);
        fh_sfib_tot_max_bot2->Fill(sfib_tot_max[1][2]);

        fh_sfib_tot_ch_max_top1->Fill(sfib_tot_ch[0][1], sfib_tot_max[0][1]);
        fh_sfib_tot_ch_max_top2->Fill(sfib_tot_ch[0][2], sfib_tot_max[0][2]);
        fh_sfib_tot_ch_max_bot1->Fill(sfib_tot_ch[1][1], sfib_tot_max[1][1]);
        fh_sfib_tot_ch_max_bot2->Fill(sfib_tot_ch[1][2], sfib_tot_max[1][2]);

        fh_sfib_tot_ch_topbot1->Fill(sfib_tot_ch[0][1], sfib_tot_ch[1][1]);

        fh_sfib_multi->Fill(sfib_multi);

        fh_sfib_clus->Fill(sfib_clus);
        fh_sfib_multi_clus->Fill(sfib_multi, sfib_clus);
        fh_sfib_multi_ch_topbot1->Fill(sfib_multi_ch[0][1], sfib_multi_ch[1][1]);
        fh_sfib_tot_top_multi1->Fill(sfib_multi_tot[0][1]);
        fh_sfib_tot_top_multi2->Fill(sfib_multi_tot[0][2]);

        // 2019.03.11 ASCII OUTPUT J.Tanaka
        if (sfib_tot_max[0][1] > 0 || sfib_tot_max[0][1] > 0)
        {
            //	cout << sfib_tot_max[0][1]<<" "<< sfib_tot_max[0][2]<<" "<< sfib_tot_max[0][3]<<endl; //Top
            /*	std::ofstream run_max202("run_max202.txt",std::ios::app);
            run_max202 << fNEvents<<" ";
            run_max202 << sfib_tot_ch[0][1]<<" "; //Top 1st
            run_max202 << sfib_tot_max[0][1]<<" "; //Top 1st
            run_max202 << sfib_tot_ch[0][2]<<" "; //Top 2nd
            run_max202 << sfib_tot_max[0][2]<<" "; //Top 2nd
            run_max202 << sfib_tot_ch[0][3]<<" "; //Top 1st
            run_max202 << sfib_tot_max[0][3]<<" "; //Top 1st
            run_max202 << sfib_tot_ch[0][4]<<" "; //Top 2nd
            run_max202 << sfib_tot_max[0][4]<<" "; //Top 2nd
            run_max202 << sfib_tot_ch[0][5]<<" "; //Top 2nd
            run_max202 << sfib_tot_max[0][5]<<" "; //Top 2nd
            run_max202 << sfib_tot_ch[1][1]<<" "; //Bot 1st
            run_max202 << sfib_tot_max[1][1]<<" "; //Bot 1st
            run_max202 << sfib_tot_ch[1][2]<<" "; //Bot 2nd
            run_max202 << sfib_tot_max[1][2]<<" "; //Bot 2nd
            run_max202 << sfib_tot_ch[1][3]<<" "; //Bot 3rd
            run_max202 << sfib_tot_max[1][3]<<" "; //Bot 3rd
            run_max202 << sfib_tot_ch[1][4]<<" "; //Bot 4th
            run_max202 << sfib_tot_max[1][4]<<" "; //Bot 4th
            run_max202 << sfib_tot_ch[1][5]<<" "; //Bot 4th
            run_max202 << sfib_tot_max[1][5]<<" "; //Bot 4th
            run_max202<<std::endl;
            run_max202.close();*/
        }
    }
    fNEvents += 1;
}

void R3BOnlineSpectraSfib::FinishEvent()
{

    for (Int_t det = 0; det < DET_MAX; det++)
    {
        if (fMappedItems.at(det))
        {
            fMappedItems.at(det)->Clear();
        }
        if (fCalItems.at(det))
        {
            fCalItems.at(det)->Clear();
        }
        if (fHitItems.at(det))
        {
            fHitItems.at(det)->Clear();
        }
    }
}

void R3BOnlineSpectraSfib::FinishTask()
{
    if (fMappedItems.at(DET_SFIB))
    {
        fh_sfib_Tot_top->Write();
        fh_sfib_Tot_top1d->Write();
        fh_sfib_Tot_top_g->Write();
        fh_sfib_Tot_top1d_g->Write();
        fh_sfib_Tot_bot->Write();
        fh_sfib_Tot_bot1d->Write();
        fh_sfib_Tot_bot_g->Write();
        fh_sfib_Tot_bot1d_g->Write();
        fh_sfib_channels->Write();
        fh_sfib_channels_g->Write();
        // fh_sfib_tofd->Write();
        // fh_sfib_pspx->Write();
        fh_sfib_tot_tvb->Write();
        fh_sfib_fmult->Write();
        fh_sfib_cmult->Write();
        fh_sfib_fmult_g->Write();
        fh_sfib_cmult_g->Write();
        fh_sfib_tot_ch_top1->Write();
        fh_sfib_tot_ch_top2->Write();
        fh_sfib_tot_ch_bot1->Write();
        fh_sfib_tot_ch_bot2->Write();
        fh_sfib_tot_max_top1->Write();
        fh_sfib_tot_max_top2->Write();
        fh_sfib_tot_max_bot1->Write();
        fh_sfib_tot_max_bot2->Write();
        fh_sfib_tot_ch_max_top1->Write();
        fh_sfib_tot_ch_max_top2->Write();
        fh_sfib_tot_ch_max_bot1->Write();
        fh_sfib_tot_ch_max_bot2->Write();
        fh_sfib_tot_ch_topbot1->Write();
        fh_sfib_multi->Write();
        fh_sfib_clus->Write();
        fh_sfib_multi_clus->Write();
        fh_sfib_multi_ch_topbot1->Write();
        fh_sfib_tot_top_multi1->Write();
        fh_sfib_tot_top_multi2->Write();
        fh_sfib_channels_top->Write();
        fh_sfib_channels_bot->Write();
        fh_sfib_channels_topvsbot->Write();
    }
}

ClassImp(R3BOnlineSpectraSfib)
