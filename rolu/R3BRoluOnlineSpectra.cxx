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
// -----                R3BRoluOnlineSpectra               -----
// -----          Created April 13th 2016 by M.Heil       -----
// -----       Updated May 11th 2022 by J.L. Rodriguez    -----
// ------------------------------------------------------------

#include "R3BRoluOnlineSpectra.h"
#include "R3BCoarseTimeStitch.h"
#include "R3BEventHeader.h"
#include "R3BLogger.h"
#include "R3BRoluCalData.h"
#include "R3BRoluMappedData.h"
#include "R3BTCalEngine.h"

#include <FairLogger.h>
#include <FairRootManager.h>
#include <FairRunOnline.h>
#include <FairRuntimeDb.h>

#include <TCanvas.h>
#include <TClonesArray.h>
#include <TH1F.h>
#include <TH2F.h>
#include <THttpServer.h>
#include <TMath.h>
#include <boost/multi_array.hpp>
#include <vector>

#define IS_NAN(x) TMath::IsNaN(x)
using namespace std;
namespace
{
    double c_range_ns = 2048 * 5;
} // namespace

R3BRoluOnlineSpectra::R3BRoluOnlineSpectra()
    : R3BRoluOnlineSpectra("RoluOnlineSpectra", 1)
{
}

R3BRoluOnlineSpectra::R3BRoluOnlineSpectra(const char* name, int iVerbose)
    : FairTask(name, iVerbose)
{
}

R3BRoluOnlineSpectra::~R3BRoluOnlineSpectra()
{
    delete fMappedItems;
    delete fCalItems;
}

InitStatus R3BRoluOnlineSpectra::Init()
{
    R3BLOG(info, "");

    FairRootManager* mgr = FairRootManager::Instance();
    R3BLOG_IF(fatal, NULL == mgr, "FairRootManager not found");

    header = dynamic_cast<R3BEventHeader*>(mgr->GetObject("EventHeader."));

    FairRunOnline* run = FairRunOnline::Instance();
    run->GetHttpServer()->Register("", this);

    // Get objects for detectors on all levels

    fMappedItems = dynamic_cast<TClonesArray*>(mgr->GetObject("RoluMapped"));
    R3BLOG_IF(fatal, NULL == fMappedItems, "RoluMapped not found");

    fCalItems = dynamic_cast<TClonesArray*>(mgr->GetObject("RoluCal"));
    R3BLOG_IF(warn, NULL == fCalItems, "RoluCal not found");

    //------------------------------------------------------------------------
    // create histograms of all detectors
    //------------------------------------------------------------------------

    // MAIN FOLDER-MWPC
    TFolder* mainfol = new TFolder("ROLU", "ROLU info");

    //------------------------------------------------------------------------
    // Rolu detector
    TCanvas* cRolu[fNofRoluDetectors];
    TString histName;
    TString canName;
    if (fMappedItems)
    {
        for (int irolucount = 0; irolucount < fNofRoluDetectors; irolucount++)
        {
            histName = Form("ROLU%d", irolucount + 1);
            cRolu[irolucount] = new TCanvas(canName.Data(), canName.Data(), 10, 10, 1010, 810);

            histName = Form("Rolu%d_channels", irolucount);
            fh_rolu_channels[irolucount] = new TH1F(histName.Data(), histName.Data(), 5, 0., 5.);
            fh_rolu_channels[irolucount]->GetXaxis()->SetTitle("Channel number");
            fh_rolu_channels[irolucount]->SetFillColor(31);

            histName = Form("Rolu%d_tot", irolucount);
            fh_rolu_tot[irolucount] = new TH2F(histName.Data(), histName.Data(), 5, 0, 5, 1500, 0., 300.);
            fh_rolu_tot[irolucount]->GetXaxis()->SetTitle("PMT number");
            fh_rolu_tot[irolucount]->GetYaxis()->SetTitle("ToT / ns");

            int color;
            for (int ichannelcount = 0; ichannelcount < 4; ichannelcount++)
            {
                color = ichannelcount + 1;
                if (color == 3)
                    color = 8;
                histName = Form("Rolu%d_tot1D_ch%d", irolucount, ichannelcount);
                fh_rolu_tot_1D[irolucount][ichannelcount] = new TH1F(histName.Data(), histName.Data(), 151, 0., 150.);
                fh_rolu_tot_1D[irolucount][ichannelcount]->GetXaxis()->SetTitle("ToT / ns");
                fh_rolu_tot_1D[irolucount][ichannelcount]->SetLineColor(color);

                histName = Form("Rolu%d_LE_raw_ch%d", irolucount, ichannelcount);
                fh_rolu_LE_raw[irolucount][ichannelcount] = new TH1F(histName.Data(), histName.Data(), 1050, 0, 10500);
                fh_rolu_LE_raw[irolucount][ichannelcount]->GetXaxis()->SetTitle("ToT / ns");
                fh_rolu_LE_raw[irolucount][ichannelcount]->SetLineColor(color);
            }

            cRolu[irolucount]->Divide(2, 2);
            cRolu[irolucount]->cd(1);
            fh_rolu_channels[irolucount]->Draw();
            gPad->SetLogy();
            cRolu[irolucount]->cd(2);
            gPad->SetLogz();
            fh_rolu_tot[irolucount]->Draw("colz");
            for (int ichannelcount = 0; ichannelcount < 4; ichannelcount++)
            {
                cRolu[irolucount]->cd(3);
                if (ichannelcount == 0)
                    fh_rolu_tot_1D[irolucount][ichannelcount]->Draw();
                else
                    fh_rolu_tot_1D[irolucount][ichannelcount]->Draw("same");
                gPad->SetLogy();
                cRolu[irolucount]->cd(4);
                gPad->SetLogy();
                if (ichannelcount == 0)
                    fh_rolu_LE_raw[irolucount][ichannelcount]->Draw();
                else
                    fh_rolu_LE_raw[irolucount][ichannelcount]->Draw("same");
            }
            mainfol->Add(cRolu[irolucount]);
        }

        run->AddObject(mainfol);
        run->GetHttpServer()->RegisterCommand("Reset_ROLU_HIST", Form("/Objects/%s/->Reset_ROLU_Histo()", GetName()));
    }

    return kSUCCESS;
}

void R3BRoluOnlineSpectra::Reset_ROLU_Histo()
{
    R3BLOG(info, "");
    if (fMappedItems)
    {
        for (int irolucount = 0; irolucount < fNofRoluDetectors; irolucount++)
        {
            fh_rolu_channels[irolucount]->Reset();
            fh_rolu_tot[irolucount]->Reset();
            for (int ichannelcount = 0; ichannelcount < 4; ichannelcount++)
            {
                fh_rolu_tot_1D[irolucount][ichannelcount]->Reset();
                fh_rolu_LE_raw[irolucount][ichannelcount]->Reset();
            }
        }
    }
}

void R3BRoluOnlineSpectra::Exec(Option_t* option)
{
    if ((fTrigger >= 0) && (header) && (header->GetTrigger() != fTrigger))
        return;

    if ((fTpat > 0) && (header) && ((header->GetTpat() & fTpat) != fTpat))
        return;

    int nPartROLU = 0;
    double tot = 0;

    if (fMappedItems)
    {
        int nHits = fMappedItems->GetEntriesFast();

        if (nHits > 0)
            nRoluEvents += 1;

        for (int ihit = 0; ihit < nHits; ihit++)
        {
            R3BRoluMappedData* hit = dynamic_cast<R3BRoluMappedData*>(fMappedItems->At(ihit));
            if (!hit)
                continue;

            // channel numbers are stored 1-based (1..n)
            fh_rolu_channels[hit->GetDetector() - 1]->Fill(hit->GetChannel());
        }
    }

    if (fCalItems)
    {
        nPartROLU = fCalItems->GetEntriesFast();
        int iDet = 0;

        for (int iPart = 0; iPart < nPartROLU; iPart++)
        {
            R3BRoluCalData* calData = dynamic_cast<R3BRoluCalData*>(fCalItems->At(iPart));
            assert(calData && "dynamic cast failed!");
            for (int chan = 0; chan < 4; chan++)
            // if (!(IS_NAN(calData->GetTimeL_ns(calData->GetTAMEXLNcha()-1))) &&
            // !(IS_NAN(calData->GetTimeL_ns(calData->GetTAMEXLNcha()-1))))
            { // TAMEX leading
                double time_L = calData->GetTimeL_ns(chan);
                double time_T = calData->GetTimeT_ns(chan);
                tot = fmod(time_T - time_L + c_range_ns + c_range_ns / 2., c_range_ns) - c_range_ns / 2.;
                if (!std::isnan(tot))
                {
                    fh_rolu_tot[calData->GetDetector() - 1]->Fill(chan + 1, tot);
                    fh_rolu_tot_1D[calData->GetDetector() - 1][chan]->Fill(tot);
                    fh_rolu_LE_raw[calData->GetDetector() - 1][chan]->Fill(calData->GetTimeL_ns(chan));
                }
            }
        }
    } // if fCallItems
    double maxentries = 0;
    for (int i = 0; i < 4; i++)
        maxentries = max(maxentries, fh_rolu_tot_1D[0][i]->GetMaximum());
    fh_rolu_tot_1D[0][0]->SetMaximum(maxentries);
}

void R3BRoluOnlineSpectra::FinishEvent()
{
    if (fMappedItems)
    {
        fMappedItems->Clear();
    }
    if (fCalItems)
    {
        fCalItems->Clear();
    }
}

void R3BRoluOnlineSpectra::FinishTask()
{
    if (fMappedItems)
    {
        for (int irolucount = 0; irolucount < fNofRoluDetectors; irolucount++)
        {
            fh_rolu_channels[irolucount]->Write();
            fh_rolu_tot[irolucount]->Write();
            for (int ichannelcount = 0; ichannelcount < 4; ichannelcount++)
            {
                fh_rolu_tot_1D[irolucount][ichannelcount]->Write();
                fh_rolu_LE_raw[irolucount][ichannelcount]->Write();
            }
        }
    }
}

ClassImp(R3BRoluOnlineSpectra)
