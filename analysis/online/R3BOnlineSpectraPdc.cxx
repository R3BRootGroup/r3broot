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
// -----              R3BOnlineSpectraPdc                 -----
// -----          Created July  2020 by M.Heil            -----
// ------------------------------------------------------------

/*
 * This task should fill histograms with detector variables which allow
 * to test the detectors online
 *
 */

#include "R3BOnlineSpectraPdc.h"

#include "R3BPdcCalData.h"
#include "R3BPdcHitData.h"
#include "R3BPdcMappedData.h"

#include "R3BEventHeader.h"
#include "R3BTCalEngine.h"

#include "R3BSamplerMappedData.h"

#include "FairLogger.h"
#include "FairRootManager.h"
#include "FairRunAna.h"
#include "FairRunOnline.h"
#include "FairRuntimeDb.h"
#include "TCanvas.h"
#include "TGaxis.h"
#include "TH1F.h"
#include "TH2F.h"
#include "THttpServer.h"
#include "TLegend.h"

#include "TClonesArray.h"
#include "TMath.h"
#include <TRandom3.h>
#include <TRandomGen.h>
#include <algorithm>
#include <array>
#include <cstdlib>
#include <ctime>
#include <fstream>
#include <iostream>
#include <sstream>
#include <vector>

#include <assert.h>
#define IS_NAN(x) TMath::IsNaN(x)
using namespace std;

R3BOnlineSpectraPdc::R3BOnlineSpectraPdc()
    : R3BOnlineSpectraPdc("OnlineSpectra", 1)
{
}

R3BOnlineSpectraPdc::R3BOnlineSpectraPdc(const char* name, Int_t iVerbose)
    : FairTask(name, iVerbose)
    , fTrigger(-1)
    , fTpat(-1)
    , fSamp(-1)
    , fNEvents(0)
{
}

R3BOnlineSpectraPdc::~R3BOnlineSpectraPdc()
{
    //	delete fhTpat ;
    //	delete fhTrigger;
    //	delete fh_SEETRAM;

    for (int i = 0; i < N_PLANE_MAX_PDC; i++)
    {
        delete fh_Pdc_Wire[i];
        delete fh_Pdc_Tot[i];
        delete fh_Pdc_ToF[i];
        delete fh_Pdc_xy[i];
        delete fh_Pdc_Time[i];
        delete fh_Pdc_Wire_vs_Events[i];
    }
}

InitStatus R3BOnlineSpectraPdc::Init()
{

    // Initialize random number:
    std::srand(std::time(0)); // use current time as seed for random generator

    LOG(info) << "R3BOnlineSpectraPdc::Init ";

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
    fMappedItems = dynamic_cast<TClonesArray*>(mgr->GetObject("PdcMapped"));
    if (NULL == fMappedItems)
        LOG(fatal) << "Branch PdcCal not found";

    fCalItems = dynamic_cast<TClonesArray*>(mgr->GetObject("PdcCal"));
    if (NULL == fCalItems)
        LOG(fatal) << "Branch PdcCal not found";

    fHitItems = dynamic_cast<TClonesArray*>(mgr->GetObject("PdcHit"));
    if (NULL == fHitItems)
        LOG(fatal) << "Branch PdcCal not found";

    //------------------------------------------------------------------------
    // create histograms of all detectors
    //------------------------------------------------------------------------

    // Trigger and Tpat
    TCanvas* cTrigg = new TCanvas("Trigg", "Triggers", 10, 10, 650, 350);
    fhTpat = new TH1F("Tpat", "Tpat", 20, 0, 20);
    fhTpat->GetXaxis()->SetTitle("Tpat value");

    fhTrigger = new TH1F("Trigger", "Trigger all", 20, 0, 20);
    fhTrigger->GetXaxis()->SetTitle("Trigger value");

    cTrigg->Divide(2, 1);
    cTrigg->cd(1);
    gPad->SetLogy();
    fhTrigger->Draw();
    cTrigg->cd(2);
    gPad->SetLogy();
    fhTpat->Draw();
    cTrigg->cd(0);

    //---------------------------------------------------------------------------------------------------
    // PDC detector

    if (fMappedItems)
    {
        TCanvas* cPdc_planes = new TCanvas("PDC_planes", "PDC planes", 10, 10, 1100, 1000);
        cPdc_planes->Divide(4, 4);

        for (Int_t j = 0; j < N_PLANE_MAX_PDC; j++)
        {
            char strName1[255];
            sprintf(strName1, "pdc_wire_plane_%d", j + 1);
            char strName2[255];
            sprintf(strName2, "PDC wires plane %d", j + 1);
            fh_Pdc_Wire[j] = new TH1F(strName1, strName2, 150, 0., 150.);
            fh_Pdc_Wire[j]->GetXaxis()->SetTitle("Wire");
            fh_Pdc_Wire[j]->GetYaxis()->SetTitle("Counts");

            char strName3[255];
            sprintf(strName3, "pdc_ToT_plane_%d", j + 1);
            char strName4[255];
            sprintf(strName4, "PDC ToT plane %d", j + 1);
            fh_Pdc_Tot[j] = new TH2F(strName3, strName4, 150, 0, 150, 100, 0., 100.);
            fh_Pdc_Tot[j]->GetXaxis()->SetTitle("Wire");
            fh_Pdc_Tot[j]->GetYaxis()->SetTitle("ToT / ns");

            char strName5[255];
            sprintf(strName5, "ped_ToF_plane_%d", j + 1);
            char strName6[255];
            sprintf(strName6, "PDC ToF plane %d", j + 1);
            fh_Pdc_ToF[j] = new TH2F(strName5, strName6, 150, 0, 150, 1000, 0, 100);
            fh_Pdc_ToF[j]->GetXaxis()->SetTitle("Wire");
            fh_Pdc_ToF[j]->GetYaxis()->SetTitle("ToF / ns");

            char strName7[255];
            sprintf(strName7, "pdc_time_plane_%d", j + 1);
            char strName8[255];
            sprintf(strName8, "PDC Time-trigger time plane %d", j + 1);
            fh_Pdc_Time[j] = new TH2F(strName7, strName8, 150, 0, 150, 1000, -1000, 1000);
            fh_Pdc_Time[j]->GetXaxis()->SetTitle("Wire");
            fh_Pdc_Time[j]->GetYaxis()->SetTitle("Time / ns");

            char strName9[255];
            sprintf(strName9, "pdc_xy_plane_%d", j + 1);
            char strName10[255];
            sprintf(strName10, "PDC xy plane %d", j + 1);
            fh_Pdc_xy[j] = new TH2F(strName9, strName10, 1500, 0, 150, 1000, 0, 100);
            fh_Pdc_xy[j]->GetXaxis()->SetTitle("x in cm");
            fh_Pdc_xy[j]->GetYaxis()->SetTitle("y in cm");

            /*
                        char strName9[255];
                        sprintf(strName9, "tofd_multihit_plane_%d", j + 1);
                        char strName10[255];
                        sprintf(strName10, "Tofd multihit plane %d", j + 1);
                        fh_tofd_multihit[j] = new TH2F(strName9, strName10, 45, 0., 45., 10, 0, 10);
                        fh_tofd_multihit[j]->GetXaxis()->SetTitle("Bar number");
                        fh_tofd_multihit[j]->GetYaxis()->SetTitle("Multihit");
            */
        }

        for (Int_t j = 0; j < N_PLANE_MAX_PDC; j++)
        {
            cPdc_planes->cd(j * N_PLANE_MAX_PDC + 1);
            fh_Pdc_Wire[j]->Draw();
            cPdc_planes->cd(j * N_PLANE_MAX_PDC + 2);
            // gPad->SetLogz();
            fh_Pdc_Tot[j]->Draw("colz");
            cPdc_planes->cd(j * N_PLANE_MAX_PDC + 3);
            // gPad->SetLogz();
            fh_Pdc_Time[j]->Draw("colz");
            cPdc_planes->cd(j * N_PLANE_MAX_PDC + 4);
            // gPad->SetLogz();
            fh_Pdc_xy[j]->Draw("colz");
        }

        cPdc_planes->cd(0);
        run->AddObject(cPdc_planes);

        run->GetHttpServer()->RegisterCommand("Reset_PDC", Form("/Tasks/%s/->Reset_PDC_Histo()", GetName()));
    }

    return kSUCCESS;
}

void R3BOnlineSpectraPdc::Reset_PDC_Histo()
{
    for (int i = 0; i < N_PLANE_MAX_PDC; i++)
    {
        fh_Pdc_Wire[i]->Reset();
        fh_Pdc_Time[i]->Reset();
        fh_Pdc_ToF[i]->Reset();
        fh_Pdc_Tot[i]->Reset();
        fh_Pdc_xy[i]->Reset();
    }
}

void R3BOnlineSpectraPdc::Exec(Option_t* option)
{
    FairRootManager* mgr = FairRootManager::Instance();
    if (NULL == mgr)
    {
        // FairLogger::GetLogger()->Fatal(MESSAGE_ORIGIN, "FairRootManager not found");
        LOG(error) << "FairRootManager not found";
        return;
    }
    /*
        if (fMappedItems)
        {
            fhTrigger->Fill(header->GetTrigger());
        }

        //   check for requested trigger (Todo: should be done globablly / somewhere else)
        if ((fTrigger >= 0) && (header) && (header->GetTrigger() != fTrigger))
            return;

        if (fMappedItems)
        {
            Int_t tpatbin;
            for (int i = 0; i < 16; i++)
            {
                tpatbin = (header->GetTpat() & (1 << i));
                if (tpatbin != 0)
                    fhTpat->Fill(i + 1);
            }
        }

        // fTpat = 1-16; fTpat_bit = 0-15
        Int_t fTpat_bit = fTpat - 1;
        Int_t itpat;
        Int_t tpatvalue;
        if (fTpat_bit >= 0)
        {
            itpat = header->GetTpat();
            tpatvalue = (itpat & (1 << fTpat_bit)) >> fTpat_bit;
            if (tpatvalue == 0)
                return;
        }
    */
    //----------------------------------------------------------------------
    // PDC
    //----------------------------------------------------------------------

    if (fMappedItems)
    {
        // cout << "MappedItems" << endl;
        auto det = fMappedItems;
        Int_t nMapped = det->GetEntriesFast();
        for (Int_t imapped = 0; imapped < nMapped; imapped++)
        {
            auto mapped = dynamic_cast<R3BPdcMappedData const*>(det->At(imapped));
            if (!mapped)
                continue; // should not happen

            Int_t const iPlane = mapped->GetPlaneId(); // 1..n
            Int_t const iWire = mapped->GetWireId();   // 1..n
            Int_t const iEdge = mapped->GetEdgeId();
            // cout << "Plane: " << iPlane << " Wire: " << iWire << endl;
            if (iPlane <= N_PLANE_MAX_PDC)
            {
                fh_Pdc_Wire[iPlane - 1]->Fill(iWire);
            }
        }
    }

    if (fCalItems)
    {
        // cout << "CalItems" << endl;
        auto det = fCalItems;
        Int_t nCals = det->GetEntriesFast();

        for (Int_t ical = 0; ical < nCals; ical++)
        {
            auto cal = dynamic_cast<R3BPdcCalData const*>(det->At(ical));
            if (!cal)
                continue; // should not happen

            Int_t plane = cal->GetPlaneId();
            Double_t wire = cal->GetWireId();
            Double_t edge = cal->GetEdgeId();
            Double_t t = cal->GetTime_ns();

            // cout << "Plane: " << plane << " Wire: " << wire << endl;
            // fh_Pdc_Time[plane-1]->Fill(wire,t);
        }
    }

    if (fHitItems)
    {
        // cout << "HitItems" << endl;
        auto det = fHitItems;
        Int_t nHits = det->GetEntriesFast();
        Double_t t0 = -10000.;
        Int_t plane = 0;
        Int_t wire = 0;
        Double_t x[4], y[4], eloss[4], t[4];
        for (Int_t i = 0; i < 4; i++)
        {
            x[i] = -10000.;
            y[i] = -10000;
            eloss[i] = 0.;
            t[i] = -10000.;
        }

        for (Int_t ihit = 0; ihit < nHits; ihit++)
        {
            auto hit = dynamic_cast<R3BPdcHitData const*>(det->At(ihit));
            if (!hit)
                continue; // should not happen
            if (t0 == -10000.)
                t0 = hit->GetTime();

            if (abs(hit->GetTime() - t0) < 20)
            {
                plane = hit->GetDetId();
                wire = hit->GetWireId();
                x[plane - 1] = hit->GetX();
                y[plane - 1] = hit->GetY();
                eloss[plane - 1] = hit->GetEloss();
                t[plane - 1] = hit->GetTime();
                // cout << "Hit plane: " << plane << " x: " << x[plane] << " y: "
                //     << y[plane] << " t: " << t[plane] << " ToT: " << eloss[plane] << endl;
                if (plane & 1)
                {
                    fh_Pdc_Tot[plane - 1]->Fill(x[plane - 1], eloss[plane - 1]);
                }
                else
                {
                    fh_Pdc_Tot[plane - 1]->Fill(y[plane - 1], eloss[plane - 1]);
                }
                fh_Pdc_Time[plane - 1]->Fill(wire, t[plane - 1]);
            }
            else
            {
                // Draw histograms
                /*
                                for(Int_t i = 0 ; i < 5 ; i++){
                                    //fh_Pdc_ToF[i]->Fill();
                                    if(i == 0) fh_Pdc_xy[i]->Fill(x[1],y[2]);
                                    if(i == 2) fh_Pdc_xy[i]->Fill(x[3],y[4]);
                                    fh_Pdc_Tot[i]->Fill(x[i+1],eloss[i+1]);
                                    //fh_Pdc_Wire_vs_Events[i]->Write();
                                }
                */
                t0 = hit->GetTime();
                ihit--;
            }
        }
    }
    fNEvents += 1;
}

void R3BOnlineSpectraPdc::FinishEvent()
{
    fMappedItems->Clear();
    fCalItems->Clear();
    fHitItems->Clear();
}

void R3BOnlineSpectraPdc::FinishTask()
{
    if (fMappedItems)
    {
        for (Int_t i = 0; i < 4; i++)
        {
            // fh_Pdc_ToF[i]->Write();
            // fh_Pdc_xy[i]->Write();
            fh_Pdc_Wire[i]->Write();
            fh_Pdc_Tot[i]->Write();
            fh_Pdc_Time[i]->Write();
            // fh_Pdc_Wire_vs_Events[i]->Write();
        }
    }
    puts("4");
}

ClassImp(R3BOnlineSpectraPdc)
