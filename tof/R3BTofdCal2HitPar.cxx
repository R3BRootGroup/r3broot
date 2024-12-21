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

// ----------------------------------------------------------------
// -----            R3BTofdCal2HitPar             -----
// -----           Created May 2016 by M.Heil             -----
// ----------------------------------------------------------------

/* Some notes:
 *
 *
 */

#include "R3BTofdCal2HitPar.h"
#include "R3BEventHeader.h"
#include "R3BLosCalData.h"
#include "R3BLosMappedData.h"
#include "R3BTofdCalData.h"
#include "R3BTofdHitPar.h"

#include "FairLogger.h"
#include "FairRootManager.h"
#include "FairRtdbRun.h"
#include "FairRunIdGenerator.h"
#include "FairRuntimeDb.h"

#include "TCanvas.h"
#include "TClonesArray.h"
#include "TF1.h"
#include "TFile.h"
#include "TGraph.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TMath.h"
#include "TProfile.h"
#include "TSpectrum.h"

#include <iostream>
#include <stdlib.h>

#define IS_NAN(x) TMath::IsNaN(x)
using namespace std;

namespace
{
    double c_range_ns = 2048 * 5;
    double c_bar_coincidence_ns = 20; // nanoseconds.
} // namespace

R3BTofdCal2HitPar::R3BTofdCal2HitPar()
    : FairTask("R3BTofdCal2HitPar", 1)
    , fCalItemsLos(NULL)
    , fUpdateRate(1000000)
    , fMinStats(100000)
    , fTrigger(-1)
    , fNofPlanes(5)
    , fPaddlesPerPlane(6)
    , fNEvents(0)
    , fCal_Par(NULL)
    , fClockFreq(1. / VFTX_CLOCK_MHZ * 1000.)
    , fTofdY(0.)
    , fTofdQ(0.)
    , fParaFile("")
{
    for (Int_t i = 0; i < N_TOFD_HIT_PLANE_MAX; i++)
    {
        fhTdiff[i] = NULL;
        fhTsync[i] = NULL;
        for (Int_t j = 0; j < N_TOFD_HIT_PADDLE_MAX; j++)
        {
            fhTotPm1[i][j] = NULL;
            fhTotPm2[i][j] = NULL;
            fhTot1vsTot2[i][j] = NULL;
            fhTot1vsPos[i][j] = NULL;
            fhTot2vsPos[i][j] = NULL;
        }
    }
}

R3BTofdCal2HitPar::R3BTofdCal2HitPar(const char* name, Int_t iVerbose)
    : FairTask(name, iVerbose)
    , fCalItemsLos(NULL)
    , fUpdateRate(1000000)
    , fMinStats(100000)
    , fTrigger(-1)
    , fNofPlanes(5)
    , fPaddlesPerPlane(6)
    , fNEvents(0)
    , fCal_Par(NULL)
    , fClockFreq(1. / VFTX_CLOCK_MHZ * 1000.)
    , fTofdY(0.)
    , fTofdQ(0.)
    , fParaFile("")
{
    for (Int_t i = 0; i < N_TOFD_HIT_PLANE_MAX; i++)
    {
        fhTdiff[i] = NULL;
        fhTsync[i] = NULL;
        for (Int_t j = 0; j < N_TOFD_HIT_PADDLE_MAX; j++)
        {
            fhTotPm1[i][j] = NULL;
            fhTotPm2[i][j] = NULL;
            fhTot1vsTot2[i][j] = NULL;
            fhTot1vsPos[i][j] = NULL;
            fhTot2vsPos[i][j] = NULL;
        }
    }
}

R3BTofdCal2HitPar::~R3BTofdCal2HitPar()
{
    for (Int_t i = 0; i < N_TOFD_HIT_PLANE_MAX; i++)
    {
        if (fhTdiff[i])
            delete fhTdiff[i];
        if (fhTsync[i])
            delete fhTsync[i];
        for (Int_t j = 0; j < N_TOFD_HIT_PADDLE_MAX; j++)
        {
            if (fhTotPm1[i][j])
                delete fhTotPm1[i][j];
            if (fhTotPm2[i][j])
                delete fhTotPm2[i][j];
            if (fhTot1vsTot2[i][j])
                delete fhTot1vsTot2[i][j];
            if (fhTot1vsPos[i][j])
                delete fhTot1vsPos[i][j];
            if (fhTot2vsPos[i][j])
                delete fhTot2vsPos[i][j];
        }
    }
    if (fCal_Par)
    {
        delete fCal_Par;
    }
}

InitStatus R3BTofdCal2HitPar::Init()
{
    FairRootManager* rm = FairRootManager::Instance();
    if (!rm)
    {
        return kFATAL;
    }

    header = dynamic_cast<R3BEventHeader*>(rm->GetObject("R3BEventHeader"));
    // may be = NULL!

    fCalData = dynamic_cast<TClonesArray*>(rm->GetObject("TofdCal"));
    if (!fCalData)
    {
        return kFATAL;
    }

    if (!fNofModules)
    {
        LOG(error) << "R3BTofdCal2HitPar::Init() Number of modules not set. ";
        return kFATAL;
    }
    fCalItemsLos = dynamic_cast<TClonesArray*>(rm->GetObject("LosCal"));
    if (NULL == fCalItemsLos)
        LOG(fatal) << "Branch LosCal not found";

    return kSUCCESS;
}

void R3BTofdCal2HitPar::SetParContainers()
{
    // container needs to be created in tcal/R3BTCalContFact.cxx AND R3BTCal needs
    // to be set as dependency in CMakelists.txt (in this case in the tof directory)
    fCal_Par = dynamic_cast<R3BTofdHitPar*>(FairRuntimeDb::instance()->getContainer("TofdHitPar"));
    if (!fCal_Par)
    {
        LOG(error) << "R3BTofdCal2HitPar::Init() Couldn't get handle on TofdHitPar. ";
    }
    //	    fCal_Par->setChanged();
}

void R3BTofdCal2HitPar::Exec(Option_t* option)
{
    if (fNEvents / 10000. == (int)fNEvents / 10000)
        cout << "Events: " << fNEvents << "         \r" << std::flush;
    // test for requested trigger (if possible)
    if ((fTrigger >= 0) && (header) && (header->GetTrigger() != fTrigger))
        return;

    Double_t timeLos = 0;
    Double_t time_r_V = nan(""), time_t_V = nan(""), time_l_V = nan(""), time_b_V = nan(""), time_rt_V = nan(""),
             time_lt_V = nan(""), time_lb_V = nan(""), time_rb_V = nan("");

    // Los detector
    if (fCalItemsLos)
    {
        Int_t nHits = fCalItemsLos->GetEntriesFast();

        for (Int_t ihit = 0; ihit < nHits; ihit++)
        {
            R3BLosCalData* calData = dynamic_cast<R3BLosCalData*>(fCalItemsLos->At(ihit));

            Int_t iDet = calData->GetDetector();
            // Int_t iCha=calData->GetChannel();

            if (!(IS_NAN(calData->GetTimeV_ns(5))))
                time_r_V = calData->GetTimeV_ns(5);
            if (!(IS_NAN(calData->GetTimeV_ns(7))))
                time_t_V = calData->GetTimeV_ns(7);
            if (!(IS_NAN(calData->GetTimeV_ns(1))))
                time_l_V = calData->GetTimeV_ns(1);
            if (!(IS_NAN(calData->GetTimeV_ns(3))))
                time_b_V = calData->GetTimeV_ns(3);
            if (!(IS_NAN(calData->GetTimeV_ns(6))))
                time_rt_V = calData->GetTimeV_ns(6);
            if (!(IS_NAN(calData->GetTimeV_ns(0))))
                time_lt_V = calData->GetTimeV_ns(0);
            if (!(IS_NAN(calData->GetTimeV_ns(2))))
                time_lb_V = calData->GetTimeV_ns(2);
            if (!(IS_NAN(calData->GetTimeV_ns(4))))
                time_rb_V = calData->GetTimeV_ns(4);
        }
        timeLos = (time_r_V + time_t_V + time_l_V + time_b_V + time_rt_V + time_lt_V + time_lb_V + time_rb_V) / 8.;
    }

    // ToFD detector

    Int_t nHits = fCalData->GetEntries();

    // Organize cals into bars.
    struct Entry
    {
        std::vector<R3BTofdCalData*> top;
        std::vector<R3BTofdCalData*> bot;
    };
    std::map<size_t, Entry> bar_map;

    for (Int_t ihit = 0; ihit < nHits; ihit++)
    {
        auto* hit = dynamic_cast<R3BTofdCalData*>(fCalData->At(ihit));
        size_t idx = hit->GetDetectorId() * fPaddlesPerPlane * hit->GetBarId();

        // std::cout << "Hits: " << hit->GetDetectorId() << ' ' << hit->GetBarId() << ' ' << hit->GetSideId() << ' '
        //          << hit->GetTimeLeading_ns() << ' ' << hit->GetTimeTrailing_ns() << '\n';

        auto ret = bar_map.insert(std::pair<size_t, Entry>(idx, Entry()));
        auto& vec = 1 == hit->GetSideId() ? ret.first->second.top : ret.first->second.bot;
        vec.push_back(hit);
    }

    // Find coincident PMT hits.
    // std::cout << "Print:\n";
    for (auto it = bar_map.begin(); bar_map.end() != it; ++it)
    {
        // for (auto it2 = it->second.top.begin(); it->second.top.end() != it2; ++it2) {
        // std::cout << "Top: " << (*it2)->GetDetectorId() << ' ' << (*it2)->GetBarId() << ' ' <<
        // (*it2)->GetTimeLeading_ns() << '\n';
        //}
        // for (auto it2 = it->second.bot.begin(); it->second.bot.end() != it2; ++it2) {
        // std::cout << "Bot: " << (*it2)->GetDetectorId() << ' ' << (*it2)->GetBarId() << ' ' <<
        // (*it2)->GetTimeLeading_ns() << '\n';
        //}
        auto const& top_vec = it->second.top;
        auto const& bot_vec = it->second.bot;
        size_t top_i = 0;
        size_t bot_i = 0;
        for (; top_i < top_vec.size() && bot_i < bot_vec.size();)
        {
            auto top = top_vec.at(top_i);
            auto bot = bot_vec.at(bot_i);
            auto top_ns = top->GetTimeLeading_ns();
            auto bot_ns = bot->GetTimeLeading_ns();

            auto dt = top_ns - bot_ns;
            // Handle wrap-around.
            auto dt_mod = fmod(dt + c_range_ns, c_range_ns);
            if (dt < 0)
            {
                // We're only interested in the short time-differences, so we
                // want to move the upper part of the coarse counter range close
                // to the lower range, i.e. we cut the middle of the range and
                // glue zero and the largest values together.
                dt_mod -= c_range_ns;
            }
            // std::cout << top_i << ' ' << bot_i << ": " << top_ns << ' ' << bot_ns << " = " << dt << ' ' <<
            // std::abs(dt_mod) << '\n';
            if (std::abs(dt_mod) < c_bar_coincidence_ns)
            {
                // Hit!

                // std::cout << "Hit!\n";
                Int_t iPlane = top->GetDetectorId(); // 1..n
                Int_t iBar = top->GetBarId();        // 1..n
                if (iPlane > fNofPlanes)             // this also errors for iDetector==0
                {
                    LOG(error) << "R3BTofdCal2HitPar::Exec() : more detectors than expected! Det: " << iPlane
                               << " allowed are 1.." << fNofPlanes;
                    continue;
                }
                if (iBar > fPaddlesPerPlane) // same here
                {
                    LOG(error) << "R3BTofdCal2HitPar::Exec() : more bars then expected! Det: " << iBar
                               << " allowed are 1.." << fPaddlesPerPlane;
                    continue;
                }

                auto top_tot = fmod(top->GetTimeTrailing_ns() - top->GetTimeLeading_ns() + c_range_ns, c_range_ns);
                auto bot_tot = fmod(bot->GetTimeTrailing_ns() - bot->GetTimeLeading_ns() + c_range_ns, c_range_ns);

                // walk corrections
                bot_ns = bot_ns - walk(bot_tot);
                top_ns = top_ns - walk(top_tot);

                auto tdiff = bot_ns - top_ns;

                // std::cout << "Saved: " << top->GetDetectorId() << ' ' << top->GetBarId() << ' '
                //<< bot->GetDetectorId() << ' ' << bot->GetBarId() << ' '
                //<< (top_ns + bot_ns)/2 << ' ' << tdiff << ' ' << top_tot << ' ' << bot_tot << '\n';

                // create histograms
                CreateHistograms(iPlane, iBar);

                // fill histograms

                fhTotPm1[iPlane - 1][iBar - 1]->Fill(bot_tot);
                fhTotPm2[iPlane - 1][iBar - 1]->Fill(top_tot);
                fhTot1vsTot2[iPlane - 1][iBar - 1]->Fill(top_tot, bot_tot);

                // Time differences of one paddle **************************
                fhTdiff[iPlane - 1]->Fill(iBar, tdiff);

                // ToF
                auto ToF = (top_ns + bot_ns) / 2 - timeLos;
                while (ToF < -c_range_ns / 2)
                    ToF += c_range_ns;
                while (ToF > c_range_ns / 2)
                    ToF -= c_range_ns;

                fhTsync[iPlane - 1]->Fill(iBar, ToF);

                // cout << "test " << (top_ns + bot_ns)/2 <<"  " << timeLos << "  " << ToF << endl;

                if (fTofdQ != 0)
                {
                    R3BTofdHitModulePar* par = fCal_Par->GetModuleParAt(iPlane, iBar);
                    if (!par)
                    {
                        LOG(info) << "R3BTofdCal2Hit::Exec : Hit par not found, Plane: " << top->GetDetectorId()
                                  << ", Bar: " << top->GetBarId();
                        continue;
                    }
                    Double_t Pos = ((bot_ns + par->GetOffset1()) - (top_ns + par->GetOffset2())) * par->GetVeff();

                    fhTot1vsPos[iPlane - 1][iBar - 1]->Fill(Pos, bot_tot);
                    fhTot2vsPos[iPlane - 1][iBar - 1]->Fill(Pos, top_tot);
                }
                ++top_i;
                ++bot_i;

                // Increment events
                fNEvents += 1;
            }
            else if (dt < 0 && dt > -c_range_ns / 2)
            {
                ++top_i;
            }
            else
            {
                ++bot_i;
            }
        }
    }
}

void R3BTofdCal2HitPar::CreateHistograms(Int_t iPlane, Int_t iBar)
{

    if (NULL == fhTdiff[iPlane - 1])
    {
        char strName1[255];
        char strName2[255];
        sprintf(strName1, "Time_Diff_Plane_%d", iPlane);
        sprintf(strName2, "Time Diff Plane %d", iPlane);
        fhTdiff[iPlane - 1] = new TH2F(strName1, strName2, 50, 0, 50, 4000, -20., 20.);
        fhTdiff[iPlane - 1]->GetXaxis()->SetTitle("Bar #");
        fhTdiff[iPlane - 1]->GetYaxis()->SetTitle("Time difference (PM1 - PM2) in ns");
    }
    if (NULL == fhTsync[iPlane - 1])
    {
        char strName[255];
        sprintf(strName, "Time_Sync_Plane_%d", iPlane);
        fhTsync[iPlane - 1] = new TH2F(strName, "", 50, 0, 50, 10000, -10, 90.);
        fhTsync[iPlane - 1]->GetXaxis()->SetTitle("Bar #");
        fhTsync[iPlane - 1]->GetYaxis()->SetTitle("ToF in ns");
    }

    if (NULL == fhTotPm1[iPlane - 1][iBar - 1])
    {
        char strName[255];
        sprintf(strName, "ToT_Plane_%d_Bar_%d_PM_1", iPlane, iBar);
        fhTotPm1[iPlane - 1][iBar - 1] = new TH1F(strName, "", 300, 0., 300.);
        fhTotPm1[iPlane - 1][iBar - 1]->GetXaxis()->SetTitle("ToT of PM1 in ns");
        fhTotPm1[iPlane - 1][iBar - 1]->GetYaxis()->SetTitle("Counts");
    }
    if (NULL == fhTotPm2[iPlane - 1][iBar - 1])
    {
        char strName[255];
        sprintf(strName, "ToT_Plane_%d_Bar_%d_PM_2", iPlane, iBar);
        fhTotPm2[iPlane - 1][iBar - 1] = new TH1F(strName, "", 300, 0., 300.);
        fhTotPm2[iPlane - 1][iBar - 1]->GetXaxis()->SetTitle("ToT of PM2 in ns");
        fhTotPm2[iPlane - 1][iBar - 1]->GetYaxis()->SetTitle("Counts");
    }
    if (NULL == fhTot1vsTot2[iPlane - 1][iBar - 1])
    {
        char strName[255];
        sprintf(strName, "Plane_%d_Bar_%d_ToT1vsToT2", iPlane, iBar);
        fhTot1vsTot2[iPlane - 1][iBar - 1] = new TH2F(strName, "", 300, 0., 300., 300, 0., 300.);
        fhTot1vsTot2[iPlane - 1][iBar - 1]->GetXaxis()->SetTitle("ToT of PM2 in ns");
        fhTot1vsTot2[iPlane - 1][iBar - 1]->GetYaxis()->SetTitle("Tot of PM1 in ns");
    }
    if (NULL == fhTot1vsPos[iPlane - 1][iBar - 1])
    {
        char strName[255];
        sprintf(strName, "Tot1_vs_Pos_Plane_%d_Bar_%d", iPlane, iBar);
        if (iPlane < 3)
            fhTot1vsPos[iPlane - 1][iBar - 1] = new TH2F(strName, "", 120, -60, 60, 400, 0., 200.);
        if (iPlane > 2)
            fhTot1vsPos[iPlane - 1][iBar - 1] = new TH2F(strName, "", 120, -60, 60, 400, 0., 200.);
        fhTot1vsPos[iPlane - 1][iBar - 1]->GetXaxis()->SetTitle("ToT of PM2 in ns");
        fhTot1vsPos[iPlane - 1][iBar - 1]->GetYaxis()->SetTitle("Tot of PM1 in ns");
    }
    if (NULL == fhTot2vsPos[iPlane - 1][iBar - 1])
    {
        char strName[255];
        sprintf(strName, "Tot2_vs_Pos_Plane_%d_Bar_%d", iPlane, iBar);
        if (iPlane < 3)
            fhTot2vsPos[iPlane - 1][iBar - 1] = new TH2F(strName, "", 120, -60, 60, 150, 50., 200.);
        if (iPlane > 2)
            fhTot2vsPos[iPlane - 1][iBar - 1] = new TH2F(strName, "", 120, -60, 60, 150, 50., 200.);
    }
}

void R3BTofdCal2HitPar::FinishEvent()
{
    if (fCalItemsLos)
    {
        fCalItemsLos->Clear();
    }
}

void R3BTofdCal2HitPar::FinishTask()
{
    for (Int_t i = 0; i < N_TOFD_HIT_PLANE_MAX; i++)
    {
        if (fhTsync[i])
            fhTsync[i]->Write();
        if (fhTdiff[i])
            fhTdiff[i]->Write();
        for (Int_t j = 0; j < N_TOFD_HIT_PADDLE_MAX; j++)
        {
            if (fhTot1vsPos[i][j])
                fhTot1vsPos[i][j]->Write();
            if (fhTot2vsPos[i][j])
                fhTot2vsPos[i][j]->Write();
            if (fhTotPm1[i][j])
                fhTotPm1[i][j]->Write();
            if (fhTotPm2[i][j])
                fhTotPm2[i][j]->Write();
            if (fhTot1vsTot2[i][j])
                fhTot1vsTot2[i][j]->Write();
            if (fhTot1vsTot2[i][j])
                fhTot1vsTot2[i][j]->Write();
        }
    }

    // Determine time offset of the 2 PMTs of one paddle. This procedure
    // assumes a sweep run in the middle of the ToF wall horizontally.
    // Since all paddles are mounted vertically one can determine the offset.
    // Half of the offset is added to PM1 and half to PM2.

    if (fTofdY == 0 && fTofdQ == 0)
    {
        calcOffset();
        calcSync();
    }
    else if (fTofdY != 0 && fTofdQ == 0)
    {
        calcVeff();
    }

    Double_t para[4];
    Double_t min = -40;
    Double_t max = 40;
    if (fTofdQ != 0)
    {
        if (fParaFile != "")
        {
            //		  char filename1[100];
            //	      sprintf(filename1,"para_paddles_run137.dat");
            ifstream a_file(fParaFile);
            if (!a_file)
            {
                cout << "Could not open file!!!" << endl;
                exit(1);
            }
            Int_t plane, bar, pm, p;

            while (!a_file.eof())
            {
                for (Int_t j = 0; j <= 3; j++)
                {
                    a_file >> plane >> bar >> pm >> p >> para[j];
                    cout << "Reading Parameter, plane " << plane << " bar " << bar << " PM " << pm << " #" << p << " : "
                         << para[j] << endl;
                }
                R3BTofdHitModulePar* par = fCal_Par->GetModuleParAt(plane, bar);
                Double_t offset1 = par->GetOffset1();
                Double_t offset2 = par->GetOffset2();
                Double_t veff = par->GetVeff();
                Double_t sync = par->GetSync();
                if (pm == 1)
                {
                    par->SetPar1a(para[0]);
                    par->SetPar1b(para[1]);
                    par->SetPar1c(para[2]);
                    par->SetPar1d(para[3]);
                }
                else
                {
                    par->SetPar2a(para[0]);
                    par->SetPar2b(para[1]);
                    par->SetPar2c(para[2]);
                    par->SetPar2d(para[3]);
                }
            }
            a_file.close();
        }
        else
        {
            for (Int_t i = 0; i < N_TOFD_HIT_PLANE_MAX; i++)
            {
                for (Int_t j = 0; j < N_TOFD_HIT_PADDLE_MAX; j++)
                {
                    if (fhTot1vsPos[i][j])
                    {
                        R3BTofdHitModulePar* par = fCal_Par->GetModuleParAt(i + 1, j + 1);

                        doubleExp(fhTot1vsPos[i][j], min, max, para);
                        Double_t offset1 = par->GetOffset1();
                        Double_t offset2 = par->GetOffset2();
                        Double_t veff = par->GetVeff();
                        Double_t sync = par->GetSync();
                        par->SetPar1a(para[0]);
                        par->SetPar1b(para[1]);
                        par->SetPar1c(para[2]);
                        par->SetPar1d(para[3]);
                    }

                    if (fhTot2vsPos[i][j])
                    {
                        R3BTofdHitModulePar* par = fCal_Par->GetModuleParAt(i + 1, j + 1);
                        doubleExp(fhTot2vsPos[i][j], min, max, para);
                        Double_t offset1 = par->GetOffset1();
                        Double_t offset2 = par->GetOffset2();
                        Double_t veff = par->GetVeff();
                        Double_t sync = par->GetSync();
                        par->SetPar2a(para[0]);
                        par->SetPar2b(para[1]);
                        par->SetPar2c(para[2]);
                        par->SetPar2d(para[3]);
                    }
                }
            }
        }
        fCal_Par->setChanged();
    }
}

void R3BTofdCal2HitPar::calcOffset()
{

    /*
       const Int_t npeaks = 10;
       Double_t sigma=1.;
       Double_t sens=2.;
       Double_t xp;
       Int_t nfound;


        for (Int_t i = 0; i < N_TOFD_HIT_PLANE_MAX; i++)
        {
            for (Int_t j = 0; j < N_TOFD_HIT_PADDLE_MAX; j++)
            {
                if (fhTdiff[i][j]) {
                    LOG(warn) << " Plane "<< i+1<<" Bar "<< j+1;
                    nfound=0;
                    TH1F *h = (TH1F*) fhTdiff[i][j]->Clone("h");
                    h->Draw();
                    TSpectrum *s = new TSpectrum(2*npeaks);
                    LOG(debug) << " Search for peaks ";
                    nfound = s->Search(h,sigma,"", sens);
                    LOG(debug) << " Found "<< nfound<<" peaks.";

                    Double_t *xpeaks = s->GetPositionX();
                    for (Int_t p=0;p<nfound;p++) {
                        xp = xpeaks[p];
                        LOG(warn) << " Position "<< xp;
                    }
                 }
             }
        }
      */

    R3BTofdHitModulePar* mpar;
    for (Int_t i = 0; i < N_TOFD_HIT_PLANE_MAX; i++)
    {
        for (Int_t j = 0; j < N_TOFD_HIT_PADDLE_MAX; j++)
        {
            mpar = new R3BTofdHitModulePar();
            Double_t offset = 0.;
            if (fhTdiff[i])
            {
                auto* histo_py = fhTdiff[i]->ProjectionY("histo_py", j + 2, j + 2, "");

                Double_t veff = mpar->GetVeff();
                Int_t binmax = histo_py->GetMaximumBin();
                offset = histo_py->GetXaxis()->GetBinCenter(binmax);
                LOG(warn) << " Plane  " << i + 1 << " Bar " << j + 1 << " Offset  " << offset;
                mpar->SetPlane(i + 1);
                mpar->SetPaddle(j + 1);
                mpar->SetOffset1(-offset / 2.);
                mpar->SetOffset2(offset / 2.);
                mpar->SetVeff(veff);
                fCal_Par->AddModulePar(mpar);
            }
        }
    }
    fCal_Par->setChanged();
}

void R3BTofdCal2HitPar::calcSync()
{

    for (Int_t i = 0; i < N_TOFD_HIT_PLANE_MAX; i++)
    {
        for (Int_t j = 0; j < N_TOFD_HIT_PADDLE_MAX; j++)
        {
            if (fhTsync[i])
            {

                auto* histo_py = fhTsync[i]->ProjectionY("histo_py", j + 2, j + 2, "");

                R3BTofdHitModulePar* par = fCal_Par->GetModuleParAt(i + 1, j + 1);
                Double_t offset1 = par->GetOffset1();
                Double_t offset2 = par->GetOffset2();
                Int_t binmax = histo_py->GetMaximumBin();
                Double_t sync = histo_py->GetXaxis()->GetBinCenter(binmax);
                par->SetOffset1(offset1);
                par->SetOffset2(offset2);
                par->SetSync(sync);
                LOG(warn) << " Plane  " << i + 1 << " Bar " << j + 1 << " Sync  " << sync;
            }
        }
    }

    fCal_Par->setChanged();
}

void R3BTofdCal2HitPar::calcVeff()
{
    for (Int_t i = 0; i < N_TOFD_HIT_PLANE_MAX; i++)
    {
        for (Int_t j = 0; j < N_TOFD_HIT_PADDLE_MAX; j++)
        {
            Double_t max = 0.;
            Double_t veff = 7.;
            if (fhTdiff[i])
            {
                R3BTofdHitModulePar* par = fCal_Par->GetModuleParAt(i + 1, j + 1);
                if (!par)
                {
                    LOG(info) << "R3BTofdCal2Hit::Exec : Hit par not found, Plane: " << i + 1 << ", Bar: " << j + 1;
                    continue;
                }

                auto* histo_py = fhTdiff[i]->ProjectionY("histo_py", j + 2, j + 2, "");

                Int_t binmax = histo_py->GetMaximumBin();
                max = histo_py->GetXaxis()->GetBinCenter(binmax);
                Double_t offset1 = par->GetOffset1();
                Double_t offset2 = par->GetOffset2();
                Double_t sync = par->GetSync();

                max = max + offset1 - offset2;
                veff = fTofdY / max; // effective speed of light in [cm/s]
                LOG(warn) << " Plane  " << i + 1 << " Bar " << j + 1 << " offset  " << par->GetOffset1();
                LOG(warn) << " Plane  " << i + 1 << " Bar " << j + 1 << " max  " << max;
                LOG(warn) << " Plane  " << i + 1 << " Bar " << j + 1 << " veff  " << veff;
                par->SetVeff(veff);
            }
        }
    }
    fCal_Par->setChanged();
}

void R3BTofdCal2HitPar::doubleExp(TH2F* histo, Double_t min, Double_t max, Double_t* para)
{
    // This fits the exponential decay of the light in a paddle. The 2 PMTs are fit with the same function but one
    // side will deliver negative attenuation parameters and the other positive.
    Double_t y[1000], x[1000];
    Int_t n = 0;
    for (Int_t j = 0; j <= 3; j++)
    {
        para[j] = 0;
    }
    TGraph* gr1 = new TGraph();
    TGraph* gr2 = new TGraph();
    TCanvas* cfit_exp = new TCanvas("cfit_exp", "fit exponential", 100, 100, 500, 500);

    cfit_exp->Clear();
    cfit_exp->Divide(1, 3);
    cfit_exp->cd(1);
    TH2F* histo1 = dynamic_cast<TH2F*>(histo->Clone());
    TH2F* histo2 = dynamic_cast<TH2F*>(histo->Clone());

    histo1->Draw("colz");

    cfit_exp->cd(2);

    for (Int_t i = 1; i < histo1->GetNbinsX() - 1; i++)
    {
        auto* histo_py = histo1->ProjectionY("histo_py", i, i, "");
        histo_py->Draw();
        x[n] = histo1->GetXaxis()->GetBinCenter(i);
        Int_t binmax = histo_py->GetMaximumBin();
        y[n] = histo_py->GetXaxis()->GetBinCenter(binmax);
        if (histo_py->GetMaximum() > 1)
            n++;
        delete histo_py;
    }

    gr1 = new TGraph(n, x, y);
    gr1->Draw("A*");

    TF1* f1 = new TF1("f1", "[0]*(exp(-[1]*(x+50))+exp(-[2]*(x+50)))+[3]", min, max);
    f1->SetParameters(110., 0.01, 1.1, 40.);
    f1->SetLineColor(2);
    gr1->Fit("f1", "", "", min, max);
    for (Int_t j = 0; j <= 3; j++)
    {
        para[j] = f1->GetParameter(j);
        cout << "Parameter: " << para[j] << endl;
    }

    // fit again but with more information and better cuts
    n = 0;
    cfit_exp->cd(3);

    for (Int_t i = 1; i < histo2->GetNbinsX(); i++)
    {
        Double_t pos = histo2->GetXaxis()->GetBinCenter(i);
        Double_t ymean = para[0] * (exp(-para[1] * (pos + 50)) + exp(-para[2] * (pos + 50))) + para[3];
        histo2->SetAxisRange(ymean - 5., ymean + 5., "Y");
        histo2->Draw("colz");
        auto* histo_py = histo2->ProjectionY("histo_py", i, i, "");
        histo_py->Draw();
        x[n] = histo2->GetXaxis()->GetBinCenter(i);
        Int_t binmax = histo_py->GetMaximumBin();
        y[n] = histo_py->GetXaxis()->GetBinCenter(binmax);
        if (histo_py->GetMaximum() > 1)
            n++;
        delete histo_py;
    }

    gr2 = new TGraph(n, x, y);
    gr2->Draw("A*");

    TF1* f2 = new TF1("f2", "[0]*(exp(-[1]*(x+50))+exp(-[2]*(x+50)))+[3]", min, max);
    f2->SetParameters(para[0], para[1], para[2], para[3]);
    f2->SetLineColor(2);
    gr2->Fit("f2", "", "", min, max);
    for (Int_t j = 0; j <= 3; j++)
    {
        para[j] = f2->GetParameter(j);
        cout << "Parameter: " << para[j] << endl;
    }

    cfit_exp->Update();
    gPad->WaitPrimitive();
    delete gr1;
    delete gr2;
    delete f1;
    delete f2;
}

Double_t R3BTofdCal2HitPar::walk(Double_t Q)
{
    Double_t y = 0;
    Double_t par1, par2, par3, par4, par5;

    Int_t voltage = 443;

    if (voltage == 443)
    {
        par1 = 5.121991e+01;
        par2 = -9.9738836e-02;
        par3 = 5.817890e+01;
        par4 = 8.208630e-03;
        par5 = -2.719159e-05;
    }

    if (voltage == 444)
    {
        par1 = 1.16631e+01;
        par2 = 3.03014e-01;
        par3 = 3.25826e+02;
        par4 = -2.04208e-01;
        par5 = 3.47279e-04;
    }

    if (voltage == 500)
    {
        par1 = 1.64344e+01;
        par2 = 2.84000e-01;
        par3 = 3.47659e+02;
        par4 = -2.70050e-01;
        par5 = 3.61515e-04;
    }

    if (voltage == 600)
    {
        par1 = 1.22606e+01;
        par2 = 3.12697e-01;
        par3 = 4.40109e+02;
        par4 = -1.86328e-01;
        par5 = 1.49519e-04;
    }

    y = -30.2 + par1 * TMath::Power(Q, par2) + par3 / Q + par4 * Q + par5 * Q * Q;

    return y;
}

Double_t R3BTofdCal2HitPar::saturation(Double_t x)
{
    Double_t kor;
    Int_t voltage = 600;

    if (voltage == 600)
    {
        if (x < 173)
        {
            kor = 0.;
        }
        else if (x > 208)
        {
            kor = -1.73665e+03 + 2.82009e+01 * 208. - 1.53846e-01 * (208. * 208.) + 2.82425e-04 * (208. * 208. * 208.);
        }
        else
        {
            kor = -1.73665e+03 + 2.82009e+01 * x - 1.53846e-01 * (x * x) + 2.82425e-04 * (x * x * x);
        }
    }

    if (voltage == 500)
    {
        if (x < 95.5)
        {
            kor = 0.;
        }
        else if (x > 124)
        {
            kor = 1.08 * x - 112.44;
        }
        else
        {
            kor = 643.257 - 16.7823 * x + 0.139822 * (x * x) - 0.000362154 * (x * x * x);
        }
    }

    if (voltage == 700)
    {
        if (x < 198)
        {
            kor = 0.;
        }
        else if (x > 298)
        {
            kor = 0.21 * x - 45.54;
        }
        else
        {
            kor = -19067 + 383.93 * x - 3.05794 * (x * x) + 0.0120429 * (x * x * x) - 2.34619e-05 * (x * x * x * x) +
                  1.81076e-08 * (x * x * x * x * x);
        }
    }

    return kor;
}

ClassImp(R3BTofdCal2HitPar)
