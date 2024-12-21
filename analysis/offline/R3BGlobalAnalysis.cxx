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
// -----                  R3BGlobalAnalysis                -----
// -----          Created April 13th 2016 by M.Heil       -----
// ------------------------------------------------------------

/*
 * This task should fill histograms with detector variables which allow
 * to test the detectors online
 *
 */

#include "R3BLosCalData.h"
#include "R3BLosHitData.h"
#include "R3BLosMappedData.h"

#include "R3BGlobalAnalysis.h"

#include "R3BSci8CalData.h"
#include "R3BSci8MappedData.h"

#include "R3BTofdCalData.h"
#include "R3BTofdHitData.h"
#include "R3BTofdMappedData.h"

#include "R3BPaddleCalData.h"

#include "R3BPspxCalData.h"
#include "R3BPspxMappedData.h"

#include "R3BEventHeader.h"
#include "R3BTCalEngine.h"

#include "R3BBunchedFiberCalData.h"
#include "R3BBunchedFiberHitData.h"
#include "R3BBunchedFiberMappedData.h"

#include "FairLogger.h"
#include "FairRootManager.h"
#include "FairRunAna.h"
#include "FairRunOnline.h"
#include "FairRuntimeDb.h"

#include "TCanvas.h"
#include "TH1F.h"
#include "TH2F.h"

#include "TClonesArray.h"
#include "TMath.h"
#include <TRandom3.h>
#include <TRandomGen.h>
#include <array>
#include <cstdlib>
#include <ctime>
#include <fstream>
#include <iostream>
#include <sstream>
#define IS_NAN(x) TMath::IsNaN(x)
using namespace std;

R3BGlobalAnalysis::R3BGlobalAnalysis()
    : R3BGlobalAnalysis("GlobalAnalysis", 1)
{
}

R3BGlobalAnalysis::R3BGlobalAnalysis(const char* name, Int_t iVerbose)
    : FairTask(name, iVerbose)
    , fTrigger(-1)
    , fNofPlanes(N_PLANE_MAX_TOFD)
    , fPaddlesPerPlane(N_PADDLE_MAX_TOFD)
    , fClockFreq(1. / VFTX_CLOCK_MHZ * 1000.)
    , fNEvents(0)
    , fFibEvents(0)
{
}

R3BGlobalAnalysis::~R3BGlobalAnalysis()
{
    for (int i = 0; i < NOF_FIB_DET; i++)
    {
        delete fh_channels_Fib[i];
        delete fh_fibers_Fib[i];
        delete fh_fiber_Fib[i];
        delete fh_mult_Fib[i];
        delete fh_time_Fib[i];
        delete fh_xpos_Fib[i];
        delete fh_ypos_Fib[i];
        delete fh_multihit_m_Fib[i];
        delete fh_multihit_s_Fib[i];
        delete fh_ToT_m_Fib[i];
        delete fh_ToT_s_Fib[i];
    }
}

InitStatus R3BGlobalAnalysis::Init()
{

    // Initialize random number:
    std::srand(std::time(0)); // use current time as seed for random generator

    LOG(info) << "R3BGlobalAnalysis::Init ";

    // try to get a handle on the EventHeader. EventHeader may not be
    // present though and hence may be null. Take care when using.

    FairRootManager* mgr = FairRootManager::Instance();
    if (NULL == mgr)
        LOG(fatal) << "FairRootManager not found";
    header = dynamic_cast<R3BEventHeader*>(mgr->GetObject("R3BEventHeader"));
    FairRunOnline* run = FairRunOnline::Instance();

    // Get objects for detectors on all levels
    assert(DET_MAX + 1 == sizeof(fDetectorNames) / sizeof(fDetectorNames[0]));
    printf("Have %d fiber detectors.\n", NOF_FIB_DET);
    for (int det = 0; det < DET_MAX; det++)
    {
        /*
                        fMappedItems.push_back((TClonesArray *)mgr->GetObject(Form("%sMapped", fDetectorNames[det])));
                        if (NULL == fMappedItems.at(det)) {
                          printf("Could not find mapped data for '%s'.\n", fDetectorNames[det]);
                        }
                        fCalItems.push_back((TClonesArray *)mgr->GetObject(Form("%sCal", fDetectorNames[det])));
                        fHitItems.push_back((TClonesArray *)mgr->GetObject(Form("%sHit", fDetectorNames[det])));
        */
        fMappedItems.push_back(dynamic_cast<TClonesArray*>(mgr->GetObject(Form("%sMapped", fDetectorNames[det]))));
        if (NULL == fMappedItems.at(det))
        {
            printf("Could not find mapped data for '%s'.\n", fDetectorNames[det]);
        }
        fCalItems.push_back(dynamic_cast<TClonesArray*>(mgr->GetObject(Form("%sCal", fDetectorNames[det]))));
        if (NULL == fCalItems.at(det))
        {
            printf("Could not find Cal data for '%s'.\n", fDetectorNames[det]);
        }
        fHitItems.push_back(dynamic_cast<TClonesArray*>(mgr->GetObject(Form("%sHit", fDetectorNames[det]))));
        if (NULL == fHitItems.at(det))
        {
            printf("Could not find hit data for '%s'.\n", fDetectorNames[det]);
        }
    }

    //------------------------------------------------------------------------
    // create histograms of all detectors
    //------------------------------------------------------------------------

    char strNameC[255];
    sprintf(strNameC, "Cave_C_position");
    fh_Cave_position = new TH2F(strNameC, "", 2100, -100., 2000., 1000, -500., 500.);

    if (fHitItems.at(DET_LOS) || fCalItems.at(DET_LOS))
    {
        char strName[255];
        sprintf(strName, "LOS_X_vs_Y_MCFD");
        //     fh_los_pos = new TH2F(strName, "", 2000, -10., 10., 2000, -10., 10.);
        fh_los_pos = new TH2F(strName, "", 2000, -10., 10., 2000, -10., 10.);

        sprintf(strName, "LosQ_vs_TofdQ");
        fhChargeLosTofD = new TH2F(strName, "", 1000, 0., 80., 1000, 0., 80.);
        fhChargeLosTofD->GetYaxis()->SetTitle("Charge LOS");
        fhChargeLosTofD->GetXaxis()->SetTitle("Charge ToFD");
    }

    //-----------------------------------------------------------------------
    // Fiber Detectors 1-NOF_FIB_DET

    char canvName[255];

    for (Int_t ifibcount = 0; ifibcount < NOF_FIB_DET; ifibcount++)
    {

        std::stringstream fibhitobj;
        std::stringstream fibmapobj;
        std::string temp;
        std::stringstream FileNameTemp;

        std::stringstream histName1, histName2, histName22, histName3, histName4, histName5, histName6, histName33,
            histName7, histName8, histName9, histName10, histName11, histName12, histName13, histName14;
        std::stringstream histTitle1, histTitle2, histTitle22, histTitle3, histTitle4, histTitle5, histTitle6,
            histTitle33, histTitle7, histTitle8, histTitle9, histTitle10, histTitle11, histTitle12, histTitle13,
            histTitle14;

        if (fCalItems.at(DET_FI_FIRST + ifibcount) || fHitItems.at(DET_FI_FIRST + ifibcount))
        {

            const char* chistName;
            const char* chistTitle;
            const char* detName;
            const char* detName2;
            const char* FileName;
            std::string tempTitle;
            std::string tempName;
            std::stringstream tempFibName;
            std::string tempFibNames;
            std::stringstream tempCanvName;

            detName = fDetectorNames[DET_FI_FIRST + ifibcount];

            cout << "I am creating canvas " << detName << endl;

            // Channels:
            fh_channels_Fib[ifibcount] =
                new TH1F(Form("%s_channels", detName), Form("%s channels", detName), 520, 0., 520.);
            fh_channels_Fib[ifibcount]->GetXaxis()->SetTitle("Channel number");
            fh_channels_Fib[ifibcount]->GetYaxis()->SetTitle("Counts");

            // Fibers:
            fh_fibers_Fib[ifibcount] =
                new TH1F(Form("%s_fibers", detName), Form("%s fibers", detName), N_FIBER_PLOT, 0., N_FIBER_PLOT);
            fh_fibers_Fib[ifibcount]->GetXaxis()->SetTitle("Fiber number");
            fh_fibers_Fib[ifibcount]->GetYaxis()->SetTitle("Counts");

            // Fiber:
            fh_fiber_Fib[ifibcount] =
                new TH1F(Form("%s_fiber", detName), Form("%s fiber", detName), N_FIBER_PLOT, 0., N_FIBER_PLOT);
            fh_fiber_Fib[ifibcount]->GetXaxis()->SetTitle("Fiber number");
            fh_fiber_Fib[ifibcount]->GetYaxis()->SetTitle("Counts");

            // Multiplicity (number of hit fibers):
            fh_mult_Fib[ifibcount] = new TH1F(Form("%s_mult", detName), Form("%s # of fibers", detName), 100, 0., 100.);
            fh_mult_Fib[ifibcount]->GetXaxis()->SetTitle("Multiplicity");
            fh_mult_Fib[ifibcount]->GetYaxis()->SetTitle("Counts");

            // Multihit MAPMT:
            fh_multihit_m_Fib[ifibcount] = new TH2F(
                Form("%s_multihit_m", detName), Form("%s multihits MAPMT", detName), 520, 0., 520., 20, 0., 20.);
            fh_multihit_m_Fib[ifibcount]->GetXaxis()->SetTitle("MAPMT channel");
            fh_multihit_m_Fib[ifibcount]->GetYaxis()->SetTitle("Multihit");

            // Multihit SAPMT:
            fh_multihit_s_Fib[ifibcount] = new TH2F(
                Form("%s_multihit_s", detName), Form("%s multihits single PMT", detName), 16, 0., 16., 20, 0., 20.);
            fh_multihit_s_Fib[ifibcount]->GetXaxis()->SetTitle("SAPMT channel");
            fh_multihit_s_Fib[ifibcount]->GetYaxis()->SetTitle("Multihit");

            // ToT MAPMT:
            fh_ToT_m_Fib[ifibcount] = new TH2F(Form("%s_tot_m", detName),
                                               Form("%s ToT of MAPMT", detName),
                                               N_FIBER_PLOT,
                                               0.,
                                               N_FIBER_PLOT,
                                               400,
                                               0.,
                                               200.);
            fh_ToT_m_Fib[ifibcount]->GetXaxis()->SetTitle("Fiber number");
            fh_ToT_m_Fib[ifibcount]->GetYaxis()->SetTitle("ToT / ns");

            // ToT SAPMT:
            fh_ToT_s_Fib[ifibcount] = new TH2F(Form("%s_tot_s", detName),
                                               Form("%s ToT of single PMT", detName),
                                               N_FIBER_PLOT,
                                               0.,
                                               N_FIBER_PLOT,
                                               400,
                                               0.,
                                               400.);
            fh_ToT_s_Fib[ifibcount]->GetXaxis()->SetTitle("Fiber number");
            fh_ToT_s_Fib[ifibcount]->GetYaxis()->SetTitle("ToT / ns");

            // Time of fiber:
            fh_time_Fib[ifibcount] = new TH2F(Form("%s_TimevsFiber", detName),
                                              Form("%s Time vs Fiber", detName),
                                              N_FIBER_PLOT,
                                              0.,
                                              N_FIBER_PLOT,
                                              20000,
                                              -1024.,
                                              1024.);
            fh_time_Fib[ifibcount]->GetXaxis()->SetTitle("Fiber number");
            fh_time_Fib[ifibcount]->GetYaxis()->SetTitle("tMAPMT-tSPMT");

            // ToF Tofd -> Fiber:
            fh_Fib_ToF[ifibcount] = new TH2F(Form("%s_tof", detName),
                                             Form("%s ToF Tofd to Fiber", detName),
                                             N_FIBER_PLOT,
                                             0.,
                                             N_FIBER_PLOT,
                                             10000,
                                             -5000.,
                                             5000.);
            fh_Fib_ToF[ifibcount]->GetYaxis()->SetTitle("ToF / ns");
            fh_Fib_ToF[ifibcount]->GetXaxis()->SetTitle("Fiber number");

            // ToF Tofd -> Fiber vs. event number:
            fh_ToF_vs_events[ifibcount] = new TH2F(Form("%s_tof_vs_events", detName),
                                                   Form("%s ToF Tofd to Fiber vs event number", detName),
                                                   1000,
                                                   0,
                                                   5e6,
                                                   2200,
                                                   -1100,
                                                   1100);
            fh_ToF_vs_events[ifibcount]->GetYaxis()->SetTitle("ToF / ns");
            fh_ToF_vs_events[ifibcount]->GetXaxis()->SetTitle("event number");

            // Not-calibrated position:
            fh_xpos_Fib[ifibcount] =
                new TH1F(Form("%s_pos", detName), Form("%s Not-calibrated position", detName), 6000, -1500., 1500.);
            fh_xpos_Fib[ifibcount]->GetXaxis()->SetTitle("x-Position");
            fh_xpos_Fib[ifibcount]->GetYaxis()->SetTitle("Counts");

            // hit fiber number vs. event number:
            fh_Fib_vs_Events[ifibcount] = new TH2F(Form("%s_fib_vs_event", detName),
                                                   Form("%s Fiber # vs. Event #", detName),
                                                   1000,
                                                   0,
                                                   5e6,
                                                   1100,
                                                   0.,
                                                   1100.);
            fh_Fib_vs_Events[ifibcount]->GetYaxis()->SetTitle("Fiber number");
            fh_Fib_vs_Events[ifibcount]->GetXaxis()->SetTitle("Event number");

            // hit MA number vs. event number:
            fh_Fibs_vs_Events[ifibcount] = new TH2F(Form("%s_fibs_vs_event", detName),
                                                    Form("%s Fiber # vs. Event #", detName),
                                                    1000,
                                                    0,
                                                    5e6,
                                                    1100,
                                                    0.,
                                                    1100.);
            fh_Fibs_vs_Events[ifibcount]->GetYaxis()->SetTitle("Fiber number");
            fh_Fibs_vs_Events[ifibcount]->GetXaxis()->SetTitle("Event number");

            // hit fiber number vs. TofD x-position:
            fh_Fibs_vs_Tofd[ifibcount] = new TH2F(Form("%s_fib_vs_TofdX", detName),
                                                  Form("%s Fiber # vs. Tofd x-pos", detName),
                                                  200,
                                                  -100,
                                                  100,
                                                  1100,
                                                  0.,
                                                  1100.);
            fh_Fibs_vs_Tofd[ifibcount]->GetYaxis()->SetTitle("Fiber number");
            fh_Fibs_vs_Tofd[ifibcount]->GetXaxis()->SetTitle("Tofd x-pos in cm");

            // hit fiber number vs. TofD y-position:
            fh_Fibs_vs_Tofd_y[ifibcount] = new TH2F(Form("%s_fib_vs_TofdY", detName),
                                                    Form("%s Fiber # vs. Tofd y-pos", detName),
                                                    800,
                                                    -100,
                                                    100,
                                                    1100,
                                                    0.,
                                                    1100.);
            fh_Fibs_vs_Tofd_y[ifibcount]->GetYaxis()->SetTitle("Fiber number");
            fh_Fibs_vs_Tofd_y[ifibcount]->GetXaxis()->SetTitle("Tofd y-pos in cm");

            // hit fiber number vs. fiber number:
            for (Int_t j = ifibcount + 1; j < NOF_FIB_DET; j++)
            {
                detName2 = fDetectorNames[DET_FI_FIRST + j];
                fh_Fib_vs_Fib[ifibcount][j] = new TH2F(Form("fib%s_vs_fib%s", detName, detName2),
                                                       Form("Fiber %s vs. Fiber %s", detName, detName2),
                                                       110,
                                                       0,
                                                       1100,
                                                       110,
                                                       0.,
                                                       1100.);
                fh_Fib_vs_Fib[ifibcount][j]->GetYaxis()->SetTitle(Form("Fiber%s", detName2));
                fh_Fib_vs_Fib[ifibcount][j]->GetXaxis()->SetTitle(Form("Fiber%s", detName));
                ;
            }

            // dx between fibers vs x
            for (Int_t j = ifibcount + 1; j < NOF_FIB_DET; j++)
            {
                detName2 = fDetectorNames[DET_FI_FIRST + j];
                fh_Fib_dx[ifibcount][j] = new TH2F(Form("fib%s_fib%s_dx", detName, detName2),
                                                   Form("dx of Fiber %s and Fiber %s", detName2, detName),
                                                   1100,
                                                   0,
                                                   1100,
                                                   1000,
                                                   -500.,
                                                   500.);
                fh_Fib_dx[ifibcount][j]->GetYaxis()->SetTitle("dx");
                fh_Fib_dx[ifibcount][j]->GetXaxis()->SetTitle(Form("Fiber%s", detName));
                ;
            }

            // Not-calibrated position:
            histName14 << detName << "_ypos";
            tempName = histName14.str();
            chistName = tempName.c_str();
            histTitle14 << detName << " Not-calibrated y-position ";
            tempTitle = histTitle14.str();
            chistTitle = tempTitle.c_str();
            fh_ypos_Fib[ifibcount] = new TH1F(chistName, chistTitle, 6000, -1500., 1500.);
            fh_ypos_Fib[ifibcount]->GetXaxis()->SetTitle("y-Position");
            fh_ypos_Fib[ifibcount]->GetYaxis()->SetTitle("Counts");
            tempName.clear();
            tempTitle.clear();

        } // end if(Mapped)

    } // end for(ifibcount)

    //---------------------------------------------------------------------------------------------------
    // TofD detector

    if (fHitItems.at(DET_TOFD) || fCalItems.at(DET_TOFD))
    {

        char strName1[255];
        char strName2[255];
        sprintf(strName1, "TofD_X_vs_Y");
        sprintf(strName2, "TofD X vs. Y");
        fh_tofd_pos = new TH2F(strName1, strName2, 20000, -100., 100., 20000, -100., 100.);

        sprintf(strName1, "TofD_Q");
        sprintf(strName2, "TofD Q");
        fh_tofd_charge = new TH1F(strName1, strName2, 200, 0., 20.);

        fh_TimePreviousEvent = new TH1F("TimePreviousEvent", "Time between 2 particles ", 300000, 0, 3000);
        fh_TimePreviousEvent->GetXaxis()->SetTitle("time / ns");
        fh_TimePreviousEvent->GetYaxis()->SetTitle("counts");

        fh_tofd_mult = new TH1F("tofd_mult", "ToFD multiplicits ", 100, 0, 100);
        fh_tofd_mult->GetXaxis()->SetTitle("multiplicity");
        fh_tofd_mult->GetYaxis()->SetTitle("counts");

        /*
                for (Int_t j = 0; j < 4; j++)
                {
                    char strName1[255];
                    sprintf(strName1, "tofd_channels_plane_%d", j+1);
                    char strName2[255];
                    sprintf(strName2, "Tofd channels plane %d", j+1);
                    fh_tofd_channels[j] = new TH1F(strName1, strName2, 90, -45., 45.);
                    fh_tofd_channels[j]->GetXaxis()->SetTitle("Channel");
                    fh_tofd_channels[j]->GetYaxis()->SetTitle("Counts");

                    char strName3[255];
                    sprintf(strName3, "tofd_ToT_plane_%d", j+1);
                    char strName4[255];
                    sprintf(strName4, "Tofd ToT plane %d", j+1);
                    fh_tofd_TotPm[j] = new TH2F(strName3, strName4, 90,-45,45,3000,0.,300.);
                    fh_tofd_TotPm[j]->GetXaxis()->SetTitle("Bar number");
                    fh_tofd_TotPm[j]->GetYaxis()->SetTitle("ToT / ns");

                    char strName7[255];
                    sprintf(strName7, "tofd_ToF_plane_%d", j+1);
                    char strName8[255];
                    sprintf(strName8, "TOFD-LOS ToF plane %d", j+1);
                    fh_tofd_ToF[j] = new TH2F(strName7, strName8, 45,0,45, 1e4, 0, 1e3);
                    fh_tofd_ToF[j]->GetXaxis()->SetTitle("Bar number");
                    fh_tofd_ToF[j]->GetYaxis()->SetTitle("ToF / ns");

                    char strName9[255];
                    sprintf(strName9, "tofd_multihit_plane_%d", j+1);
                    char strName10[255];
                    sprintf(strName10, "Tofd multihit plane %d", j+1);
                    fh_tofd_multihit[j] = new TH2F(strName9, strName10, 45, 0., 45., 10, 0, 10);
                    fh_tofd_multihit[j]->GetXaxis()->SetTitle("Bar number");
                    fh_tofd_multihit[j]->GetYaxis()->SetTitle("Multihit");



                    for(Int_t p = 0;p<N_PADDLE_MAX_PTOF; p++){
                        char strName11[255];
                        sprintf(strName11, "tofd_ToT_vs_pos_paddle_%d_plane_%d", p+1,j+1);
                        char strName12[255];
                        sprintf(strName12, "Tofd ToT vs. Pos bar %d plane %d", p+1, j+1);
                        fh_tofd_TotvsPos[j][p] = new TH2F(strName11, strName12, 100, -50., 50., 300, 0, 300);
                        fh_tofd_TotvsPos[j][p]->GetXaxis()->SetTitle("position");
                        fh_tofd_TotvsPos[j][p]->GetYaxis()->SetTitle("ToT / ns");
                    }



                    if(j<3)
                    {
                        int jk = j + 1;
                        char strName5[255];
                        sprintf(strName5, "tofd_ToF_plane_%d_%d", jk, jk+1);
                        char strName6[255];
                        sprintf(strName6, "Tofd ToF plane %d and %d ", jk, jk+1);
                        fh_tofd_dt[j] = new TH2F(strName5, strName6, 50,0,50,8000, -40., 40);
                        fh_tofd_dt[j]->GetXaxis()->SetTitle("Bar number");
                        fh_tofd_dt[j]->GetYaxis()->SetTitle("dt / ns");
                    }
                }



                char strName11[255];
                sprintf(strName11, "tofd_ToT_vs_pos_paddle_plane");
                char strName12[255];
                sprintf(strName12, "Tofd ToT vs. Pos bar plane");
                fh_tofd_TotvsPos = new TH2F(strName11, strName12, 1000, -50., 50., 3000, 0, 300);
                fh_tofd_TotvsPos->GetXaxis()->SetTitle("position");
                fh_tofd_TotvsPos->GetYaxis()->SetTitle("ToT / ns");

        */
    }

    // -------------------------------------------------------------------------

    return kSUCCESS;
}

void R3BGlobalAnalysis::Exec(Option_t* option)
{
    if (fNEvents / 1000000. == (int)fNEvents / 1000000)
        cout << "Events: " << fNEvents << flush << '\r';

    FairRootManager* mgr = FairRootManager::Instance();
    if (NULL == mgr)
        LOG(fatal) << "FairRootManager not found";
    // check for requested trigger (Todo: should be done globablly / somewhere else)
    if ((fTrigger >= 0) && (header) && (header->GetTrigger() != fTrigger))
        return;

    //----------------------------------------------------------------------
    // LOS detector
    //----------------------------------------------------------------------
    Double_t timeLos[16] = { 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0. };
    Double_t LosQ = 0.;

    if (fHitItems.at(DET_LOS))
    {
        auto det = fHitItems.at(DET_LOS);
        Int_t nHits = det->GetEntriesFast();

        for (Int_t ihit = 0; ihit < nHits; ihit++)
        {
            R3BLosHitData* hitData = dynamic_cast<R3BLosHitData*>(det->At(ihit));
            timeLos[ihit] = hitData->fTime_ns;

            fh_los_pos->Fill(hitData->fX_cm, hitData->fY_cm);

            fh_Cave_position->Fill(0., hitData->fX_cm);

            if (hitData->fZ > LosQ)
                LosQ = hitData->fZ;
        }
    }

    //----------------------------------------------------------------------
    // TOFD
    //----------------------------------------------------------------------
    Double_t timeTofd[16] = { 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0. };
    Double_t xTofd[16] = { 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0. };
    Double_t yTofd[16] = { 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0. };
    Double_t qTofd[16] = { 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0. };
    Double_t TofdQ = 0.;
    Double_t TofdX = 0.;
    Double_t TofdY = 0.;

    if (fHitItems.at(DET_TOFD))
    {
        auto det = fHitItems.at(DET_TOFD);
        Int_t nHits = det->GetEntriesFast();
        //	    cout<<"ToFD hits: "<<nHits<<endl;
        unsigned long long time0 = header->GetTimeStamp();
        double_t time1 = -1.;

        if (nHits > 0)
        {
            fh_tofd_mult->Fill(nHits);

            for (Int_t ihit = 0; ihit < nHits; ihit++)
            {
                R3BTofdHitData* hit = dynamic_cast<R3BTofdHitData*>(det->At(ihit));
                if (!hit)
                    continue; // should not happen
                if (ihit > 15)
                {
                    // cout << "Error, more than 16 hits" << endl;
                    continue;
                }
                timeTofd[ihit] = hit->GetTime();
                xTofd[ihit] = hit->GetX();
                yTofd[ihit] = hit->GetY();
                qTofd[ihit] = hit->GetEloss();
                if (qTofd[ihit] > TofdQ)
                {
                    TofdQ = qTofd[ihit];
                    TofdX = xTofd[ihit];
                    TofdY = yTofd[ihit];
                }
                //		    cout<<"ToFD: "<<ihit<<" x: "<< xTofd[ihit] << " y: " << yTofd[ihit]
                //		        << " q: "<< qTofd[ihit] << " t: "<< timeTofd[ihit] << endl;

                fh_tofd_pos->Fill(xTofd[ihit], yTofd[ihit]);
                fh_Cave_position->Fill(1000., xTofd[ihit]);
                fh_tofd_charge->Fill(qTofd[ihit]);

                double_t time2 = 0.;
                if (time1 < 0 && qTofd[ihit] > 6.)
                    time1 = timeTofd[ihit];
                if (qTofd[ihit] > 6.)
                    time2 = timeTofd[ihit];

                //			cout << "Time0: "<<time0<<" time1: "<< time1 <<" time2: "<< time2<<endl;

                if (time2 > time1 && time1 > 0. && time2 > 0.)
                {
                    //				cout<<"Time Test "<<time0<<"  "<<time1<< "   "<< time2 <<"  " <<time_previous_event
                    //<<endl;
                    fh_TimePreviousEvent->Fill(time2 - time1);
                    time1 = time2;
                    //				t0_prev=time0;
                    //				time_previous_event=double(time0)+time2-time1;

                    //				cout << "Time0: "<<time0<<" time1: "<< time1 <<" time2: "<< time2<<endl;
                    //				cout << "Result: "<<((double(time0)-time_previous_event)+time2-time1)<<endl;
                    //				cout << "time_previous_event " <<time_previous_event<<endl;
                }
            }
        }
    }

    // fhChargeLosTofD->Fill(TofdQ,LosQ/2.);

    //----------------------------------------------------------------------
    // Fiber detectors
    //----------------------------------------------------------------------
    Double_t dtime = 0.0 / 0.0;
    Int_t test = 0;
    Int_t counter = 0;
    for (Int_t ifibcount = 0; ifibcount < NOF_FIB_DET; ifibcount++)
    {
        Int_t iFib = 0;
        FibMax[ifibcount] = 0;
        auto detMapped = fMappedItems.at(DET_FI_FIRST + ifibcount);
        auto detHit = fHitItems.at(DET_FI_FIRST + ifibcount);

        if (detMapped)
        {
            Int_t nHits = detMapped->GetEntriesFast();
            std::vector<UInt_t> mapmt_num(512);
            std::vector<UInt_t> spmt_num(16);
            for (Int_t ihit = 0; ihit < nHits; ihit++)
            {
                R3BBunchedFiberMappedData* hit = dynamic_cast<R3BBunchedFiberMappedData*>(detMapped->At(ihit));
                if (!hit)
                    continue;

                // channel numbers are stored 1-based (1..n)
                Int_t iCha = hit->GetChannel(); // 1..

                if (hit->IsMAPMT() && hit->IsLeading())
                {
                    fh_channels_Fib[ifibcount]->Fill(iCha); // Fill which clockTDC channel has events
                    ++mapmt_num.at(hit->GetChannel() - 1);  // multihit of a given clockTDC channel
                }

                if (!hit->IsMAPMT() && hit->IsLeading())
                {
                    ++spmt_num.at(hit->GetChannel() - 1); // multihit of a given PADI channel
                }
            }

            for (int i = 0; i < 512; ++i)
            {
                auto m = mapmt_num.at(i);
                if (m > 0)
                    fh_multihit_m_Fib[ifibcount]->Fill(i, m); // multihit of a given clockTDC channel
            }

            for (int i = 0; i < 16; ++i)
            {
                auto s = spmt_num.at(i);

                if (s > 0)
                    fh_multihit_s_Fib[ifibcount]->Fill(i, s); // multihit of a given PADI channel
            }
        }

        if (detHit)
        {

            Int_t nHits = detHit->GetEntriesFast();
            //            cout<<"Fiber Hits: "<<nHits<<endl;

            std::vector<UInt_t> fiber_num(2049);

            Double_t xposfib = 0. / 0.;
            Double_t yposfib = 0. / 0.;
            Double_t totMax = 0.;
            Double_t tfib = 0. / 0., tof_fib_s = 0. / 0., tof_fib_m = 0. / 0.;
            Double_t randx;
            Int_t iFibMax = 0;
            Double_t totMax_MA = 0.;
            Int_t iFibMax_MA = 0;

            Double_t spmtMax = nan("");
            Double_t mapmtMax = nan("");
            Double_t tofMax = nan("");
            Double_t y;
            Double_t tof[14] = { 0., 0., 0., 0., 60., 60., 0., 0., 0., 140., 140., 0., 140., 140. };
            Double_t z[14] = { -51., -49., 0., 0., 45., 50., 0., 0., 0., 700., 670., 0., 630., 600. };
            Double_t x[14] = { 0., 0., 0., 0., 0., 0., 0., 0., 0., -200., -180., 0., -150., -130. };

            if (counter == 0 && nHits > 0)
            {
                fFibEvents++;
                counter = 1;
            }

            for (Int_t ihit = 0; ihit < nHits; ihit++)
            {
                Double_t tMAPMT = 0. / 0.;
                Double_t tSPMT = 0. / 0.;

                R3BBunchedFiberHitData* hit = dynamic_cast<R3BBunchedFiberHitData*>(detHit->At(ihit));
                if (!hit)
                    continue;

                iFib = hit->GetFiberId(); // 1..
                ++fiber_num.at(iFib);     // multihit of a given fiber
                auto m_fib = fiber_num.at(iFib);

                // times
                tMAPMT = hit->GetMAPMTTime_ns();
                tSPMT = hit->GetSPMTTime_ns();

                dtime = tMAPMT - tSPMT;

                // "Push" two times in the same clock cycle:
                while (dtime < -1024)
                {
                    tSPMT = tSPMT - 2048.;
                    dtime = tMAPMT - tSPMT;
                }
                while (dtime > 1024.)
                {
                    tSPMT = tSPMT + 2048.;
                    dtime = tMAPMT - tSPMT;
                }

                // "Push" the Fib times in the same cycle with Tofd:
                if (timeTofd[0] > 0. && !(IS_NAN(timeTofd[0])))
                {
                    while (tSPMT - timeTofd[0] < -1024.)
                    {
                        tMAPMT = tMAPMT + 2048.;
                        tSPMT = tSPMT + 2048.;
                    }
                    while (tSPMT - timeTofd[0] > 1024.)
                    {
                        tMAPMT = tMAPMT - 2048.;
                        tSPMT = tSPMT - 2048.;
                    }
                }

                // find fiber with most light
                //      Double_t ToT = sqrt(hit->GetSPMTToT_ns()*hit->GetMAPMTToT_ns());
                Double_t ToT = (3. * hit->GetMAPMTToT_ns() + hit->GetSPMTToT_ns()) / 4.;
                Double_t ToT_MA = hit->GetMAPMTToT_ns();

                // Not-calibrated ToF:
                //      tfib = (tMAPMT + tSPMT) / 2.;
                tfib = tSPMT;
                if (tfib > 0. && !(IS_NAN(tfib)) && timeTofd[0] > 0. && !(IS_NAN(timeTofd[0])))
                    tof_fib_s = tfib - timeTofd[0];
                if (tMAPMT > 0. && !(IS_NAN(tMAPMT)) && timeTofd[0] > 0. && !(IS_NAN(timeTofd[0])))
                    tof_fib_m = tMAPMT - timeTofd[0];

                // if(fNEvents<10000 && ifibcount == 13){
                if (fNEvents < 10000 && iFib > 900 && iFib < 1020 && ifibcount == 13)
                {
                    test = 0;
                    cout << "Det: " << ifibcount << " Fiber: " << iFib << " MAToT: " << hit->GetMAPMTToT_ns()
                         << " SToT: " << hit->GetSPMTToT_ns() << " sqrt: " << ToT << " tof: " << tof_fib_s
                         << " tM: " << tof_fib_m << endl;
                }

                //                if(ToT>totMax && abs(tof_fib_s-tof[ifibcount])<20. && ToT<1000.) {
                if (ToT > totMax && ToT < 1000.)
                {
                    totMax = ToT;
                    iFibMax = iFib;
                    spmtMax = hit->GetSPMTToT_ns();
                    mapmtMax = hit->GetMAPMTToT_ns();
                    // tofMax=tof_fib_m;
                    tofMax = tof_fib_s;
                    // Not-calibrated position:
                    randx = (std::rand() / (float)RAND_MAX);
                    xposfib = (-n_fiber[ifibcount] / 2. + iFib + (0.5 - randx));
                    yposfib = tMAPMT - tSPMT;
                }

                // if(ToT_MA>totMax_MA && abs(tof_fib_s-tof[ifibcount])<20. && ToT_MA<1000.) {
                if (ToT_MA > totMax_MA)
                {
                    totMax_MA = ToT_MA;
                    iFibMax_MA = iFib;
                }

            } // end for(ihit)

            //      if(abs(tof_fib_s-tof[ifibcount])<20.)
            if (totMax > 0)
            {
                if (TofdQ > 2. && TofdQ < 3.)
                {
                    fh_fibers_Fib[ifibcount]->Fill(iFibMax);
                    fh_fiber_Fib[ifibcount]->Fill(iFibMax);
                    fh_ToT_s_Fib[ifibcount]->Fill(iFibMax, spmtMax);
                    fh_ToT_m_Fib[ifibcount]->Fill(iFibMax, mapmtMax);
                    fh_time_Fib[ifibcount]->Fill(iFibMax, yposfib);
                    fh_Fib_ToF[ifibcount]->Fill(iFibMax, tofMax);
                    fh_xpos_Fib[ifibcount]->Fill(xposfib);
                    fh_ypos_Fib[ifibcount]->Fill(yposfib);
                    fh_Fibs_vs_Events[ifibcount]->Fill(fFibEvents, iFibMax);
                    fh_Fib_vs_Events[ifibcount]->Fill(fFibEvents, iFibMax_MA);
                    fh_ToF_vs_events[ifibcount]->Fill(fFibEvents, tofMax);
                    FibMax[ifibcount] = iFibMax;
                }
            }

            if (nHits > 0)
                fh_mult_Fib[ifibcount]->Fill(nHits);

            if (TofdQ > 2. && TofdQ < 3.)
            {
                fh_Fibs_vs_Tofd[ifibcount]->Fill(TofdX, iFibMax);
                // cout<<"test "<<TofdX<<"  "<<iFibMax<<endl;
            }
            if (TofdX > 1.4 && TofdX < 2.8)
            {
                fh_Fibs_vs_Tofd_y[ifibcount]->Fill(TofdY, iFibMax);
            }
            fh_Cave_position->Fill(z[ifibcount], x[ifibcount] + xposfib);

            if (test)
                cout << "fiber max: " << iFibMax << endl;

        } // end if(aHit[ifibcount])

        if (ifibcount == 13 && test)
        {
            cout << "   " << endl;
            cout << "new event! ******************************************************" << endl;
        }

    } // end for(ifibcount)

    for (Int_t i = 0; i < NOF_FIB_DET; i++)
    {
        for (Int_t j = i + 1; j < NOF_FIB_DET; j++)
        {
            if (fHitItems.at(i + DET_FI_FIRST) && fHitItems.at(j + DET_FI_FIRST))
            {
                if (TofdQ > 0. && TofdQ < 100.)
                {
                    fh_Fib_vs_Fib[i][j]->Fill(FibMax[i], FibMax[j]);
                    if (FibMax[i] > 0 && FibMax[j] > 0)
                        fh_Fib_dx[i][j]->Fill(FibMax[i], FibMax[j] - FibMax[i]);
                }
            }
        }
    }

    fNEvents += 1;
}

void R3BGlobalAnalysis::FinishEvent()
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

void R3BGlobalAnalysis::FinishTask()
{
    fh_Cave_position->Write();
    if (fMappedItems.at(DET_LOS))
    {
        fhChargeLosTofD->Write();
        fh_los_pos->Write();
    }
    if (fHitItems.at(DET_TOFD))
    {
        fh_tofd_pos->Write();
        fh_tofd_charge->Write();
        fh_TimePreviousEvent->Write();
        fh_tofd_mult->Write();
    }

    for (Int_t ifibcount = 0; ifibcount < NOF_FIB_DET; ifibcount++)
    {
        if (fMappedItems.at(ifibcount + DET_FI_FIRST))
        {
            fh_channels_Fib[ifibcount]->Write();
            fh_fibers_Fib[ifibcount]->Write();
            fh_fiber_Fib[ifibcount]->Write();
            fh_mult_Fib[ifibcount]->Write();
            fh_time_Fib[ifibcount]->Write();
            fh_xpos_Fib[ifibcount]->Write();
            fh_ypos_Fib[ifibcount]->Write();
            fh_multihit_m_Fib[ifibcount]->Write();
            fh_multihit_s_Fib[ifibcount]->Write();
            fh_ToT_m_Fib[ifibcount]->Write();
            fh_ToT_s_Fib[ifibcount]->Write();
            fh_Fib_vs_Events[ifibcount]->Write();
            fh_Fibs_vs_Events[ifibcount]->Write();
            fh_Fibs_vs_Tofd[ifibcount]->Write();
            fh_Fibs_vs_Tofd_y[ifibcount]->Write();
            fh_Fib_ToF[ifibcount]->Write();
            fh_ToF_vs_events[ifibcount]->Write();
        }
    }

    for (Int_t i = 0; i < NOF_FIB_DET; i++)
    {
        for (Int_t j = i + 1; j < NOF_FIB_DET; j++)
        {
            if (fHitItems.at(i + DET_FI_FIRST) && fHitItems.at(j + DET_FI_FIRST))
            {
                fh_Fib_vs_Fib[i][j]->Write();
                fh_Fib_dx[i][j]->Write();
            }
        }
    }
}

ClassImp(R3BGlobalAnalysis)
