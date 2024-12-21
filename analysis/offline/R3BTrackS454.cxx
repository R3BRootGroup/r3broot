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
// -----                  R3BTrackS454                -----
// -----          Created April 13th 2016 by M.Heil       -----
// ------------------------------------------------------------

/*
 * This task should fill histograms with detector variables which allow
 * to test the detectors online
 *
 */

#include "R3BCalifaMappedData.h"
#include "R3BLosCalData.h"
#include "R3BLosHitData.h"
#include "R3BLosMappedData.h"

#include "R3BBeamMonitorMappedData.h"

#include "R3BTrackS454.h"

#include "R3BSci8CalData.h"
#include "R3BSci8MappedData.h"

#include "R3BTofdCalData.h"
#include "R3BTofdHitData.h"
#include "R3BTofdMappedData.h"

#include "R3BRoluCalData.h"
#include "R3BRoluMappedData.h"

#include "R3BPaddleCalData.h"

#include "R3BPspxCalData.h"
#include "R3BPspxMappedData.h"

#include "R3BEventHeader.h"
#include "R3BTCalEngine.h"

#include "R3BBunchedFiberCalData.h"
#include "R3BBunchedFiberHitData.h"
#include "R3BBunchedFiberMappedData.h"

#include "R3BMCTrack.h"
#include "R3BTrack.h"

#include "FairLogger.h"
#include "FairRootManager.h"
#include "FairRunAna.h"
#include "FairRunOnline.h"
#include "FairRuntimeDb.h"

#include "TCanvas.h"
#include "TH1F.h"
#include "TH2F.h"

#include "TCutG.h"
#include "tracker_routines.h"

#include "TClonesArray.h"
#include "TFile.h"
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

R3BTrackS454::R3BTrackS454()
    : R3BTrackS454("Track", 1)
{
}

R3BTrackS454::R3BTrackS454(const char* name, Int_t iVerbose)
    : FairTask(name, iVerbose)
    , fTrigger(-1)
    , fTpat(-1)
    , fCuts(0)
    , fGraphCuts(0)
    , fGhost(0)
    , fPairs(0)
    , fB(-1672)
    , fSimu(0)
    , fNEvents(0)
    , fTrackItems(new TClonesArray("R3BTrack"))
    , fNofTrackItems()
{
}

R3BTrackS454::~R3BTrackS454()
{
    for (int i = 0; i < NOF_FIB_DET; i++)
    {
        delete fh_xy_Fib[i];
        delete fh_mult_Fib[i];
        delete fh_ToT_Fib[i];
    }
    delete fTrackItems;
}

InitStatus R3BTrackS454::Init()
{

    // Initialize random number:
    std::srand(std::time(0)); // use current time as seed for random generator

    LOG(info) << "R3BTrackS454::Init ";

    // try to get a handle on the EventHeader. EventHeader may not be
    // present though and hence may be null. Take care when using.

    FairRootManager* mgr = FairRootManager::Instance();
    if (NULL == mgr)
        LOG(fatal) << "FairRootManager not found";

    header = dynamic_cast<R3BEventHeader*>(mgr->GetObject("R3BEventHeader"));
    FairRunOnline* run = FairRunOnline::Instance();

    // Get objects for detectors on all levels
    fMCTrack = dynamic_cast<TClonesArray*>(mgr->GetObject("MCTrack"));
    if (fMCTrack)
        mgr->Register("MCTrack", "Monte Carlo Tracks", fMCTrack, kTRUE);

    assert(DET_MAX + 1 == sizeof(fDetectorNames) / sizeof(fDetectorNames[0]));
    printf("Have %d fiber detectors.\n", NOF_FIB_DET);
    fMappedItems.push_back(dynamic_cast<TClonesArray*>(mgr->GetObject(Form("%sMappedData", fDetectorNames[0]))));
    if (NULL == fMappedItems.at(0))
    {
        printf("Could not find mapped data for '%s'.\n", fDetectorNames[0]);
    }
    fCalItems.push_back(dynamic_cast<TClonesArray*>(mgr->GetObject(Form("%sCrystalCalData", fDetectorNames[0]))));
    if (NULL == fCalItems.at(0))
    {
        printf("Could not find Cal data for '%s'.\n", fDetectorNames[0]);
    }
    fHitItems.push_back(dynamic_cast<TClonesArray*>(mgr->GetObject(Form("%sHit", fDetectorNames[0]))));
    if (NULL == fHitItems.at(0))
    {
        printf("Could not find hit data for '%s'.\n", fDetectorNames[0]);
    }
    for (int det = 1; det < DET_MAX; det++)
    {
        fMappedItems.push_back(dynamic_cast<TClonesArray*>(mgr->GetObject(Form("%sMapped", fDetectorNames[det]))));
        if (NULL == fMappedItems.at(det))
        {
            printf("Could not find mapped data for '%s'.\n", fDetectorNames[det]);
        }
        if (det == 9)
            maxevent = mgr->CheckMaxEventNo();
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

    mgr->Register("Track", "Land", fTrackItems, kTRUE);

    //------------------------------------------------------------------------
    // graphical cuts
    //------------------------------------------------------------------------
    if (fGraphCuts)
    {
        cut_fi11_fi3a = NULL;
        cut_fi12_fi3b = NULL;
        cut_Fi13vsFi11 = NULL;
        cut_Fi10vsFi12 = NULL;

        TFile* f = TFile::Open("fiberCuts.root", "read");
        cut_Fi13vsFi11 = dynamic_cast<TCutG*>(f->Get("cut_dx_Fi13_Fi11"));
        cut_Fi10vsFi12 = dynamic_cast<TCutG*>(f->Get("cut_dx_Fi10_Fi12"));

        TFile* f3a = TFile::Open("myCutsFi3a.root", "read");
        cut_fi11_fi3a = dynamic_cast<TCutG*>(f3a->Get("cut_fi11_fi3a"));

        TFile* f3b = TFile::Open("myCutsFi3b.root", "read");
        cut_fi12_fi3b = dynamic_cast<TCutG*>(f3b->Get("cut_fi12_fi3b"));
    }
    //------------------------------------------------------------------------
    // create histograms of all detectors
    //------------------------------------------------------------------------

    char strNameC[255];
    sprintf(strNameC, "Cave_C_position");
    fh_Cave_position = new TH2F(strNameC, "", 210, -100., 200., 1000, -500., 500.);

    //-----------------------------------------------------------------------
    // BeamMonitor

    // get the theoretical calib factors for SEETRAM
    Double_t fexp = float(fsens_SEE + 9);
    Double_t fpow = float(pow(10., fexp));
    calib_SEE = 135641.7786 * fpow;
    LOG(debug) << fsens_SEE << ", " << fexp << ", " << fpow << ", " << calib_SEE << endl;

    fh_Tpat = new TH1F("Tpat", "Tpat", 20, 0, 20);
    fh_Tpat->GetXaxis()->SetTitle("Tpat value");

    fh_Trigger = new TH1F("Trigger", "Trigger all", 20, 0, 20);
    fh_Trigger->GetXaxis()->SetTitle("Trigger value");

    fh_IC = new TH1F("IC", "IC ", 1000, 0, 1000);
    fh_IC->GetXaxis()->SetTitle("spill number");
    fh_IC->GetYaxis()->SetTitle("IC counts");

    fh_SEE = new TH1F("SEETRAM", "SEETRAM ", 1000, 0, 1000);
    fh_SEE->GetXaxis()->SetTitle("spill number");
    fh_SEE->GetYaxis()->SetTitle("SEETRAM counts");

    fh_TOFDOR = new TH1F("TOFDOR", "TOFDOR ", 1000, 0, 1000);
    fh_TOFDOR->GetXaxis()->SetTitle("spill number");
    fh_TOFDOR->GetYaxis()->SetTitle("TOFDOR counts");

    //-----------------------------------------------------------------------
    // compare against MC Simulations

    fh_target_px = new TH1F("target_px", "target px ", 2000, -500., 500);
    fh_target_px->GetXaxis()->SetTitle("px / MeV/c");
    fh_target_px->GetYaxis()->SetTitle("counts");

    fh_target_py = new TH1F("target_py", "target py ", 2000, -500., 500);
    fh_target_py->GetXaxis()->SetTitle("py / MeV/c");
    fh_target_py->GetYaxis()->SetTitle("counts");

    fh_target_pz = new TH1F("target_pz", "target pz ", 30000, 0., 30000.);
    fh_target_pz->GetXaxis()->SetTitle("pz / MeV/c");
    fh_target_pz->GetYaxis()->SetTitle("counts");

    fh_target_p = new TH1F("target_p", "target p ", 30000, 0., 30000.);
    fh_target_p->GetXaxis()->SetTitle("pz / MeV/c");
    fh_target_p->GetYaxis()->SetTitle("counts");

    fh_px_He = new TH1F("px_He", " px He", 2000, -500., 500);
    fh_px_He->GetXaxis()->SetTitle("px / MeV/c");
    fh_px_He->GetYaxis()->SetTitle("counts");

    fh_py_He = new TH1F("py_He", " py He", 2000, -500., 500);
    fh_py_He->GetXaxis()->SetTitle("py / MeV/c");
    fh_py_He->GetYaxis()->SetTitle("counts");

    fh_pz_He = new TH1F("pz_He", " pz He", 30000, 0., 30000.);
    fh_pz_He->GetXaxis()->SetTitle("pz / MeV/c");
    fh_pz_He->GetYaxis()->SetTitle("counts");

    fh_p_He = new TH1F("p_He", " p He", 30000, 0., 30000.);
    fh_p_He->GetXaxis()->SetTitle("pz / MeV/c");
    fh_p_He->GetYaxis()->SetTitle("counts");

    fh_px_C = new TH1F("px_C", " px C", 2000, -500., 500);
    fh_px_C->GetXaxis()->SetTitle("px / MeV/c");
    fh_px_C->GetYaxis()->SetTitle("counts");

    fh_py_C = new TH1F("py_C", " py C", 2000, -500., 500);
    fh_py_C->GetXaxis()->SetTitle("py / MeV/c");
    fh_py_C->GetYaxis()->SetTitle("counts");

    fh_pz_C = new TH1F("pz_C", " pz C", 30000, 0., 30000.);
    fh_pz_C->GetXaxis()->SetTitle("pz / MeV/c");
    fh_pz_C->GetYaxis()->SetTitle("counts");

    fh_p_C = new TH1F("p_C", " p C", 30000, 0., 30000.);
    fh_p_C->GetXaxis()->SetTitle("pz / MeV/c");
    fh_p_C->GetYaxis()->SetTitle("counts");

    fh_target_xy = new TH2F("target_xy", "target xy ", 100, -0.5, 0.5, 100, -0.5, 0.5);
    fh_target_xy->GetXaxis()->SetTitle("x / cm");
    fh_target_xy->GetYaxis()->SetTitle("y / cm");

    fh_chi2 = new TH1F("chi2", "chi2 ", 1000, 0., 100);
    fh_chi2->GetXaxis()->SetTitle("Chi2");
    fh_chi2->GetYaxis()->SetTitle("counts");

    fh_dx = new TH1F("tracker_dx", "tracker dx ", 200, -5., 5);
    fh_dx->GetXaxis()->SetTitle("dx / cm");
    fh_dx->GetYaxis()->SetTitle("counts");

    fh_dy = new TH1F("tracker_dy", "tracker dy ", 200, -5., 5.);
    fh_dy->GetXaxis()->SetTitle("dy / cm");
    fh_dy->GetYaxis()->SetTitle("counts");

    fh_dz = new TH1F("tracker_dz", "tracker dz ", 200, -5., 5.);
    fh_dz->GetXaxis()->SetTitle("dz / cm");
    fh_dz->GetYaxis()->SetTitle("counts");

    fh_dpx = new TH1F("tracker_dpx", "tracker dpx ", 200, -10, 10);
    fh_dpx->GetXaxis()->SetTitle("dpx / percent");
    fh_dpx->GetYaxis()->SetTitle("counts");

    fh_dpy = new TH1F("tracker_dpy", "tracker dpy ", 200, -10, 10);
    fh_dpy->GetXaxis()->SetTitle("dpy / percent");
    fh_dpy->GetYaxis()->SetTitle("counts");

    fh_dpz = new TH1F("tracker_dpz", "tracker dpz ", 200, -10, 10);
    fh_dpz->GetXaxis()->SetTitle("dpz / percent");
    fh_dpz->GetYaxis()->SetTitle("counts");

    fh_dp = new TH1F("tracker_dp", "tracker dp ", 200, -10, 10);
    fh_dp->GetXaxis()->SetTitle("dp / percent");
    fh_dp->GetYaxis()->SetTitle("counts");

    fh_thetax_dpx = new TH2F("tracker_thethax_dpx", "tracker theta_x vs dpx ", 100, -50, 50, 200, -100., 100.);
    fh_thetax_dpx->GetXaxis()->SetTitle("dpx / percent");
    fh_thetax_dpx->GetYaxis()->SetTitle("theta_x / mrad");

    fh_dpy_dpx = new TH2F("tracker_dpy_dpx", "tracker dpy vs dpx ", 200, -100, 100, 200, -100., 100.);
    fh_dpy_dpx->GetXaxis()->SetTitle("dpx / percent");
    fh_dpy_dpx->GetYaxis()->SetTitle("dpy / percent");

    fh_thetay_dpy = new TH2F("tracker_thethay_dpy", "tracker theta_y vs dpy ", 100, -50, 50, 200, -100., 100.);
    fh_thetay_dpy->GetXaxis()->SetTitle("dpy / percent");
    fh_thetay_dpy->GetYaxis()->SetTitle("theta_y / mrad");

    fh_thetax_dpy = new TH2F("tracker_thethax_dpy", "tracker theta_x vs dpy ", 100, -50, 50, 200, -100., 100.);
    fh_thetax_dpy->GetXaxis()->SetTitle("dpy / percent");
    fh_thetax_dpy->GetYaxis()->SetTitle("theta_x / mrad");

    fh_thetay_dpx = new TH2F("tracker_thethay_dpx", "tracker theta_y vs dpx ", 100, -50, 50, 200, -100., 100.);
    fh_thetay_dpx->GetXaxis()->SetTitle("dpx / percent");
    fh_thetay_dpx->GetYaxis()->SetTitle("theta_y / mrad");

    fh_thetax_dpx_abs =
        new TH2F("tracker_thethax_dpx_abs", "tracker theta_x vs dpx abs ", 100, -50, 50, 200, -100., 100.);
    fh_thetax_dpx_abs->GetXaxis()->SetTitle("dpx");
    fh_thetax_dpx_abs->GetYaxis()->SetTitle("theta_x / mrad");

    fh_x_dpx = new TH2F("tracker_x_dpx", "tracker x vs dpx ", 100, -50, 50, 200, -10., 10.);
    fh_x_dpx->GetXaxis()->SetTitle("dpx / percent");
    fh_x_dpx->GetYaxis()->SetTitle("x / cm");

    fh_y_dpy = new TH2F("tracker_y_dpy", "tracker y vs dpy ", 100, -50, 50, 200, -10., 10.);
    fh_y_dpy->GetXaxis()->SetTitle("dpy / percent");
    fh_y_dpy->GetYaxis()->SetTitle("y / cm");

    // histograms for track hits
    for (Int_t i = 0; i < ndet; i++)
    {
        fh_xy[i] = new TH2F(Form("xy_%i", i), Form("xy of Det %i", i), 600, -30, 30, 1200, -60., 60.);
        fh_xy[i]->GetXaxis()->SetTitle("x / cm");
        fh_xy[i]->GetYaxis()->SetTitle("y / cm");

        fh_p_vs_x[i] = new TH2F(Form("PvsX%i", i), Form("pz vs. x of Det %i", i), 1200, -60, 60, 4000, 0., 40000.);
        fh_p_vs_x[i]->GetXaxis()->SetTitle("x position / cm");
        fh_p_vs_x[i]->GetYaxis()->SetTitle("p / MeV/c");

        fh_p_vs_x_test[i] =
            new TH2F(Form("PvsX_test%i", i), Form("p vs. x of Det %i test", i), 1200, -60, 60, 4000, 0., 40000.);
        fh_p_vs_x_test[i]->GetXaxis()->SetTitle("x position / cm");
        fh_p_vs_x_test[i]->GetYaxis()->SetTitle("p / MeV/c");
    }
    //-----------------------------------------------------------------------
    // Fiber Detectors 1-NOF_FIB_DET

    char canvName[255];
    UInt_t Nmax = 1e7;
    for (Int_t ifibcount = 0; ifibcount < NOF_FIB_DET; ifibcount++)
    {
        if (fCalItems.at(DET_FI_FIRST + ifibcount) || fHitItems.at(DET_FI_FIRST + ifibcount))
        {

            const char* detName;
            const char* detName2;
            detName = fDetectorNames[DET_FI_FIRST + ifibcount];

            LOG(debug) << "I am creating canvas " << detName << endl;

            // xy:
            fh_xy_Fib[ifibcount] =
                new TH2F(Form("%s_xy", detName), Form("%s xy", detName), 600, -30., 30., 200, -1000., 1000.);
            fh_xy_Fib[ifibcount]->GetXaxis()->SetTitle("x / cm ");
            fh_xy_Fib[ifibcount]->GetYaxis()->SetTitle("y / cm");

            fh_xy_Fib_ac[ifibcount] = new TH2F(
                Form("%s_xy_ac", detName), Form("%s xy after cuts", detName), 600, -30., 30., 200, -100., 100.);
            fh_xy_Fib_ac[ifibcount]->GetXaxis()->SetTitle("x / cm ");
            fh_xy_Fib_ac[ifibcount]->GetYaxis()->SetTitle("y / cm");

            // Multiplicity (number of hit fibers):
            fh_mult_Fib[ifibcount] = new TH1F(Form("%s_mult", detName), Form("%s # of fibers", detName), 500, 0., 500.);
            fh_mult_Fib[ifibcount]->GetXaxis()->SetTitle("Multiplicity");
            fh_mult_Fib[ifibcount]->GetYaxis()->SetTitle("Counts");

            fh_mult_Fib_ac[ifibcount] =
                new TH1F(Form("%s_mult_ac", detName), Form("%s # of fibers after cuts", detName), 500, 0., 500.);
            fh_mult_Fib_ac[ifibcount]->GetXaxis()->SetTitle("Multiplicity");
            fh_mult_Fib_ac[ifibcount]->GetYaxis()->SetTitle("Counts");

            // ToT MAPMT:
            fh_ToT_Fib[ifibcount] =
                new TH2F(Form("%s_tot_m", detName), Form("%s ToT of MAPMT", detName), 600, -30., 30, 400, 0., 400.);
            fh_ToT_Fib[ifibcount]->GetXaxis()->SetTitle("Fiber x / cm");
            fh_ToT_Fib[ifibcount]->GetYaxis()->SetTitle("ToT / ns");

            fh_ToT_Fib_ac[ifibcount] = new TH2F(Form("%s_tot_m_ac", detName),
                                                Form("%s ToT of MAPMT after cuts", detName),
                                                600,
                                                -30.,
                                                30,
                                                400,
                                                0.,
                                                400.);
            fh_ToT_Fib_ac[ifibcount]->GetXaxis()->SetTitle("Fiber x / cm");
            fh_ToT_Fib_ac[ifibcount]->GetYaxis()->SetTitle("ToT / ns");

            // ToF Tofd -> Fiber:
            fh_Fib_ToF[ifibcount] = new TH2F(
                Form("%s_tof", detName), Form("%s ToF Tofd to Fiber", detName), 600, -30., 30, 10000, -1100., 1100.);
            fh_Fib_ToF[ifibcount]->GetYaxis()->SetTitle("ToF / ns");
            fh_Fib_ToF[ifibcount]->GetXaxis()->SetTitle("x / cm");

            fh_Fib_ToF_ac[ifibcount] = new TH2F(Form("%s_tof_ac", detName),
                                                Form("%s ToF Tofd to Fiber after cuts", detName),
                                                600,
                                                -30.,
                                                30,
                                                10000,
                                                -1000.,
                                                1000.);
            fh_Fib_ToF_ac[ifibcount]->GetYaxis()->SetTitle("ToF / ns");
            fh_Fib_ToF_ac[ifibcount]->GetXaxis()->SetTitle("x / cm");

            // Time:
            fh_Fib_Time[ifibcount] =
                new TH2F(Form("%s_time", detName), Form("%s Time", detName), 600, -30., 30, 4000, -4000., 4000.);
            fh_Fib_Time[ifibcount]->GetYaxis()->SetTitle("Time / ns");
            fh_Fib_Time[ifibcount]->GetXaxis()->SetTitle("x / cm");

            fh_Fib_Time_ac[ifibcount] = new TH2F(
                Form("%s_time_ac", detName), Form("%s Time after cuts", detName), 600, -30., 30, 4000, -4000., 4000.);
            fh_Fib_Time_ac[ifibcount]->GetYaxis()->SetTitle("Time / ns");
            fh_Fib_Time_ac[ifibcount]->GetXaxis()->SetTitle("x / cm");

            // ToF Tofd -> Fiber vs. event number:
            fh_ToF_vs_Events[ifibcount] = new TH2F(Form("%s_tof_vs_events", detName),
                                                   Form("%s ToF Tofd to Fiber vs event number", detName),
                                                   10000,
                                                   0,
                                                   Nmax,
                                                   2200,
                                                   -5100,
                                                   5100);
            fh_ToF_vs_Events[ifibcount]->GetYaxis()->SetTitle("ToF / ns");
            fh_ToF_vs_Events[ifibcount]->GetXaxis()->SetTitle("event number");

            fh_ToF_vs_Events_ac[ifibcount] = new TH2F(Form("%s_tof_vs_events_ac", detName),
                                                      Form("%s ToF Tofd to Fiber vs event number after cuts", detName),
                                                      10000,
                                                      0,
                                                      Nmax,
                                                      2200,
                                                      -5100,
                                                      5100);
            fh_ToF_vs_Events_ac[ifibcount]->GetYaxis()->SetTitle("ToF / ns");
            fh_ToF_vs_Events_ac[ifibcount]->GetXaxis()->SetTitle("event number");

            // hit fiber number vs. event number:
            fh_Fib_vs_Events[ifibcount] = new TH2F(Form("%s_fib_vs_event", detName),
                                                   Form("%s Fiber # vs. Event #", detName),
                                                   10000,
                                                   0,
                                                   Nmax,
                                                   600,
                                                   -30.,
                                                   30.);
            fh_Fib_vs_Events[ifibcount]->GetYaxis()->SetTitle("Fiber x / cm");
            fh_Fib_vs_Events[ifibcount]->GetXaxis()->SetTitle("Event number");

            fh_Fib_vs_Events_ac[ifibcount] = new TH2F(Form("%s_fib_vs_event_ac", detName),
                                                      Form("%s Fiber # vs. Event # after cuts", detName),
                                                      10000,
                                                      0,
                                                      Nmax,
                                                      600,
                                                      -30.,
                                                      30.);
            fh_Fib_vs_Events_ac[ifibcount]->GetYaxis()->SetTitle("Fiber x / cm");
            fh_Fib_vs_Events_ac[ifibcount]->GetXaxis()->SetTitle("Event number");

            // hit fiber number vs. TofD position:
            fh_Fibs_vs_Tofd[ifibcount] = new TH2F(Form("%s_fib_vs_TofdX", detName),
                                                  Form("%s Fiber # vs. Tofd x-pos", detName),
                                                  200,
                                                  -100,
                                                  100,
                                                  600,
                                                  -30.,
                                                  30.);
            fh_Fibs_vs_Tofd[ifibcount]->GetYaxis()->SetTitle("Fiber x / cm");
            fh_Fibs_vs_Tofd[ifibcount]->GetXaxis()->SetTitle("Tofd x / cm");

            fh_Fibs_vs_Tofd_ac[ifibcount] = new TH2F(Form("%s_fib_vs_TofdX_ac", detName),
                                                     Form("%s Fiber # vs. Tofd x-pos after cuts", detName),
                                                     200,
                                                     -100,
                                                     100,
                                                     600,
                                                     -30.,
                                                     30.);
            fh_Fibs_vs_Tofd_ac[ifibcount]->GetYaxis()->SetTitle("Fiber x / cm");
            fh_Fibs_vs_Tofd_ac[ifibcount]->GetXaxis()->SetTitle("Tofd x / cm");

            // hit fiber vs. fiber position:

        } // end if(Mapped)

    } // end for(ifibcount)
    fh_Fib13_vs_Fib11 = new TH2F("fib13_vs_fib11", "Fiber 13 vs. Fiber 11", 1000, -50, 50, 1000, -50., 50.);
    fh_Fib13_vs_Fib11->GetYaxis()->SetTitle("Fiber13");
    fh_Fib13_vs_Fib11->GetXaxis()->SetTitle("Fiber11");

    fh_Fib11_vs_Fib3a = new TH2F("fib11_vs_fib3a", "Fiber 11 vs. Fiber 3a", 1000, -50, 50, 1000, -50., 50.);
    fh_Fib11_vs_Fib3a->GetYaxis()->SetTitle("Fiber11");
    fh_Fib11_vs_Fib3a->GetXaxis()->SetTitle("Fiber3a");

    fh_Fib10_vs_Fib12 = new TH2F("fib10_vs_fib12", "Fiber 10 vs. Fiber 12", 1000, -50, 50, 1000, -50., 50.);
    fh_Fib10_vs_Fib12->GetYaxis()->SetTitle("Fiber10");
    fh_Fib10_vs_Fib12->GetXaxis()->SetTitle("Fiber12");

    fh_Fib12_vs_Fib3b = new TH2F("fib12_vs_fib3b", "Fiber 12 vs. Fiber 3b", 1000, -50, 50, 1000, -50., 50.);
    fh_Fib12_vs_Fib3b->GetYaxis()->SetTitle("Fiber12");
    fh_Fib12_vs_Fib3b->GetXaxis()->SetTitle("Fiber3b");

    // dx between fibers vs x
    fh_Fib13_vs_Fib11_dx = new TH2F("fib13_fib11_dx", "dx of Fiber 13 and Fiber 11", 1000, -50, 50, 1000, -50., 50.);
    fh_Fib13_vs_Fib11_dx->GetYaxis()->SetTitle("xFi13 - xFi11 / cm");
    fh_Fib13_vs_Fib11_dx->GetXaxis()->SetTitle("x Fi11 / cm");

    fh_Fib11_vs_Fib3a_dx = new TH2F("fib11_fib3a_dx", "dx of Fiber 11 and Fiber 3a", 1000, -50, 50, 1000, -50., 50.);
    fh_Fib11_vs_Fib3a_dx->GetYaxis()->SetTitle("xFi11 - xFi3a / cm");
    fh_Fib11_vs_Fib3a_dx->GetXaxis()->SetTitle("x Fi3a / cm");

    fh_Fib10_vs_Fib12_dx = new TH2F("fib10_fib12_dx", "dx of Fiber 10 and Fiber 12", 1000, -50, 50, 1000, -50., 50.);
    fh_Fib10_vs_Fib12_dx->GetYaxis()->SetTitle("xFi10 - xFi12 / cm");
    fh_Fib10_vs_Fib12_dx->GetXaxis()->SetTitle("x Fi12 / cm");

    fh_Fib12_vs_Fib3b_dx = new TH2F("fib12_fib3b_dx", "dx of Fiber 12 and Fiber 3b", 1000, -50, 50, 1000, -50., 50.);
    fh_Fib12_vs_Fib3b_dx->GetYaxis()->SetTitle("xFi12 - xFi3b / cm");
    fh_Fib12_vs_Fib3b_dx->GetXaxis()->SetTitle("x Fi3b / cm");

    fh_Fib13_vs_Fib11_back =
        new TH2F("fib13_vs_fib11_back", "Fiber 13 vs. Fiber 11 back", 1000, -50, 50, 1000, -50., 50.);
    fh_Fib13_vs_Fib11_back->GetYaxis()->SetTitle("Fiber13");
    fh_Fib13_vs_Fib11_back->GetXaxis()->SetTitle("Fiber11");

    fh_Fib11_vs_Fib3a_back =
        new TH2F("fib11_vs_fib3a_back", "Fiber 11 vs. Fiber 3a back", 1000, -50, 50, 1000, -50., 50.);
    fh_Fib11_vs_Fib3a_back->GetYaxis()->SetTitle("Fiber11");
    fh_Fib11_vs_Fib3a_back->GetXaxis()->SetTitle("Fiber3a");

    fh_Fib10_vs_Fib12_back =
        new TH2F("fib10_vs_fib12_back", "Fiber 10 vs. Fiber 12 back", 1000, -50, 50, 1000, -50., 50.);
    fh_Fib10_vs_Fib12_back->GetYaxis()->SetTitle("Fiber10");
    fh_Fib10_vs_Fib12_back->GetXaxis()->SetTitle("Fiber12");

    fh_Fib12_vs_Fib3b_back =
        new TH2F("fib12_vs_fib3b_back", "Fiber 12 vs. Fiber 3b back", 1000, -50, 50, 1000, -50., 50.);
    fh_Fib12_vs_Fib3b_back->GetYaxis()->SetTitle("Fiber12");
    fh_Fib12_vs_Fib3b_back->GetXaxis()->SetTitle("Fiber3b");

    // dx between fibers vs x
    fh_Fib13_vs_Fib11_dx_back =
        new TH2F("fib13_fib11_dx_back", "dx of Fiber 13 and Fiber 11 back", 1000, -50, 50, 1000, -50., 50.);
    fh_Fib13_vs_Fib11_dx_back->GetYaxis()->SetTitle("xFi13 - xFi11 / cm");
    fh_Fib13_vs_Fib11_dx_back->GetXaxis()->SetTitle("x Fi11 / cm");

    fh_Fib11_vs_Fib3a_dx_back =
        new TH2F("fib11_fib3a_dx_back", "dx of Fiber 11 and Fiber 3a back", 1000, -50, 50, 1000, -50., 50.);
    fh_Fib11_vs_Fib3a_dx_back->GetYaxis()->SetTitle("xFi11 - xFi3a / cm");
    fh_Fib11_vs_Fib3a_dx_back->GetXaxis()->SetTitle("x Fi3a / cm");

    fh_Fib10_vs_Fib12_dx_back =
        new TH2F("fib10_fib12_dx_back", "dx of Fiber 10 and Fiber 12 back", 1000, -50, 50, 1000, -50., 50.);
    fh_Fib10_vs_Fib12_dx_back->GetYaxis()->SetTitle("xFi10 - xFi12 / cm");
    fh_Fib10_vs_Fib12_dx_back->GetXaxis()->SetTitle("x Fi12 / cm");

    fh_Fib12_vs_Fib3b_dx_back =
        new TH2F("fib12_fib3b_dx_back", "dx of Fiber 12 and Fiber 3b back", 1000, -50, 50, 1000, -50., 50.);
    fh_Fib12_vs_Fib3b_dx_back->GetYaxis()->SetTitle("xFi12 - xFi3b / cm");
    fh_Fib12_vs_Fib3b_dx_back->GetXaxis()->SetTitle("x Fi3b / cm");

    //---------------------------------------------------------------------------------------------------
    // TofD detector

    if (fHitItems.at(DET_TOFD) || fCalItems.at(DET_TOFD))
    {

        // xy:
        fh_xy_tofd = new TH2F("tofd_xy", "tofd xy", 200, -100., 100., 200, -100., 100.);
        fh_xy_tofd->GetXaxis()->SetTitle("x / cm ");
        fh_xy_tofd->GetYaxis()->SetTitle("y / cm");

        fh_xy_tofd_ac = new TH2F("tofd_xy_ac", "tofd xy after cuts", 200, -100., 100., 200, -100., 100.);
        fh_xy_tofd_ac->GetXaxis()->SetTitle("x / cm ");
        fh_xy_tofd_ac->GetYaxis()->SetTitle("y / cm");

        fh_tofd_charge = new TH1F("tofd_Q", "Charge of Tofd", 200, 0., 20.);
        fh_tofd_charge->GetXaxis()->SetTitle("x / cm ");
        fh_tofd_charge->GetYaxis()->SetTitle("y / cm");

        fh_tofd_charge_ac = new TH1F("tofd_Q_ac", "Charge of Tofd after cuts", 200, 0., 20.);
        fh_tofd_charge_ac->GetXaxis()->SetTitle("x / cm ");
        fh_tofd_charge_ac->GetYaxis()->SetTitle("y / cm");

        fh_tofd_mult = new TH1F("tofd_mult", "ToFD multiplicits ", 100, 0, 100);
        fh_tofd_mult->GetXaxis()->SetTitle("multiplicity");
        fh_tofd_mult->GetYaxis()->SetTitle("counts");

        fh_tofd_mult_ac = new TH1F("tofd_mult_ac", "ToFD multiplicits after cuts", 100, 0, 100);
        fh_tofd_mult_ac->GetXaxis()->SetTitle("multiplicity");
        fh_tofd_mult_ac->GetYaxis()->SetTitle("counts");

        fh_TimePreviousEvent = new TH1F("TimePreviousEvent", "Time between 2 particles ", 3000, 0, 3000);
        fh_TimePreviousEvent->GetXaxis()->SetTitle("time / ns");
        fh_TimePreviousEvent->GetYaxis()->SetTitle("counts");

        fh_tofd_time = new TH1F("tofd_time", "Tofd times ", 40000, -2000, 2000);
        fh_tofd_time->GetXaxis()->SetTitle("time / ns");
        fh_tofd_time->GetYaxis()->SetTitle("counts");

        fh_tofd_time_ac = new TH1F("tofd_time_ac", "Tofd times after cut", 40000, -2000, 2000);
        fh_tofd_time_ac->GetXaxis()->SetTitle("time / ns");
        fh_tofd_time_ac->GetYaxis()->SetTitle("counts");

        fh_tofd_q2_vs_q1 = new TH2F("tofd_q2_vs_q1", "tofd q2 vs. q1", 500, 0., 50., 500, 0., 50.);
        fh_tofd_q2_vs_q1->GetXaxis()->SetTitle("q1");
        fh_tofd_q2_vs_q1->GetYaxis()->SetTitle("q2");

        fh_tofd_q2_vs_q1_ac = new TH2F("tofd_q2_vs_q1_ac", "tofd q2 vs. q1 after cut", 500, 0., 50., 500, 0., 50.);
        fh_tofd_q2_vs_q1_ac->GetXaxis()->SetTitle("q1");
        fh_tofd_q2_vs_q1_ac->GetYaxis()->SetTitle("q2");
    }

    if (fMappedItems.at(DET_CALIFA))
    {
        fh_califa_energy = new TH2F("fh_califa_energy", "Califa E vs crystal id", 2000, 0, 2000, 1000, 0., 1000.);
        fh_califa_energy->GetYaxis()->SetTitle("Energy / MeV");
        fh_califa_energy->GetXaxis()->SetTitle("Crystal #");
    }

    // Analysis
    fh_chiy_vs_chix = new TH2F("chiy_vs_chix", "chi y vs. q1", 500, 0., 500., 500, 0., 500.);
    fh_chiy_vs_chix->GetXaxis()->SetTitle("chi x");
    fh_chiy_vs_chix->GetYaxis()->SetTitle("chi y");

    fh_theta26 = new TH1F("theta26", "theta 26 ", 500, 0., 5);
    fh_theta26->GetXaxis()->SetTitle("angle / degree");
    fh_theta26->GetYaxis()->SetTitle("counts");

    fh_theta26_simu = new TH1F("theta26_simu", "theta 26 simulation", 500, 0., 5);
    fh_theta26_simu->GetXaxis()->SetTitle("angle / degree");
    fh_theta26_simu->GetYaxis()->SetTitle("counts");

    fh_Erel = new TH1F("Erel", "Erel ", 5000, 0., 100);
    fh_Erel->GetXaxis()->SetTitle("Erel / MeV");
    fh_Erel->GetYaxis()->SetTitle("counts");

    fh_Erel_simu = new TH1F("Erel_simu", "Erel simulation", 6000, -10., 50);
    fh_Erel_simu->GetXaxis()->SetTitle("Erel / MeV");
    fh_Erel_simu->GetYaxis()->SetTitle("counts");

    fh_dErel_vs_x = new TH2F("dErel_vs_x", "delta Erel vs. x", 200, 0., 200., 100, 0., 10.);
    fh_dErel_vs_x->GetXaxis()->SetTitle("TofD x / cm");
    fh_dErel_vs_x->GetYaxis()->SetTitle("Delta Erel / MeV");

    fh_dErel_vs_y = new TH2F("dErel_vs_y", "delta Erel vs. y", 200, -100., 100., 100, -5., 5.);
    fh_dErel_vs_y->GetXaxis()->SetTitle("TofD y / cm");
    fh_dErel_vs_y->GetYaxis()->SetTitle("Delta Erel / MeV");

    // -------------------------------------------------------------------------
    // Rene's tracker.

    if (tracker)
        init_from_cpp_();

    return kSUCCESS;
}

void R3BTrackS454::Exec(Option_t* option)
{
    if (fNEvents / 10000. == (int)fNEvents / 10000)
        std::cout << "\rEvents: " << fNEvents << " / " << maxevent << " (" << (int)(fNEvents * 100. / maxevent)
                  << " %) "
                  << " Tofd: " << counterTofd << " tracked: " << counter2 << " chix: " << counter3
                  << " chiy: " << counter4 << std::flush;
    // cout << "Event: " << fNEvents << endl;
    fNEvents += 1;

    Bool_t debug2 = false;

    FairRootManager* mgr = FairRootManager::Instance();
    if (NULL == mgr)
        LOG(fatal) << "FairRootManager not found";

    if (header)
    {
        time = header->GetTimeStamp();
        //		if (time > 0) cout << "header time: " << time << endl;
        if (time_start == 0 && time > 0)
        {
            time_start = time;
            fNEvents_start = fNEvents;
            // cout << "Start event number " << fNEvents_start << endl;
        }

        if (header->GetTrigger() == 12)
        {
            // spill start in nsec
            // cout << "spill start" << endl;
            num_spills++;
            // cout << "Spill: " << num_spills << endl;
        }
        if (header->GetTrigger() == 13)
        {
            // spill end  in nsec
            // cout << "spill stop" << endl;
        }

        fh_Trigger->Fill(header->GetTrigger());
        //   check for requested trigger (Todo: should be done globablly / somewhere else)
        if ((fTrigger >= 0) && (header) && (header->GetTrigger() != fTrigger))
        {
            counterWrongTrigger++;
            return;
        }

        Int_t tpatbin;
        for (int i = 0; i < 16; i++)
        {
            tpatbin = (header->GetTpat() & (1 << i));
            if (tpatbin != 0)
                fh_Tpat->Fill(i + 1);
        }

        // fTpat = 1-16; fTpat_bit = 0-15
        Int_t fTpat_bit = fTpat - 1;
        Int_t itpat;
        Int_t tpatvalue;
        if (fTpat_bit >= 0)
        {
            itpat = header->GetTpat();
            tpatvalue = (itpat && (1 < fTpat_bit)) >> fTpat_bit;
            if (tpatvalue == 0)
            {
                counterWrongTpat++;
                return;
            }
        }
    }

    if (fMappedItems.at(DET_BMON))
    {
        unsigned long IC;
        unsigned long SEETRAM_raw;
        Double_t SEETRAM = 0.;
        unsigned long TOFDOR = 0;

        auto detBmon = fMappedItems.at(DET_BMON);
        Int_t nHitsbm = detBmon->GetEntriesFast();
        // cout<<"Bmon hits: "<<nHitsbm<<endl;

        for (Int_t ihit = 0; ihit < nHitsbm; ihit++)
        {
            R3BBeamMonitorMappedData* hit = dynamic_cast<R3BBeamMonitorMappedData*>(detBmon->At(ihit));
            if (!hit)
                continue;

            IC = hit->GetIC(); // negative values if offset not high enough
            counts_IC += (double)IC;

            SEETRAM_raw = hit->GetSEETRAM(); // raw counts
            if (SEETRAM_raw > 0)
            {
                SEETRAM = (double)SEETRAM_raw * calib_SEE; // calibrated SEETRAM counts
                counts_SEE += SEETRAM;

                TOFDOR = hit->GetTOFDOR(); // only positive values possible
                counts_TofD += TOFDOR;

                if (see_start == 0)
                {
                    see_start = SEETRAM;
                    ic_start = IC;
                    tofdor_start = TOFDOR;
                }

                // cout << "time " << time << endl;
                // cout << "IC   " << IC << "  " << counts_IC << "  " << endl;
                // cout << "SEE  " << SEETRAM_raw << "  " << counts_SEE << "  " << SEETRAM << endl;
                // cout << "number of 16O: " << SEETRAM - see_start << "  " << see_start << endl;
                // cout << "TofD " << TOFDOR << "  " << counts_TofD << "  " << endl;
            }

            // IC:
            Int_t yIC = IC - ic_start;
            fh_IC->Fill(num_spills, yIC);

            // SEETRAM:
            Int_t ySEE = SEETRAM - see_start;
            fh_SEE->Fill(num_spills, ySEE);
            // Double_t ySEE_part = (SEETRAM-see_mem)*fNorm*1.e+3-see_offset*calib_SEE;

            // TOFDOR:
            Int_t yTOFDOR = TOFDOR - tofdor_start;
            fh_TOFDOR->Fill(num_spills, yTOFDOR);
        }
    }

    Bool_t RoluCut = false;
    if (fMappedItems.at(DET_ROLU))
    {
        // rolu
        auto detRolu = fMappedItems.at(DET_ROLU);
        Int_t nHitsRolu = detRolu->GetEntriesFast();
        // cout<<"ROLU hits: "<<nHitsRolu<<endl;

        for (Int_t ihit = 0; ihit < nHitsRolu; ihit++)
        {
            R3BRoluMappedData* hitRolu = dynamic_cast<R3BRoluMappedData*>(detRolu->At(ihit));
            if (!hitRolu)
                continue;

            // channel numbers are stored 1-based (1..n)
            Int_t iDet = hitRolu->GetDetector(); // 1..
            Int_t iCha = hitRolu->GetChannel();  // 1..
            RoluCut = true;
        }
    }
    if (RoluCut)
    {
        // cout << "ROLU cut applied !!!" << endl;
        counterRolu++;
        return;
    }

    Bool_t CalifaCluster = false;
    if (fMappedItems.at(DET_CALIFA))
    {
        // CALIFA
        auto detCalifa = fMappedItems.at(DET_CALIFA);
        Int_t nHitsCalifa = detCalifa->GetEntriesFast();
        // cout<<"Califa hits: "<<nHitsCalifa<<endl;

        for (Int_t ihit = 0; ihit < nHitsCalifa; ihit++)
        {
            R3BCalifaMappedData* hitCalifa = dynamic_cast<R3BCalifaMappedData*>(detCalifa->At(ihit));
            if (!hitCalifa)
                continue;

            Int_t Crystal = hitCalifa->GetCrystalId();
            Int_t Energy = hitCalifa->GetEnergy();
            // cout << "Califa: " << Crystal << " Energy: " << Energy << endl;
            if (Energy > 0)
            {
                fh_califa_energy->Fill(Crystal, Energy);
                CalifaCluster = true;
            }
        }
    }
    if (CalifaCluster)
    {
        counterCalifa++;
        //		return;
    }
    if (fMCTrack)
    {
        // read in Monte Carlo Track parameter

        Int_t nHitsMCTrack = fMCTrack->GetEntriesFast();
        // cout << "MCTrack hits: " << nHitsMCTrack << endl;

        for (Int_t l = 0; l < nHitsMCTrack; l++)
        {
            R3BMCTrack* aTrack = dynamic_cast<R3BMCTrack*>(fMCTrack->At(l));

            Int_t PID = aTrack->GetPdgCode();
            Int_t mother = aTrack->GetMotherId();
            LOG(debug) << "PID " << PID << endl;
            if (mother < 0)
            {
                if (PID == 1000020040)
                {
                    // 4He
                    XHes = aTrack->GetStartX();
                    YHes = aTrack->GetStartY();
                    ZHes = aTrack->GetStartZ();
                    THes = aTrack->GetStartT();

                    pHexs = aTrack->GetPx() * 1000.;
                    pHeys = aTrack->GetPy() * 1000.;
                    pHezs = aTrack->GetPz() * 1000.;
                    pHes = sqrt((pHexs * pHexs) + (pHeys * pHeys) + (pHezs * pHezs));

                    LOG(debug) << "******************************************" << endl;
                    LOG(debug) << "Track In 4He"
                               << "x " << XHes << " y " << YHes << " z " << ZHes << endl;
                    LOG(debug) << "px " << pHexs << " py " << pHeys << " z " << pHezs << endl;
                }
                if (PID == 1000060120)
                {
                    // 12C
                    XCs = aTrack->GetStartX();
                    YCs = aTrack->GetStartY();
                    ZCs = aTrack->GetStartZ();
                    TCs = aTrack->GetStartT();

                    pCxs = aTrack->GetPx() * 1000.;
                    pCys = aTrack->GetPy() * 1000.;
                    pCzs = aTrack->GetPz() * 1000.;
                    pCs = sqrt((pCxs * pCxs) + (pCys * pCys) + (pCzs * pCzs));

                    LOG(debug) << "******************************************" << endl;
                    LOG(debug) << "Track In 12C"
                               << "x " << XCs << " y " << YCs << " z " << ZCs << endl;
                    LOG(debug) << "px " << pCxs << " py " << pCys << " z " << pCzs << endl;
                }
                if (PID == 1000080160)
                {
                    // 12C
                    Xf = aTrack->GetStartX();
                    Yf = aTrack->GetStartY();
                    Zf = aTrack->GetStartZ();

                    Pxf = aTrack->GetPx() * 1000.;
                    Pyf = aTrack->GetPy() * 1000.;
                    Pzf = aTrack->GetPz() * 1000.;
                    Pf_tot = sqrt((Pxf * Pxf) + (Pyf * Pyf) + (Pzf * Pzf));

                    LOG(debug) << "******************************************" << endl;
                    LOG(debug) << "Track In 16O"
                               << "x " << Xf << " y " << Yf << " z " << Zf << endl;
                    LOG(debug) << "px " << Pxf << " py " << Pyf << " z " << Pzf << endl;
                }
            }
        }
    }

    Double_t xTest = 0.;
    Double_t yTest = 0.;

    Int_t max = 10000;
    Int_t detector[max];
    Double_t xdet[max];
    Double_t ydet[max];
    Double_t zdet[max];
    Int_t qdet[max + 2];
    qdet[max] = 2;
    qdet[max + 1] = 6;
    Double_t xFi13[max];
    Double_t yFi13[max];
    Double_t qFi13[max];
    Double_t tFi13[max];
    Bool_t fFi13[max];
    Double_t xFi12[max];
    Double_t yFi12[max];
    Double_t qFi12[max];
    Double_t tFi12[max];
    Bool_t fFi12[max];
    Double_t xFi11[max];
    Double_t yFi11[max];
    Double_t qFi11[max];
    Double_t tFi11[max];
    Bool_t fFi11[max];
    Double_t xFi10[max];
    Double_t yFi10[max];
    Double_t qFi10[max];
    Double_t tFi10[max];
    Bool_t fFi10[max];
    Double_t xFi3a[max];
    Double_t yFi3a[max];
    Double_t qFi3a[max];
    Double_t tFi3a[max];
    Bool_t fFi3a[max];
    Double_t xFi3b[max];
    Double_t yFi3b[max];
    Double_t qFi3b[max];
    Double_t tFi3b[max];
    Bool_t fFi3b[max];
    for (int i = 0; i < max; i++)
    {
        fFi13[i] = false;
        fFi12[i] = false;
        fFi11[i] = false;
        fFi10[i] = false;
        fFi3a[i] = false;
        fFi3b[i] = false;
    }

    countdet = 0;

    Double_t track[12]; // x,y,z, px, py, pz
    Double_t chi[6];    // x,y,z, px, py, pz

    Int_t n_det = 10;
    if (fGhost)
        n_det = 11;

    Double_t x[n_det];
    Double_t y[n_det];
    Double_t z[n_det];
    Double_t q[n_det];
    Double_t t[n_det];
    Double_t x1[n_det];
    Double_t y1[n_det];
    Double_t z1[n_det];
    Double_t q1[n_det];
    Double_t t1[n_det];
    Double_t x2[n_det];
    Double_t y2[n_det];
    Double_t z2[n_det];
    Double_t q2[n_det];
    Double_t t2[n_det];
    Double_t xMax[n_det];
    Double_t yMax[n_det];
    Double_t zMax[n_det];
    Double_t qMax[n_det];
    Double_t tMax[n_det];

    Double_t xTrack[n_det];
    Double_t yTrack[n_det];
    Double_t zTrack[n_det];
    Double_t qTrack[n_det];

    Bool_t pat1[2 * n_det];
    Bool_t pat2[2 * n_det];

    Int_t id, id1, id2;

    Int_t det = 0;
    Int_t det1 = 0;
    Int_t det2 = 0;

    Double_t target[3];
    Double_t track1_best[6]; // x,y,z, px, py, pz
    Double_t track2_best[6]; // x,y,z, px, py, pz
    Double_t chi_best[6];    // chi2, chi2_red

    Double_t chi2;
    Double_t chi2_best;

    target[0] = 0.;
    target[1] = 0.;
    target[2] = 0.;

    for (int i = 0; i < n_det; i++)
    {
        x[i] = 0.;
        y[i] = 0.;
        z[i] = 0.;
        q[i] = 0.;
        t[i] = -1000.;

        x1[i] = 0.;
        y1[i] = 0.;
        z1[i] = 0.;
        q1[i] = 0.;
        t1[i] = -1000.;

        x2[i] = 0.;
        y2[i] = 0.;
        z2[i] = 0.;
        q2[i] = 0.;
        t2[i] = -1000.;

        xMax[i] = -1000.;
        yMax[i] = -1000.;
        zMax[i] = -1000.;
        qMax[i] = -1000.;
        tMax[i] = -1000.;
    }
    pHex = 0.;
    pHey = 0.;
    pHez = 0.;
    pCx = 0.;
    pCy = 0.;
    pCz = 0.;

    // is also number of ifibcount
    Int_t fi3a = 0;
    Int_t fi3b = 1;
    Int_t fi10 = 2;
    Int_t fi11 = 3;
    Int_t fi12 = 4;
    Int_t fi13 = 5;
    Int_t tofd1r = 6;
    Int_t tofd1l = 7;
    Int_t tofd2r = 8;
    Int_t tofd2l = 9;
    Int_t ghost = 10;

    Double_t tof = 0.;
    Bool_t pair = false;
    Bool_t single = false;
    Double_t tStart = 0.;
    Bool_t first = true;
    Bool_t alpha = false;
    Bool_t carbon = false;
    Double_t x_4He = 0.;
    Double_t y_4He = 0.;
    Double_t z_4He = 0.;
    Double_t x_12C = 0.;
    Double_t y_12C = 0.;
    Double_t z_12C = 0.;
    Bool_t debug_in = false;

    for (Int_t i = 0; i < 10; i++)
    {
        tPrev[i] = -1000.;
        detPrev[i] = -1;
    }

    Bool_t coinc = false;

    auto detTofd = fHitItems.at(DET_TOFD);
    Int_t nHits = detTofd->GetEntriesFast();

    if (nHits > 0)
    {
        fh_tofd_mult->Fill(nHits);
        counterTofd++;
        if (debug_in)
        {
            cout << "********************************" << endl;
            cout << "ToFD hits: " << nHits << endl;
        }
    }
    if (fB != -1672 && nHits > 100)
        return;
    if (nHits > 100)
        return;

    Int_t multTofd = 0;

    // loop over ToFD
    for (Int_t ihit = 0; ihit < nHits; ihit++)
    {

        R3BTofdHitData* hitTofd = dynamic_cast<R3BTofdHitData*>(detTofd->At(ihit));
        pair = false;

        if (IS_NAN(hitTofd->GetTime()))
            continue;
        if (debug_in)
        {
            cout << "Hit " << ihit << " of " << nHits << " charge " << hitTofd->GetEloss() << " time "
                 << hitTofd->GetTime() << endl;
        }
        Double_t ttt = hitTofd->GetTime();
        fh_tofd_time->Fill(ttt);
        if (fCuts && (ttt < -100. || ttt > 100.) && !fSimu)
        { // trigger window -1500, 1500
            if (debug_in)
                cout << "No trigger particle!" << endl;
            continue;
        }
        Double_t qqq = hitTofd->GetEloss(); // / 1.132;
        fh_tofd_charge->Fill(qqq);

        Double_t xxx = hitTofd->GetX();
        Double_t yyy = hitTofd->GetY();
        Double_t y_corr = 0.;
        Double_t randx;
        //        randx = (std::rand() / (float)RAND_MAX) * 2.8 - 1.4;
        randx = (std::rand() / (float)RAND_MAX) - 0.5;
        fh_xy_tofd->Fill(xxx + randx, yyy);
        // first looking for the right charge
        if (fB == -1102 && !fSimu)
        {
            //            if (qqq < 10. || qqq > 14.)
            //                continue;
            if (xxx >= 46. && xxx < 47.)
                y_corr = 0.08;
            if (xxx >= 47. && xxx < 48.)
                y_corr = 0.04;
            if (xxx >= 49. && xxx < 50.)
                y_corr = 0.05;
        }
        if (fB == -1292 && !fSimu)
        {
            //            if (qqq < 8.5 || qqq > 10.)
            //                continue;
            if (xxx > 30. && xxx < 31.)
                y_corr = 0.03;
            if (xxx >= 32. && xxx < 33.)
                y_corr = -0.16;
            if (xxx >= 33. && xxx < 34.)
                y_corr = 0.0;
            if (xxx >= 35. && xxx < 36.)
                y_corr = 0.12;
        }
        if (fB == -1482 && !fSimu)
        {
            if (qqq < 1. || qqq > 9.)
                continue;
            if (xxx >= 15. && xxx < 16.)
                y_corr = -0.01;
            if (xxx >= 16. && xxx < 17.)
                y_corr = 0.00;
            if (xxx >= 18. && xxx < 19.)
                y_corr = -0.04;
        }
        if (fB == -1862 && !fSimu)
        {
            //            if (qqq < 6.5 || qqq > 8.)
            //                continue;
            if (xxx >= -17. && xxx < -16.)
                y_corr = 0.12;
            if (xxx >= -16. && xxx < -15.)
                y_corr = 0.08;
            if (xxx >= -15. && xxx < -14.)
                y_corr = 0.0;
        }
        if (fB == -2052 && !fSimu)
        {
            //            if (qqq < 6. || qqq > 7.)
            //                continue;
            if (xxx >= -34. && xxx < -33.)
                y_corr = 0.07;
            if (xxx >= -33. && xxx < -32.)
                y_corr = 0.07;
            if (xxx >= -32. && xxx < -31.)
                y_corr = -0.06;
            if (xxx >= -31. && xxx < -30.)
                y_corr = -0.05;
        }
        if (fB == -2242 && !fSimu)
        {
            //            if (qqq < 5.8 || qqq > 6.6)
            //                continue;
            if (xxx >= -49. && xxx < -48.)
                y_corr = 0.05;
            if (xxx >= -48. && xxx < -47.)
                y_corr = -0.035;
            if (xxx >= -47. && xxx < -46.)
                y_corr = 0.04;
        }
        if (fB == -1672)
        {

            // zurück ändern
            //            if (!fPairs && (qqq < 1.5 || qqq > 2.5))

            if (!fPairs && (qqq < 7.5 || qqq > 8.5))
            {
                if (debug_in)
                    cout << "Not the right charge!" << endl;
                continue;
            }
            if (fPairs && !(qqq > 1.5 && qqq < 2.5) && !(qqq > 5.5 && qqq < 6.5))
            {
                if (debug_in)
                    cout << "Not the right charge!" << endl;
                continue;
            }
            y_corr = 0.0;
        }

        for (int i = 0; i < n_det; i++)
        {
            xMax[i] = -1000.;
            yMax[i] = -1000.;
            zMax[i] = -1000.;
            qMax[i] = -1000.;
            tMax[i] = -1000.;
        }

        id2 = hitTofd->GetDetId();
        if (hitTofd->GetX() / 100. <= 0.)
        {
            // tof rechts
            if (id2 == 1)
            {
                det2 = tofd1r;
            }
            else if (id2 == 2)
            {
                det2 = tofd2r;
            }
        }
        else
        {
            // tof links
            if (id2 == 1)
            {
                det2 = tofd1l;
            }
            else if (id2 == 2)
            {
                det2 = tofd2l;
            }
        }

        x2[det2] = hitTofd->GetX() / 100.;
        y2[det2] = hitTofd->GetY() / 100. + y_corr;
        if (y2[det] < -0.4 || y2[det] > 0.4)
        {
            if (debug_in)
                cout << "Wrong y-position!" << endl;
            continue;
        }

        // Achtung, nur kurzfristig verwenden
        /*
                if (id2 > 1)
                    continue;
                if ( x2[det2] * 100. < 7. || x2[det2] * 100. > 9.)
                    continue;
        */
        fh_xy_tofd_ac->Fill(x2[det2] * 100. + randx, y2[det2] * 100.);

        //        cout << "Test "<< hitTofd->GetY() / 100. << "  corr  " << y_corr << endl;
        //        if(fB != - 1672) y2[det2] = 0.;
        z2[det2] = 0.;
        q2[det2] = hitTofd->GetEloss(); // / 1.132;
        fh_tofd_charge_ac->Fill(q2[det2]);
        // Achtung, ändern
        if (!fPairs)
            q2[det2] = 8.;
        t2[det2] = hitTofd->GetTime();
        fh_tofd_time_ac->Fill(t2[det2]);

        /*        cout << "ToFD Hit"
                     << " x: " << x2[det2] << " y: " << y2[det2] << " q: " << q2[det2] << " t: " << t2[det2] << " ID "
           << id2
                     << endl;
        */
        // register hits for tracker as long a time is in the coincidence window
        if ((abs(t2[det2] - t1[det1]) < 5.) || first)
        {
            if (debug_in)
                cout << "2 particle within 5 ns   " << first << endl;
            // register point for tracker
            detector[countdet] = det2;
            xdet[countdet] = x2[det2];
            ydet[countdet] = y2[det2];
            zdet[countdet] = z2[det2];
            qdet[countdet] = (int)(q2[det2] + 0.5);
            if (debug_in)
            {
                cout << "registered"
                     << " x: " << x2[det2] << " y: " << y2[det2] << " q: " << q2[det2] << " t: " << t2[det2] << " ID "
                     << id2 << endl;
            }
            if (abs(qdet[countdet] - 2.) < 0.5)
            {
                alpha = true;
                /*                cout << "alpha: " << alpha << endl;
                                cout << "registered"
                                     << " x: " << x2[det2] << " y: " << y2[det2] << " q: " << q2[det2] << " t: " <<
                   t2[det2] << " ID "
                                     << id2 << endl;
                */
                x_4He = xdet[countdet];
                y_4He = ydet[countdet];
                z_4He = zdet[countdet];
            }
            if (abs(qdet[countdet] - 6.) < 0.5)
            {
                carbon = true;
                /*                cout << "carbon: " << carbon << endl;
                                cout << "registered"
                                     << " x: " << x2[det2] << " y: " << y2[det2] << " q: " << q2[det2] << " t: " <<
                   t2[det2] << " ID "
                                     << id2 << endl;
                */
                x_12C = xdet[countdet];
                y_12C = ydet[countdet];
                z_12C = zdet[countdet];
            }
            countdet++;
            single = true;
            first = false;
            tStart = t2[det2];

            det1 = det2;
            x1[det1] = x2[det2];
            y1[det1] = y2[det2];
            z1[det1] = 0.;
            q1[det1] = q2[det2];
            t1[det1] = t2[det2];
            id1 = id2;

            // since we had a coincidence, continue with next event, if not last event.
            if (ihit < nHits - 1)
                continue;
        }

        if (!single)
            continue;

        if (fPairs && !(alpha && carbon))
            continue;
        alpha = false;
        carbon = false;
        //        cout << "Found pair!!!!!! " << endl;

        delta = sqrt((x_12C - x_4He) * (x_12C - x_4He) + (y_12C - y_4He) * (y_12C - y_4He));
        //		cout << "Delta: " << delta << endl;

        if (debug2)
            cout << "found good event, xpos " << x2[det2] * 100. << " counts tofd " << hits1 << endl;
        hits1++;

        single = false;
        first = true;

        if (ihit < nHits - 1)
            ihit--;

        det1 = det2;
        x1[det1] = x2[det2];
        y1[det1] = y2[det2];
        z1[det1] = 0.;
        q1[det1] = q2[det2];
        t1[det1] = t2[det2];
        id1 = id2;

        multTofd++;

        counterTofdMulti++;

        // only consider fiber with maximum energy deposit
        Bool_t maxWerte = false;
        if (!fPairs && fB != -1672)
            maxWerte = true;

        // cut in ToT for Fibers
        Double_t cutQ = 0.;
        if (!fPairs || fB != -1672)
            // check if cut can be applied
            cutQ = 0.;

        if (debug_in)
            cout << "start fiber analysis" << endl;

        // loop over fiber 13
        auto detHit13 = fHitItems.at(DET_FI13);
        Int_t nHits13 = detHit13->GetEntriesFast();
        LOG(debug) << "Fi13 hits: " << nHits13 << endl;

        Int_t mult13 = 0;
        for (Int_t ihit13 = 0; ihit13 < nHits13; ihit13++)
        {
            det = fi13;
            R3BBunchedFiberHitData* hit13 = dynamic_cast<R3BBunchedFiberHitData*>(detHit13->At(ihit13));
            x1[det] = hit13->GetX() / 100.;
            y1[det] = hit13->GetY() / 100.;
            z1[det] = 0.;
            q1[det] = hit13->GetEloss();
            /*
                        if (q1[det] > 9.)
                            q1[det] = 8.;
                        else if (q1[det] > 4.)
                            q1[det] = 6.;
                        else if (q1[det] > 0.)
                            q1[det] = 2.;
            */
            t1[det] = hit13->GetTime();
            tof = tStart - t1[det];

            // Fill histograms before cuts
            fh_Fib_ToF[det]->Fill(x1[det] * 100., tof);
            fh_xy_Fib[det]->Fill(x1[det] * 100., y1[det] * 100.);
            fh_mult_Fib[det]->Fill(nHits13);
            fh_ToT_Fib[det]->Fill(x1[det] * 100., q1[det]);
            fh_Fibs_vs_Tofd[det]->Fill(x1[tofd1r] * 100. + randx, x1[det] * 100.);
            fh_Fib_vs_Events[det]->Fill(fNEvents, x1[det] * 100.);
            fh_ToF_vs_Events[det]->Fill(fNEvents, tof);
            fh_Fib_Time[det]->Fill(x1[det] * 100., t1[det]);

            if (debug2)
                cout << "Fi13 bc: " << ihit13 << " x1: " << x1[det] << " y1: " << y1[det] << " q1: " << q1[det]
                     << " t1: " << t1[det] << endl;

            hits13bc++;

            // Cuts on Fi13
            //            if (fCuts && x1[det] * 100. < -24.5)
            if (fCuts && x1[det] * 100. < -25.75)
                continue;
            if (fCuts && (y1[det] < -0.4 || y1[det] > 0.4))
                continue;
            if (fCuts && !fPairs && q1[det] < cutQ)
                continue;
            if (fCuts && (tof < -20 || tof > 20) && !fSimu)
                continue;
            //            if (fGraphCuts && !cut_Fi13vsTofd->IsInside(x1[tofd1r] * 100., x1[det] * 100.))
            //                continue;

            hits13++;

            xFi13[mult13] = x1[det];
            yFi13[mult13] = y1[det];
            qFi13[mult13] = q1[det];
            tFi13[mult13] = tof;
            mult13++;
            if (mult13 > 100)
                continue;

            if (q1[det] > qMax[det])
            {
                qMax[det] = q1[det];
                xMax[det] = x1[det];
                yMax[det] = y1[det];
                zMax[det] = z1[det];
                tMax[det] = t1[det];
            }

            // Fill histograms
            fh_Fib_ToF_ac[det]->Fill(x1[det] * 100., tof);
            fh_xy_Fib_ac[det]->Fill(x1[det] * 100., y1[det] * 100.);
            fh_mult_Fib_ac[det]->Fill(mult13);
            fh_ToT_Fib_ac[det]->Fill(x1[det] * 100., q1[det]);
            fh_Fibs_vs_Tofd_ac[det]->Fill(x1[tofd1r] * 100. + randx, x1[det] * 100.);
            fh_Fib_vs_Events_ac[det]->Fill(fNEvents, x1[det] * 100.);
            fh_ToF_vs_Events_ac[det]->Fill(fNEvents, tof);
            fh_Fib_Time_ac[det]->Fill(x1[det] * 100., t1[det]);

            if (debug2)
                cout << "Fi13: " << ihit13 << " x1: " << x1[det] << " y1: " << y1[det] << " q1: " << q1[det]
                     << " t1: " << t1[det] << endl;
            /*
                        if (!maxWerte)
                        {
                            detector[countdet] = det;
                            xdet[countdet] = x1[det];
                            ydet[countdet] = y1[det];
                            zdet[countdet] = z1[det];
                            qdet[countdet] = 0;
                            countdet++;
                        }
            */
        }

        if (mult13 > 0 && maxWerte)
        {
            detector[countdet] = fi13;
            xdet[countdet] = xMax[fi13];
            ydet[countdet] = yMax[fi13];
            zdet[countdet] = zMax[fi13];
            qdet[countdet] = 0;
            countdet++;
        }

        // loop over fiber 11
        auto detHit11 = fHitItems.at(DET_FI11);
        Int_t nHits11 = detHit11->GetEntriesFast();
        LOG(debug) << "Fi11 hits: " << nHits11 << endl;
        Int_t mult11 = 0;
        for (Int_t ihit11 = 0; ihit11 < nHits11; ihit11++)
        {
            det = fi11;
            R3BBunchedFiberHitData* hit11 = dynamic_cast<R3BBunchedFiberHitData*>(detHit11->At(ihit11));
            x1[det] = hit11->GetX() / 100.;
            y1[det] = hit11->GetY() / 100.;
            z1[det] = 0.;
            q1[det] = hit11->GetEloss();
            /*
                        if (q1[det] > 9.)
                            q1[det] = 8.;
                        else if (q1[det] > 4.)
                            q1[det] = 6.;
                        else if (q1[det] > 0.)
                            q1[det] = 2.;
            */
            t1[det] = hit11->GetTime();
            tof = tStart - t1[det];

            // Fill histograms before cuts
            fh_Fib_ToF[det]->Fill(x1[det] * 100., tof);
            fh_xy_Fib[det]->Fill(x1[det] * 100., y1[det] * 100);
            fh_mult_Fib[det]->Fill(nHits11);
            fh_ToT_Fib[det]->Fill(x1[det] * 100., q1[det]);
            fh_Fibs_vs_Tofd[det]->Fill(x1[tofd1r] * 100. + randx, x1[det] * 100.);
            fh_Fib_vs_Events[det]->Fill(fNEvents, x1[det] * 100.);
            fh_ToF_vs_Events[det]->Fill(fNEvents, tof);
            fh_Fib_Time[det]->Fill(x1[det] * 100., t1[det]);

            if (debug2)
                cout << "Fi11 bc: " << ihit11 << " x1: " << x1[det] << " y1: " << y1[det] << " q1: " << q1[det]
                     << " t1: " << t1[det] << endl;
            hits11bc++;

            // Cuts on Fi11
            //            if (fCuts && x1[det] * 100. < -24.4)
            if (fCuts && x1[det] * 100. < -25.75)
                continue;
            if (fCuts && !fPairs && q1[det] < cutQ)
                continue;
            if (fCuts && (y1[det] < -0.4 || y1[det] > 0.4))
                continue;
            if (fCuts && (tof < -20 || tof > 20) && !fSimu)
                continue;
            //            if (fGraphCuts && !cut_Fi13vsTofd->IsInside(x1[tofd1r] * 100., x1[det] * 100.))
            //                continue;

            hits11++;

            xFi11[mult11] = x1[det];
            yFi11[mult11] = y1[det];
            qFi11[mult11] = q1[det];
            tFi11[mult11] = tof;
            mult11++;
            if (mult11 > 100)
                continue;

            if (q1[det] > qMax[det])
            {
                qMax[det] = q1[det];
                xMax[det] = x1[det];
                yMax[det] = y1[det];
                zMax[det] = z1[det];
                tMax[det] = t1[det];
            }

            // Fill histograms
            fh_Fib_ToF_ac[det]->Fill(x1[det] * 100., tof);
            fh_xy_Fib_ac[det]->Fill(x1[det] * 100., y1[det] * 100);
            fh_mult_Fib_ac[det]->Fill(mult11);
            fh_ToT_Fib_ac[det]->Fill(x1[det] * 100., q1[det]);
            fh_Fibs_vs_Tofd_ac[det]->Fill(x1[tofd1r] * 100. + randx, x1[det] * 100.);
            fh_Fib_vs_Events_ac[det]->Fill(fNEvents, x1[det] * 100.);
            fh_ToF_vs_Events_ac[det]->Fill(fNEvents, tof);
            fh_Fib_Time_ac[det]->Fill(x1[det] * 100., t1[det]);

            if (debug2)
                cout << "Fi11: " << ihit11 << " x1: " << x1[det] << " y1: " << y1[det] << " q1: " << q1[det]
                     << " t1: " << t1[det] << endl;
            /*
                        if (!maxWerte)
                        {
                            detector[countdet] = det;
                            xdet[countdet] = x1[det];
                            ydet[countdet] = y1[det];
                            zdet[countdet] = z1[det];
                            qdet[countdet] = 0;
                            countdet++;
                        }
            */
        }

        if (mult11 > 0 && maxWerte)
        {
            detector[countdet] = fi11;
            xdet[countdet] = xMax[fi11];
            ydet[countdet] = yMax[fi11];
            zdet[countdet] = zMax[fi11];
            qdet[countdet] = 0;
            countdet++;
        }

        // loop over fiber 10
        auto detHit10 = fHitItems.at(DET_FI10);
        Int_t nHits10 = detHit10->GetEntriesFast();
        LOG(debug) << "Fi10 hits: " << nHits10 << endl;
        Int_t mult10 = 0;
        for (Int_t ihit10 = 0; ihit10 < nHits10; ihit10++)
        {
            det = fi10;
            R3BBunchedFiberHitData* hit10 = dynamic_cast<R3BBunchedFiberHitData*>(detHit10->At(ihit10));
            x1[det] = hit10->GetX() / 100.;
            y1[det] = hit10->GetY() / 100.;
            z1[det] = 0.;
            q1[det] = hit10->GetEloss();
            /*
                        if (q1[det] > 9.)
                            q1[det] = 8.;
                        else if (q1[det] > 4.)
                            q1[det] = 6.;
                        else if (q1[det] > 0.)
                            q1[det] = 2.;
            */
            t1[det] = hit10->GetTime();
            tof = tStart - t1[det];

            // Fill histograms before cuts
            fh_Fib_ToF[det]->Fill(x1[det] * 100., tof);
            fh_xy_Fib[det]->Fill(x1[det] * 100., y1[det] * 100.);
            fh_mult_Fib[det]->Fill(nHits10);
            fh_ToT_Fib[det]->Fill(x1[det] * 100., q1[det]);
            fh_Fibs_vs_Tofd[det]->Fill(x1[tofd1l] * 100. + randx, x1[det] * 100.);
            fh_Fib_vs_Events[det]->Fill(fNEvents, x1[det] * 100.);
            fh_ToF_vs_Events[det]->Fill(fNEvents, tof);
            fh_Fib_Time[det]->Fill(x1[det] * 100., t1[det]);

            hits10bc++;

            // Cuts on Fi10
            if (fCuts && !fPairs && q1[det] < cutQ)
                continue;
            //            if (fCuts && x1[det] * 100. < -24.1)
            if (fCuts && x1[det] * 100. < -22.85)
                continue;
            if (fCuts && (y1[det] < -0.4 || y1[det] > 0.4))
                continue;
            if (fCuts && (tof < -20 || tof > 20) && !fSimu)
                continue;

            hits10++;

            xFi10[mult10] = x1[det];
            yFi10[mult10] = y1[det];
            qFi10[mult10] = q1[det];
            tFi10[mult10] = tof;
            mult10++;
            if (mult10 > 100)
                continue;

            if (q1[det] > qMax[det])
            {
                qMax[det] = q1[det];
                xMax[det] = x1[det];
                yMax[det] = y1[det];
                zMax[det] = z1[det];
                tMax[det] = t1[det];
            }

            // Fill histograms
            fh_Fib_ToF_ac[det]->Fill(x1[det] * 100., tof);
            fh_xy_Fib_ac[det]->Fill(x1[det] * 100., y1[det] * 100.);
            fh_mult_Fib_ac[det]->Fill(mult10);
            fh_ToT_Fib_ac[det]->Fill(x1[det] * 100., q1[det]);
            fh_Fibs_vs_Tofd_ac[det]->Fill(x1[tofd1l] * 100. + randx, x1[det] * 100.);
            fh_Fib_vs_Events_ac[det]->Fill(fNEvents, x1[det] * 100.);
            fh_ToF_vs_Events_ac[det]->Fill(fNEvents, tof);
            fh_Fib_Time_ac[det]->Fill(x1[det] * 100., t1[det]);

            if (debug2)
                cout << "Fi10: " << ihit10 << " x1: " << x1[det] << " y1: " << y1[det] << " q1: " << q1[det]
                     << " t1: " << t1[det] << endl;
            /*
                        if (!maxWerte)
                        {
                            detector[countdet] = det;
                            xdet[countdet] = x1[det];
                            ydet[countdet] = y1[det];
                            zdet[countdet] = z1[det];
                            qdet[countdet] = 0;
                            countdet++;
                        }
            */
        }

        if (mult10 > 0 && maxWerte)
        {
            detector[countdet] = fi10;
            xdet[countdet] = xMax[fi10];
            ydet[countdet] = yMax[fi10];
            zdet[countdet] = zMax[fi10];
            qdet[countdet] = 0;
            countdet++;
        }

        // loop over fiber 12
        auto detHit12 = fHitItems.at(DET_FI12);
        Int_t nHits12 = detHit12->GetEntriesFast();
        LOG(debug) << "Fi12 hits: " << nHits12 << endl;
        Int_t mult12 = 0;
        for (Int_t ihit12 = 0; ihit12 < nHits12; ihit12++)
        {
            det = fi12;
            R3BBunchedFiberHitData* hit12 = dynamic_cast<R3BBunchedFiberHitData*>(detHit12->At(ihit12));
            x1[det] = hit12->GetX() / 100.;
            y1[det] = hit12->GetY() / 100.;
            z1[det] = 0.;
            q1[det] = hit12->GetEloss();
            /*
                        if (q1[det] > 9.)
                            q1[det] = 8.;
                        else if (q1[det] > 4.)
                            q1[det] = 6.;
                        else if (q1[det] > 0.)
                            q1[det] = 2.;
            */
            t1[det] = hit12->GetTime();
            tof = tStart - t1[det];

            // Fill histograms before cuts
            fh_Fib_ToF[det]->Fill(x1[det] * 100., tof);
            fh_xy_Fib[det]->Fill(x1[det] * 100., y1[det] * 100.);
            fh_mult_Fib[det]->Fill(nHits12);
            fh_ToT_Fib[det]->Fill(x1[det] * 100., q1[det]);
            fh_Fibs_vs_Tofd[det]->Fill(x1[tofd1l] * 100. + randx, x1[det] * 100.);
            fh_Fib_vs_Events[det]->Fill(fNEvents, x1[det] * 100.);
            fh_ToF_vs_Events[det]->Fill(fNEvents, tof);
            fh_Fib_Time[det]->Fill(x1[det] * 100., t1[det]);

            hits12bc++;

            // Cuts on Fi12
            if (fCuts && !fPairs && q1[det] < cutQ)
                continue;
            //            if (fCuts && x1[det] * 100. < -24.4)
            if (fCuts && x1[det] * 100. < -23.45)
                continue;
            if (fCuts && (y1[det] < -0.4 || y1[det] > 0.4))
                continue;
            if (fCuts && (tof < -20 || tof > 20) && !fSimu)
                continue;

            hits12++;

            xFi12[mult12] = x1[det];
            yFi12[mult12] = y1[det];
            qFi12[mult12] = q1[det];
            tFi12[mult12] = tof;
            mult12++;
            if (mult12 > 100)
                continue;

            if (q1[det] > qMax[det])
            {
                qMax[det] = q1[det];
                xMax[det] = x1[det];
                yMax[det] = y1[det];
                zMax[det] = z1[det];
                tMax[det] = t1[det];
            }

            // Fill histograms
            fh_Fib_ToF_ac[det]->Fill(x1[det] * 100., tof);
            fh_xy_Fib_ac[det]->Fill(x1[det] * 100., y1[det] * 100.);
            fh_mult_Fib_ac[det]->Fill(mult12);
            fh_ToT_Fib_ac[det]->Fill(x1[det] * 100., q1[det]);
            fh_Fibs_vs_Tofd_ac[det]->Fill(x1[tofd1l] * 100. + randx, x1[det] * 100.);
            fh_Fib_vs_Events_ac[det]->Fill(fNEvents, x1[det] * 100.);
            fh_ToF_vs_Events_ac[det]->Fill(fNEvents, tof);
            fh_Fib_Time_ac[det]->Fill(x1[det] * 100., t1[det]);

            if (debug2)
                cout << "Fi12: " << ihit12 << " x1: " << x1[det] << " y1: " << y1[det] << " q1: " << q1[det]
                     << " t1: " << t1[det] << endl;
            /*
                        if (!maxWerte)
                        {
                            detector[countdet] = det;
                            xdet[countdet] = x1[det];
                            ydet[countdet] = y1[det];
                            zdet[countdet] = z1[det];
                            qdet[countdet] = 0;
                            countdet++;
                        }
            */
        }

        if (mult12 > 0 && maxWerte)
        {
            detector[countdet] = fi12;
            xdet[countdet] = xMax[fi12];
            ydet[countdet] = yMax[fi12];
            zdet[countdet] = zMax[fi12];
            qdet[countdet] = 0;
            countdet++;
        }

        // loop over fiber 3a
        auto detHit3a = fHitItems.at(DET_FI3A);
        Int_t nHits3a = detHit3a->GetEntriesFast();
        LOG(debug) << "Fi3a hits: " << nHits3a << endl;
        Int_t mult3a = 0;
        fh_mult_Fib[fi3a]->Fill(nHits3a);
        for (Int_t ihit3a = 0; ihit3a < nHits3a; ihit3a++)
        {
            det = fi3a;
            R3BBunchedFiberHitData* hit3a = dynamic_cast<R3BBunchedFiberHitData*>(detHit3a->At(ihit3a));
            x1[det] = hit3a->GetX() / 100.;
            y1[det] = hit3a->GetY() / 100.;
            z1[det] = 0.;
            q1[det] = hit3a->GetEloss();
            /*
                        if (q1[det] > 5.)
                            q1[det] = 8.;
                        else if (q1[det] > 2.)
                            q1[det] = 6.;
                        else if (q1[det] > 0.)
                            q1[det] = 2.;
            */
            t1[det] = hit3a->GetTime();
            tof = tStart - t1[det];

            // Fill histograms before cuts
            fh_Fib_ToF[det]->Fill(x1[det] * 100., tof);
            fh_xy_Fib[det]->Fill(x1[det] * 100., y1[det] * 100.);
            fh_mult_Fib[det]->Fill(nHits3a);
            fh_ToT_Fib[det]->Fill(x1[det] * 100., q1[det]);
            fh_Fibs_vs_Tofd[det]->Fill(x1[tofd1r] * 100. + randx, x1[det] * 100.);
            fh_Fib_vs_Events[det]->Fill(fNEvents, x1[det] * 100.);
            fh_ToF_vs_Events[det]->Fill(fNEvents, tof);
            fh_Fib_Time[det]->Fill(x1[det] * 100., t1[det]);

            // Cuts on Fi3a
            //            if (fCuts && y1[det] * 100. > 50.)
            //                continue;
            //            if (fCuts && y1[det] *100. < -50.)
            //                continue;
            // if (fCuts && !fPairs && q1[det] < cutQ)
            // continue;
            if (fCuts && (y1[det] < -0.2 || y1[det] > 0.2))
                continue;
            if (fCuts && (tof < -20 || tof > 20) && !fSimu)
                continue;

            xFi3a[mult3a] = x1[det];
            yFi3a[mult3a] = y1[det];
            qFi3a[mult3a] = q1[det];
            tFi3a[mult3a] = tof;
            mult3a++;
            if (mult3a > 100)
                continue;

            if (q1[det] > qMax[det])
            {
                qMax[det] = q1[det];
                xMax[det] = x1[det];
                yMax[det] = y1[det];
                zMax[det] = z1[det];
                tMax[det] = t1[det];
            }

            // Fill histograms
            fh_Fib_ToF_ac[det]->Fill(x1[det] * 100., tof);
            fh_xy_Fib_ac[det]->Fill(x1[det] * 100., y1[det] * 100.);
            fh_mult_Fib_ac[det]->Fill(mult3a);
            fh_ToT_Fib_ac[det]->Fill(x1[det] * 100., q1[det]);
            fh_Fibs_vs_Tofd_ac[det]->Fill(x1[tofd1r] * 100. + randx, x1[det] * 100.);
            fh_Fib_vs_Events_ac[det]->Fill(fNEvents, x1[det] * 100.);
            fh_ToF_vs_Events_ac[det]->Fill(fNEvents, tof);
            fh_Fib_Time_ac[det]->Fill(x1[det] * 100., t1[det]);

            if (debug2)
                cout << "Fi3a " << ihit3a << " x1: " << x1[det] << " y1: " << y1[det] << " q1: " << q1[det]
                     << " t1: " << tof << endl;
            /*
                        if (!maxWerte)
                        {
                            detector[countdet] = det;
                            xdet[countdet] = x1[det];
                            ydet[countdet] = y1[det];
                            zdet[countdet] = z1[det];
                            qdet[countdet] = 0;
                            countdet++;
                        }
            */
        }

        if (mult3a > 0 && maxWerte)
        {
            detector[countdet] = fi3a;
            xdet[countdet] = xMax[fi3a];
            ydet[countdet] = yMax[fi3a];
            zdet[countdet] = zMax[fi3a];
            qdet[countdet] = 0;
            countdet++;
        }

        // loop over fiber 3b
        auto detHit3b = fHitItems.at(DET_FI3B);
        Int_t nHits3b = detHit3b->GetEntriesFast();
        LOG(debug) << "Fi3b hits: " << nHits3b << endl;
        Int_t mult3b = 0;
        fh_mult_Fib[fi3b]->Fill(nHits3b);

        for (Int_t ihit3b = 0; ihit3b < nHits3b; ihit3b++)
        {
            det = fi3b;
            R3BBunchedFiberHitData* hit3b = dynamic_cast<R3BBunchedFiberHitData*>(detHit3b->At(ihit3b));
            x1[det] = hit3b->GetX() / 100.;
            y1[det] = hit3b->GetY() / 100.;
            z1[det] = 0.;
            q1[det] = hit3b->GetEloss();
            /*
                        if (q1[det] > 5.)
                            q1[det] = 8.;
                        else if (q1[det] > 2.)
                            q1[det] = 6.;
                        else if (q1[det] > 0.)
                            q1[det] = 2.;
            */
            t1[det] = hit3b->GetTime();
            tof = tStart - t1[det];

            // Fill histograms before cuts
            fh_Fib_ToF[det]->Fill(x1[det] * 100., tof);
            fh_xy_Fib[det]->Fill(x1[det] * 100., y1[det] * 100.);
            fh_mult_Fib[det]->Fill(nHits3b);
            fh_ToT_Fib[det]->Fill(x1[det] * 100., q1[det]);
            fh_Fibs_vs_Tofd[det]->Fill(x1[tofd1l] * 100. + randx, x1[det] * 100.);
            fh_Fib_vs_Events[det]->Fill(fNEvents, x1[det] * 100.);
            fh_ToF_vs_Events[det]->Fill(fNEvents, tof);
            fh_Fib_Time[det]->Fill(x1[det] * 100., t1[det]);

            // Cuts on Fi3b
            //            if (fCuts && y1[det] * 100. > 50.)
            //                continue;
            //            if (fCuts && y1[det] * 100. < -50.)
            //                continue;
            // if (fCuts && !fPairs && q1[det] < cutQ)
            // continue;
            if (fCuts && (y1[det] < -0.2 || y1[det] > 0.2))
                continue;
            if (fCuts && (tof < -20 || tof > 20) && !fSimu)
                continue;

            xFi3b[mult3b] = x1[det];
            yFi3b[mult3b] = y1[det];
            qFi3b[mult3b] = q1[det];
            tFi3b[mult3b] = tof;
            mult3b++;
            if (mult3b > 100)
                continue;

            if (q1[det] > qMax[det])
            {
                qMax[det] = q1[det];
                xMax[det] = x1[det];
                yMax[det] = y1[det];
                zMax[det] = z1[det];
                tMax[det] = t1[det];
            }

            // Fill histograms
            fh_Fib_ToF_ac[det]->Fill(x1[det] * 100., tof);
            fh_xy_Fib_ac[det]->Fill(x1[det] * 100., y1[det] * 100.);
            fh_mult_Fib_ac[det]->Fill(mult3b);
            fh_ToT_Fib_ac[det]->Fill(x1[det] * 100., q1[det]);
            fh_Fibs_vs_Tofd_ac[det]->Fill(x1[tofd1l] * 100. + randx, x1[det] * 100.);
            fh_Fib_vs_Events_ac[det]->Fill(fNEvents, x1[det] * 100.);
            fh_ToF_vs_Events_ac[det]->Fill(fNEvents, tof);
            fh_Fib_Time_ac[det]->Fill(x1[det] * 100., t1[det]);

            if (debug2)
                cout << "Fi3b " << ihit3b << " x1: " << x1[det] << " y1: " << y1[det] << " q1: " << q1[det]
                     << " t1: " << tof << endl;
            /*
                        if (!maxWerte)
                        {
                            detector[countdet] = det;
                            xdet[countdet] = x1[det];
                            ydet[countdet] = y1[det];
                            zdet[countdet] = z1[det];
                            qdet[countdet] = 0;
                            countdet++;
                        }
            */
        }

        if (mult3b > 0 && maxWerte)
        {
            detector[countdet] = fi3b;
            xdet[countdet] = xMax[fi3b];
            ydet[countdet] = yMax[fi3b];
            zdet[countdet] = zMax[fi3b];
            qdet[countdet] = 0;
            countdet++;
        }

        if (fGhost && (mult3a + mult3b) < 1)
        {
            det = ghost;
            detector[countdet] = det;
            xdet[countdet] = 0.;
            ydet[countdet] = 0.;
            zdet[countdet] = 0.;
            qdet[countdet] = 0;
            countdet++;
        }

        // Plots of correlations of Fiber detectors and register events for tracker
        for (Int_t i = 0; i < mult11; i++)
        {
            for (Int_t j = 0; j < mult13; j++)
            {
                if (xFi11[i] > -100 && xFi13[j] > -100)
                {
                    fh_Fib13_vs_Fib11->Fill(xFi11[i] * 100., xFi13[j] * 100.);
                    fh_Fib13_vs_Fib11_dx->Fill(xFi11[i] * 100., xFi13[j] * 100. - xFi11[i] * 100.);
                    if (!maxWerte)
                    {
                        if (!fFi11[i])
                        {
                            detector[countdet] = fi11;
                            xdet[countdet] = xFi11[i];
                            ydet[countdet] = yFi11[i];
                            zdet[countdet] = 0.;
                            qdet[countdet] = 0; // qFi11[i];
                            countdet++;
                            fFi11[i] = true;
                        }
                        if (!fFi13[j])
                        {
                            detector[countdet] = fi13;
                            xdet[countdet] = xFi13[j];
                            ydet[countdet] = yFi13[j];
                            zdet[countdet] = 0.;
                            qdet[countdet] = 0; // qFi13[j];
                            countdet++;
                            fFi13[j] = true;
                        }
                    }
                }
            }
            for (Int_t j = 0; j < mult3a; j++)
            {
                if (xFi3a[j] > -100 && xFi11[i] > -100)
                {
                    if (fGraphCuts && !cut_fi11_fi3a->IsInside(xFi3a[j] * 100., xFi11[i] * 100.))
                        continue;
                    fh_Fib11_vs_Fib3a->Fill(xFi3a[j] * 100., xFi11[i] * 100.);
                    fh_Fib11_vs_Fib3a_dx->Fill(xFi3a[j] * 100., xFi11[i] * 100. - xFi3a[j] * 100.);
                    if (!maxWerte)
                    {
                        if (!fFi3a[j])
                        {
                            detector[countdet] = fi3a;
                            xdet[countdet] = xFi3a[j];
                            ydet[countdet] = yFi3a[j];
                            zdet[countdet] = 0.;
                            qdet[countdet] = 0; // qFi3a[j];
                            countdet++;
                            fFi3a[j] = true;
                        }
                    }
                }
            }
        }

        for (Int_t i = 0; i < mult12; i++)
        {
            for (Int_t j = 0; j < mult10; j++)
            {
                if (xFi12[i] > -100 && xFi10[j] > -100)
                {
                    fh_Fib10_vs_Fib12->Fill(xFi12[i] * 100., xFi10[j] * 100.);
                    fh_Fib10_vs_Fib12_dx->Fill(xFi12[i] * 100., xFi10[j] * 100. - xFi12[i] * 100.);
                    if (!maxWerte)
                    {
                        if (!fFi12[i])
                        {
                            detector[countdet] = fi12;
                            xdet[countdet] = xFi12[i];
                            ydet[countdet] = yFi12[i];
                            zdet[countdet] = 0.;
                            qdet[countdet] = 0; // qFi12[i];
                            countdet++;
                            fFi12[i] = true;
                        }
                        if (!fFi10[j])
                        {
                            detector[countdet] = fi10;
                            xdet[countdet] = xFi10[j];
                            ydet[countdet] = yFi10[j];
                            zdet[countdet] = 0.;
                            qdet[countdet] = 0; // qFi10[j];
                            countdet++;
                            fFi10[j] = true;
                        }
                    }
                }
            }
            for (Int_t j = 0; j < mult3b; j++)
            {
                if (xFi3b[j] > -100 && xFi12[i] > -100)
                {
                    if (fGraphCuts && !cut_fi12_fi3b->IsInside(xFi3b[j] * 100., xFi12[i] * 100.))
                        continue;
                    fh_Fib12_vs_Fib3b->Fill(xFi3b[j] * 100., xFi12[i] * 100.);
                    fh_Fib12_vs_Fib3b_dx->Fill(xFi3b[j] * 100., xFi12[i] * 100. - xFi3b[j] * 100.);
                    if (!maxWerte)
                    {
                        if (!fFi3b[j])
                        {
                            detector[countdet] = fi3b;
                            xdet[countdet] = xFi3b[j];
                            ydet[countdet] = yFi3b[j];
                            zdet[countdet] = 0.;
                            qdet[countdet] = 0; // qFi3b[j];
                            countdet++;
                            fFi3b[j] = true;
                        }
                    }
                }
            }
        }

        if (countdet > 50)
        {
            cout << "Too many hits!!!" << endl;
            continue;
        }

        // here call tracker
        chi2 = 1.E100;
        Bool_t debug = false;
        if (debug && ((mult10 > 0 && mult12 > 0 && mult3b > 0) || (mult11 > 0 && mult13 > 0 && mult3a > 0)))
        {
            cout << "# of points" << countdet << endl;
            for (Int_t i = 0; i < countdet; i++)
            {
                cout << "#" << i << " Det: " << detector[i] << " x: " << xdet[i] << " y: " << ydet[i]
                     << " q: " << qdet[i] << endl;
                if (detector[i] < 6)
                    qdet[i] = 0;
            }

            for (Int_t i = 0; i < ndet; i++)
            {
                LOG(debug2) << "Max Det: " << i << " max x: " << xMax[i] << " max y: " << yMax[i]
                            << " max q: " << qMax[i] << endl;
            }
        }

        for (Int_t i = 0; i < countdet; i++)
        {
            if (detector[i] < 6)
                qdet[i] = 0;
        }

        Bool_t writeFile = false;
        if (writeFile && counter1 < 200000 &&
            ((mult10 > 0 && mult12 > 0 && mult3b > 0) || (mult11 > 0 && mult13 > 0 && mult3a > 0)))
        {
            ofstream myfile;
            myfile.open("events_v35.dat", ios::out | ios::app);
            counter1++;
            cout << "Counter: " << counter1 << endl;
            myfile << countdet << "  \n";
            for (Int_t i = 0; i < countdet; i++)
            {

                myfile << detector[i] << "    " << xdet[i] << "    " << ydet[i] << " \n";
            }
            myfile.close();
        }
        // here call tracker
        chi2 = 1.E100;
        counter2++;
        if (tracker && fPairs && ((mult10 > 0 && mult12 > 0 && mult3b > 0) || (mult11 > 0 && mult13 > 0 && mult3a > 0)))
        {
            // double track
            counter2++;
            qdet[max] = 2;
            qdet[max + 1] = 6;
            Bool_t det_coord = true;
            Bool_t st = true;
            multi_track_extended_output_from_cpp_(
                &max, &countdet, &det_coord, &st, target, detector, qdet, xdet, ydet, zdet, track, chi, pat1, pat2);

            chi2 = chi[4] + chi[5];
            fh_chiy_vs_chix->Fill(chi[0], chi[1]);
            fh_chiy_vs_chix->Fill(chi[2], chi[3]);
            fh_chi2->Fill(chi2);

            if (debug)
            {
                cout << "track1: " << track[0] << "  " << track[1] << "  " << track[2] << "  " << track[3] << "  "
                     << track[4] << "  " << track[5] << endl;
                cout << "track2: " << track[6] << "  " << track[7] << "  " << track[8] << "  " << track[9] << "  "
                     << track[10] << "  " << track[11] << endl;

                cout << "chi: " << chi[0] << "  " << chi[1] << "  " << chi[2] << "  " << chi[3] << "  " << chi[4]
                     << "  " << chi[5] << endl;

                cout << "******************************************" << endl;
                cout << "Track In 4He"
                     << "x " << XHes << " y " << YHes << " z " << ZHes << endl;
                cout << "px " << pHexs << " py " << pHeys << " z " << pHezs << endl;

                cout << "Track In 12C"
                     << "x " << XCs << " y " << YCs << " z " << ZCs << endl;
                cout << "px " << pCxs << " py " << pCys << " z " << pCzs << endl;
            }

            if (chi[0] < 1e10)
                counter3++;
            if (chi[1] < 1e10)
                counter4++;

            if (chi[0] < 1e10 && chi[1] < 1e10)
            {
                counterTracker++;
                Output1(track, chi);
            }
            // we have a hit
            for (Int_t i = 0; i < ndet; i++)
            {
                xTrack[i] = -1000.;
                yTrack[i] = -1000.;
                zTrack[i] = -1000.;
                qTrack[i] = -1000.;
            }
            Int_t charge = 0;
            if (debug)
                cout << "# of points back" << countdet << endl;
            for (Int_t i = 0; i < countdet; i++)
            {

                if (debug)
                {
                    cout << "back #" << i << " Det: " << detector[i] << " x: " << xdet[i] << " y: " << ydet[i]
                         << " q: " << qdet[i] << endl;
                }
                if (qdet[i] == 2)
                {
                    xTrack[detector[i]] = xdet[i];
                    yTrack[detector[i]] = ydet[i];
                    zTrack[detector[i]] = zdet[i];
                    qTrack[detector[i]] = qdet[i];
                }
                // plot hits of the track
                if (detector[i] != ghost)
                {
                    fh_xy[detector[i]]->Fill(xdet[i] * 100., ydet[i] * 100.);
                    fh_p_vs_x[detector[i]]->Fill(xdet[i] * 100. + randx, track[5]);
                    fh_p_vs_x_test[detector[i]]->Fill(
                        xdet[i] * 100., sqrt(track[3] * track[3] + track[4] * track[4] + track[5] * track[5]));
                }
            }
            // Plots of correlations of Fiber detectors

            fh_Fib13_vs_Fib11_back->Fill(xTrack[3] * 100., xTrack[5] * 100.);
            fh_Fib13_vs_Fib11_dx_back->Fill(xTrack[3] * 100., xTrack[5] * 100. - xTrack[3] * 100.);
            fh_Fib11_vs_Fib3a_back->Fill(xTrack[0] * 100., xTrack[3] * 100.);
            fh_Fib11_vs_Fib3a_dx_back->Fill(xTrack[0] * 100., xTrack[3] * 100. - xTrack[0] * 100.);
            fh_Fib10_vs_Fib12_back->Fill(xTrack[4] * 100., xTrack[2] * 100.);
            fh_Fib10_vs_Fib12_dx_back->Fill(xTrack[4] * 100., xTrack[2] * 100. - xTrack[4] * 100.);
            fh_Fib12_vs_Fib3b_back->Fill(xTrack[1] * 100., xTrack[4] * 100.);
            fh_Fib12_vs_Fib3b_dx_back->Fill(xTrack[1] * 100., xTrack[4] * 100. - xTrack[1] * 100.);

            if (chi[0] < 1e10 && chi[1] < 1e10)
            {
                // sorted, lowest Charge first
                new ((*fTrackItems)[fNofTrackItems++])
                    R3BTrack(track[0], track[1], track[2], track[3], track[4], track[5], 2, 2, chi[0], chi[1], 0);
                new ((*fTrackItems)[fNofTrackItems++])
                    R3BTrack(track[6], track[7], track[8], track[9], track[10], track[11], 6, 2, chi[2], chi[3], 0);
            }
        }
        if (tracker && !fPairs && ((mult10 > 0 && mult12 > 0) || (mult11 > 0 && mult13 > 0)))
        {
            // single track
            qdet[max] = 8;
            qdet[max + 1] = 8;
            counter2++;
            Bool_t det_coord = true;
            Bool_t st = false;
            // multi_track_from_cpp_(
            //    &max, &countdet, &det_coord, &st, target, detector, qdet, xdet, ydet, zdet, track, chi);
            multi_track_extended_output_from_cpp_(
                &max, &countdet, &det_coord, &st, target, detector, qdet, xdet, ydet, zdet, track, chi, pat1, pat2);

            cout << "back from tracker!" << endl;

            chi2 = chi[0] + chi[1];
            fh_chiy_vs_chix->Fill(chi[0], chi[1]);
            fh_chi2->Fill(chi2);

            if (chi[0] < 1.e10)
                counter3++;
            if (chi[1] < 1.e10)
                counter4++;
            if (chi[0] < 1.e10 && chi[1] < 1.e10)
            {
                // fill histograms
                Output2(track, chi);
            }

            if (debug)
            {
                cout << "track1: " << track[0] << "  " << track[1] << "  " << track[2] << endl;
                cout << "track1: " << track[3] << "  " << track[4] << "  " << track[5] << endl;
                cout << "chi: " << chi[0] << "  " << chi[1] << endl;
            }

            if (chi[0] < 1.e10 && chi[1] < 1.e10)
            {
                counterTracker++;
                // we have a hit
                for (Int_t i = 0; i < ndet; i++)
                {
                    xTrack[i] = -1000.;
                    yTrack[i] = -1000.;
                    zTrack[i] = -1000.;
                    qTrack[i] = -1000.;
                }
                Int_t charge = 0;
                LOG(debug2) << "# of points back" << countdet << endl;
                for (Int_t i = 0; i < countdet; i++)
                {

                    // cout << "back #" << i << " Det: " << detector[i] << " x: " << xdet[i]
                    //     << " y: " << ydet[i] << " q: " << qdet[i] << endl;
                    xTrack[detector[i]] = xdet[i];
                    yTrack[detector[i]] = ydet[i];
                    zTrack[detector[i]] = zdet[i];
                    qTrack[detector[i]] = qdet[i];
                    if (qdet[i] > charge)
                        charge = qdet[i];
                }
                // plot hits of the track
                for (Int_t i = 0; i < ndet; i++)
                {
                    fh_xy[i]->Fill(xTrack[i] * 100., yTrack[i] * 100.);
                    fh_p_vs_x[i]->Fill(xTrack[i] * 100., track[5]);
                    fh_p_vs_x_test[i]->Fill(xTrack[i] * 100.,
                                            sqrt(track[3] * track[3] + track[4] * track[4] + track[5] * track[5]));
                }
                // Plots of correlations of Fiber detectors

                fh_Fib13_vs_Fib11_back->Fill(xTrack[3] * 100., xTrack[5] * 100.);
                fh_Fib13_vs_Fib11_dx_back->Fill(xTrack[3] * 100., xTrack[5] * 100. - xTrack[3] * 100.);
                fh_Fib11_vs_Fib3a_back->Fill(xTrack[0] * 100., xTrack[3] * 100.);
                fh_Fib11_vs_Fib3a_dx_back->Fill(xTrack[0] * 100., xTrack[3] * 100. - xTrack[0] * 100.);
                fh_Fib10_vs_Fib12_back->Fill(xTrack[4] * 100., xTrack[2] * 100.);
                fh_Fib10_vs_Fib12_dx_back->Fill(xTrack[4] * 100., xTrack[2] * 100. - xTrack[4] * 100.);
                fh_Fib12_vs_Fib3b_back->Fill(xTrack[1] * 100., xTrack[4] * 100.);
                fh_Fib12_vs_Fib3b_dx_back->Fill(xTrack[1] * 100., xTrack[4] * 100. - xTrack[1] * 100.);

                // store hits in track level
                new ((*fTrackItems)[fNofTrackItems++])
                    R3BTrack(track[0], track[1], track[2], track[3], track[4], track[5], charge, 2, chi[0], chi[1], 0);
            }
        }
        chi2_best = 1.E10;
        pHex = 0.;
        pHey = 0.;
        pHez = 0.;
        pCx = 0.;
        pCy = 0.;
        pCz = 0.;

        if (fMCTrack)
        {
            Double_t costh26s =
                (pHexs * pCxs + pHeys * pCys + pHezs * pCzs) /
                (sqrt(pHexs * pHexs + pHeys * pHeys + pHezs * pHezs) * sqrt(pCxs * pCxs + pCys * pCys + pCzs * pCzs));

            Double_t theta_26s = acos(costh26s) * 180. / 3.14159; // opening angle

            fh_theta26_simu->Fill(theta_26s);

            Double_t eHes = sqrt(pow(pHes, 2) + pow(mHe, 2)) - mHe;
            Double_t eCs = sqrt(pow(pCs, 2) + pow(mC, 2)) - mC;

            Double_t Erelas = sqrt(pow((mHe + mC + eHes + eCs), 2) - pow(pHes, 2) - pow(pCs, 2) -
                                   2 * pHes * pCs * cos(theta_26s * 3.1415 / 180.)) -
                              (mHe + mC); // Erel

            Double_t m_invs =
                sqrt(mHe * mHe + mC * mC + 2. * sqrt(pCs * pCs + mC * mC) * sqrt(pHes * pHes + mHe * mHe) -
                     2. * pHes * pCs * cos(theta_26s * 3.1415 / 180.));
            Double_t Erelbs = m_invs - mHe - mC;

            fh_Erel_simu->Fill(Erelbs);
            LOG(debug) << "Theta 26 simu: " << theta_26s << " Erel simu: " << Erelas << " " << Erelbs << endl;
        }

        for (Int_t i = 0; i < 12; i++)
        {
            track[i] = 0.;
        }
        for (int i = 0; i < n_det; i++)
        {
            x[i] = 0.;
            y[i] = 0.;
            z[i] = 0.;
            q[i] = 0.;
            t[i] = -1000.;

            x1[i] = 0.;
            y1[i] = 0.;
            z1[i] = 0.;
            q1[i] = 0.;
            t1[i] = -1000.;

            x2[i] = 0.;
            y2[i] = 0.;
            z2[i] = 0.;
            q2[i] = 0.;
            t2[i] = -1000.;

            xMax[i] = -1000.;
            yMax[i] = -1000.;
            zMax[i] = -1000.;
            qMax[i] = -1000.;
            tMax[i] = -1000.;
        }

    } // end ToFD loop

    if (multTofd > 0)
        fh_tofd_mult_ac->Fill(multTofd);
}
void R3BTrackS454::Output1(Double_t track[12], Double_t chi[6])
{

    pHex = track[3];
    pHey = track[4];
    pHez = track[5];

    pCx = track[9];
    pCy = track[10];
    pCz = track[11];

    LOG(debug) << "He: " << pHex << "  " << pHey << "  " << pHez << endl;
    LOG(debug) << "C: " << pCx << "  " << pCy << "  " << pCz << endl;

    fh_target_xy->Fill(track[0] * 100., track[1] * 100.);
    fh_px_He->Fill(pHex);
    fh_py_He->Fill(pHey);
    fh_pz_He->Fill(pHez);

    fh_px_C->Fill(pCx);
    fh_py_C->Fill(pCy);
    fh_pz_C->Fill(pCz);

    // Calculate angle between alphs and C
    if (pCz == 0 || pHez == 0)
        return;

    Double_t costh26 = (pHex * pCx + pHey * pCy + pHez * pCz) /
                       (sqrt(pHex * pHex + pHey * pHey + pHez * pHez) * sqrt(pCx * pCx + pCy * pCy + pCz * pCz));

    Double_t theta_26 = acos(costh26) * 180. / 3.14159; // opening angle

    fh_theta26->Fill(theta_26);

    Double_t pHe = sqrt(pow(pHex, 2) + pow(pHey, 2) + pow(pHez, 2));
    Double_t pC = sqrt(pow(pCx, 2) + pow(pCy, 2) + pow(pCz, 2));
    Double_t eHe = sqrt(pow(pHe, 2) + pow(mHe, 2)) - mHe;
    Double_t eC = sqrt(pow(pC, 2) + pow(mC, 2)) - mC;

    fh_p_He->Fill(pHe);
    fh_p_C->Fill(pC);

    Double_t Erela =
        sqrt(pow((mHe + mC + eHe + eC), 2) - pow(pHe, 2) - pow(pC, 2) - 2 * pHe * pC * cos(theta_26 * 3.1415 / 180.)) -
        (mHe + mC); // Erel

    Double_t m_inv = sqrt(mHe * mHe + mC * mC + 2. * sqrt(pC * pC + mC * mC) * sqrt(pHe * pHe + mHe * mHe) -
                          2. * pHe * pC * cos(theta_26 * 3.1415 / 180.));
    Double_t Erelb = m_inv - mHe - mC;

    fh_Erel->Fill(Erelb);
    fh_dErel_vs_x->Fill(delta * 100., Erela);
    //	fh_dErel_vs_y->Fill(yTest * 100., Erela - 1.);
    cout << "Theta 26: " << theta_26 << " Erel: " << Erela << " " << Erelb << endl;
}

void R3BTrackS454::Output2(Double_t track_parameter[12], Double_t chi_single_parameter[6])
{
    // compare

    LOG(debug) << "******************************************" << endl;
    LOG(debug) << "chi_sqingle_parameter " << chi_single_parameter[0] << "  " << chi_single_parameter[1] << endl;
    LOG(debug) << "xyz: " << track_parameter[0] * 100. << "  " << track_parameter[1] * 100. << "  "
               << track_parameter[2] * 100. << endl;
    LOG(debug) << "p: " << track_parameter[3] << "  " << track_parameter[4] << "  " << track_parameter[5] << endl;
    Double_t p_tot = sqrt(track_parameter[3] * track_parameter[3] + track_parameter[4] * track_parameter[4] +
                          track_parameter[5] * track_parameter[5]);

    fh_target_xy->Fill(track_parameter[0] * 100., track_parameter[1] * 100.);
    fh_target_px->Fill(track_parameter[3]);
    fh_target_py->Fill(track_parameter[4]);
    fh_target_pz->Fill(track_parameter[5]);
    fh_target_p->Fill(p_tot);
    fh_chi2->Fill(chi_single_parameter[0]);

    // Comparison with simulation
    Double_t thetax_simu = Pxf / Pzf * 1000.;
    Double_t thetay_simu = Pyf / Pzf * 1000.;

    fh_thetax_dpx->Fill((Pxf - track_parameter[3]) / Pxf * 100., thetax_simu);
    fh_thetax_dpx_abs->Fill((Pxf - track_parameter[3]), thetax_simu);
    fh_thetay_dpy->Fill((Pyf - track_parameter[4]) / Pyf * 100., thetay_simu);

    fh_thetax_dpy->Fill((Pyf - track_parameter[4]) / Pyf * 100., thetax_simu);
    fh_thetay_dpx->Fill((Pxf - track_parameter[3]) / Pxf * 100., thetay_simu);
    fh_dpy_dpx->Fill((Pxf - track_parameter[3]) / Pxf * 100., (Pyf - track_parameter[4]) / Pyf * 100.);

    fh_x_dpx->Fill((Pxf - track_parameter[3]) / Pxf * 100., Xf);
    fh_y_dpy->Fill((Pyf - track_parameter[4]) / Pyf * 100., Yf);
    fh_dx->Fill(Xf - track_parameter[0] * 100.);
    fh_dy->Fill(Yf - track_parameter[1] * 100.);
    fh_dz->Fill(Zf - track_parameter[2] * 100.);
    fh_dpx->Fill((Pxf - track_parameter[3]) / Pxf * 100.);
    fh_dpy->Fill((Pyf - track_parameter[4]) / Pyf * 100.);
    fh_dpz->Fill((Pzf - track_parameter[5]) / Pzf * 100.);
    fh_dp->Fill((Pf_tot - p_tot) / Pf_tot * 100.);
}
void R3BTrackS454::FinishEvent()
{

    fTrackItems->Clear();
    fNofTrackItems = 0;

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

void R3BTrackS454::FinishTask()
{

    cout << "Statistics:" << endl;
    cout << "Events: " << fNEvents << endl;
    cout << "Wrong Trigger: " << counterWrongTrigger << endl;
    cout << "Wrong Tpat: " << counterWrongTpat << endl;
    cout << "ROLU veto: " << counterRolu << endl;
    cout << "Califa veto: " << counterCalifa << endl;
    cout << "TofD: " << counterTofd << endl;
    cout << "TofD multi: " << counterTofdMulti << endl;
    cout << "Tracker: " << counterTracker << endl;

    cout << "Hits TofD " << hits1 << endl;
    cout << "Eff. Fi10 min: " << hits10 << "  " << hits10 / hits1 << endl;
    cout << "Eff. Fi11 min: " << hits11 << "  " << hits11 / hits1 << endl;
    cout << "Eff. Fi12 min: " << hits12 << "  " << hits12 / hits1 << endl;
    cout << "Eff. Fi13 min: " << hits13 << "  " << hits13 / hits1 << endl;

    cout << "Eff. Fi10 max: " << hits10bc << "  " << hits10bc / hits1 << endl;
    cout << "Eff. Fi11 max: " << hits11bc << "  " << hits11bc / hits1 << endl;
    cout << "Eff. Fi12 max: " << hits12bc << "  " << hits12bc / hits1 << endl;
    cout << "Eff. Fi13 max: " << hits13bc << "  " << hits13bc / hits1 << endl;

    fh_Tpat->Write();
    fh_Trigger->Write();
    fh_Cave_position->Write();
    fh_TOFDOR->Write();
    fh_SEE->Write();
    fh_IC->Write();

    fh_target_xy->Write();
    fh_target_px->Write();
    fh_target_py->Write();
    fh_target_pz->Write();
    fh_target_p->Write();

    fh_px_He->Write();
    fh_py_He->Write();
    fh_pz_He->Write();
    fh_p_He->Write();
    fh_px_C->Write();
    fh_py_C->Write();
    fh_pz_C->Write();
    fh_p_C->Write();

    fh_chi2->Write();

    fh_dx->Write();
    fh_dy->Write();
    fh_dz->Write();
    fh_dpx->Write();
    fh_dpy->Write();
    fh_dpz->Write();
    fh_dp->Write();
    fh_thetax_dpx->Write();
    fh_thetax_dpx_abs->Write();
    fh_thetay_dpy->Write();
    fh_thetax_dpy->Write();
    fh_thetay_dpx->Write();
    fh_x_dpx->Write();
    fh_y_dpy->Write();
    fh_dpy_dpx->Write();

    if (fMappedItems.at(DET_CALIFA))
    {
        fh_califa_energy->Write();
    }

    if (fHitItems.at(DET_TOFD))
    {
        fh_xy_tofd->Write();
        fh_xy_tofd_ac->Write();
        fh_tofd_charge->Write();
        fh_tofd_charge_ac->Write();
        fh_TimePreviousEvent->Write();
        fh_tofd_time->Write();
        fh_tofd_time_ac->Write();
        fh_tofd_mult->Write();
        fh_tofd_mult_ac->Write();
        fh_tofd_q2_vs_q1->Write();
        fh_tofd_q2_vs_q1_ac->Write();
    }

    for (Int_t ifibcount = 0; ifibcount < NOF_FIB_DET; ifibcount++)
    {
        if (fCalItems.at(ifibcount + DET_FI_FIRST) || fHitItems.at(ifibcount + DET_FI_FIRST))
        {
            fh_xy_Fib[ifibcount]->Write();
            fh_xy_Fib_ac[ifibcount]->Write();
            fh_mult_Fib[ifibcount]->Write();
            fh_mult_Fib_ac[ifibcount]->Write();
            fh_ToT_Fib[ifibcount]->Write();
            fh_ToT_Fib_ac[ifibcount]->Write();
            fh_Fib_vs_Events[ifibcount]->Write();
            fh_Fib_vs_Events_ac[ifibcount]->Write();
            fh_Fibs_vs_Tofd[ifibcount]->Write();
            fh_Fibs_vs_Tofd_ac[ifibcount]->Write();
            fh_Fib_ToF[ifibcount]->Write();
            fh_Fib_ToF_ac[ifibcount]->Write();
            fh_ToF_vs_Events[ifibcount]->Write();
            fh_ToF_vs_Events_ac[ifibcount]->Write();
            fh_Fib_Time[ifibcount]->Write();
            fh_Fib_Time_ac[ifibcount]->Write();
        }
    }

    fh_Fib13_vs_Fib11->Write();
    fh_Fib13_vs_Fib11_dx->Write();
    fh_Fib11_vs_Fib3a->Write();
    fh_Fib11_vs_Fib3a_dx->Write();
    fh_Fib10_vs_Fib12->Write();
    fh_Fib10_vs_Fib12_dx->Write();
    fh_Fib12_vs_Fib3b->Write();
    fh_Fib12_vs_Fib3b_dx->Write();

    fh_Fib13_vs_Fib11_back->Write();
    fh_Fib13_vs_Fib11_dx_back->Write();
    fh_Fib11_vs_Fib3a_back->Write();
    fh_Fib11_vs_Fib3a_dx_back->Write();
    fh_Fib10_vs_Fib12_back->Write();
    fh_Fib10_vs_Fib12_dx_back->Write();
    fh_Fib12_vs_Fib3b_back->Write();
    fh_Fib12_vs_Fib3b_dx_back->Write();

    fh_theta26->Write();
    fh_Erel->Write();
    fh_theta26_simu->Write();
    fh_Erel_simu->Write();
    fh_chiy_vs_chix->Write();
    fh_dErel_vs_x->Write();
    fh_dErel_vs_y->Write();

    for (Int_t i = 0; i < ndet; i++)
    {
        fh_xy[i]->Write();
        fh_p_vs_x[i]->Write();
        fh_p_vs_x_test[i]->Write();
    }
}

ClassImp(R3BTrackS454)
