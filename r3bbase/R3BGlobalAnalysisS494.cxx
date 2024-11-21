/******************************************************************************
 *   Copyright (C) 2019 GSI Helmholtzzentrum für Schwerionenforschung GmbH    *
 *   Copyright (C) 2019 Members of R3B Collaboration                          *
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
// -----                  R3BGlobalAnalysisS494                -----
// -----          Created April 13th 2016 by M.Heil       -----
// ------------------------------------------------------------

/*
 * This task should fill histograms with detector variables which allow
 * to test the detectors online
 *
 */

#include "R3BCalifaClusterData.h"
#include "R3BCalifaCrystalCalData.h"
#include "R3BCalifaMappedData.h"

#include "R3BGlobalAnalysisS494.h"

#include "R3BTofdHitData.h"

#include "R3BFiberMAPMTHitData.h"

#include "R3BEventHeader.h"
#include "R3BTCalEngine.h"

#include "R3BFiberMAPMTCalData.h"
#include "R3BFiberMAPMTHitData.h"
#include "R3BFiberMAPMTMappedData.h"

#include "R3BBunchedFiberCalData.h"
#include "R3BBunchedFiberHitData.h"
#include "R3BBunchedFiberMappedData.h"

#include "R3BMCTrack.h"
#include "R3BTrack.h"
#include "R3BTrackData.h"

#include "FairLogger.h"
#include "FairRootManager.h"
#include "FairRunAna.h"
#include "FairRunOnline.h"
#include "FairRuntimeDb.h"
#include "R3BLogger.h"

#include "TCanvas.h"
#include "TH1F.h"
#include "TH2F.h"

#include "TFile.h"

#include "tracker_routines.h"

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

R3BGlobalAnalysisS494::R3BGlobalAnalysisS494()
    : R3BGlobalAnalysisS494("GlobalAnalysis", 1)
{
}

R3BGlobalAnalysisS494::R3BGlobalAnalysisS494(const char* name, Int_t iVerbose)
    : FairTask(name, iVerbose)
    , fTrigger(-1)
    , fTpat(-1)
    , fCuts(0)
    , fGhost(0)
    , fPairs(1)
    , fSimu(0)
    , fGraphCuts(0)
    , fB(-1710)
    , fcut_chiX(5000)
    , fcut_chiY(5000)
    , fThetaGrazing(0.79)
    , fxfibcut(0.1)
    , fyfibcut(0.1)
    , fvis(0)
    , fEvsE("")
    , fCalifaCut("")
    , fNEvents(0)
    , fMappedItemsCalifa(NULL)
    , fCalItemsCalifa(NULL)
    , fHitItemsCalifa(NULL)
    , fCalifaMappedItems(new TClonesArray("R3BCalifaMappedData"))
    , fCalifaCalItems(new TClonesArray("R3BCalifaCrystalCalData"))
    , fCalifaHitItems(new TClonesArray("R3BCalifaClusterData"))
    , fNofCalifaMappedItems(0)
    , fNofCalifaCalItems(0)
    , fNofCalifaHitItems(0)
{
}

R3BGlobalAnalysisS494::~R3BGlobalAnalysisS494() {}

InitStatus R3BGlobalAnalysisS494::Init()
{
    // Initialize random number:
    std::srand(std::time(0)); // use current time as seed for random generator

    cout << "R3BGlobalAnalysisS494::Init ";

    // try to get a handle on the EventHeader. EventHeader may not be
    // present though and hence may be null. Take care when using.

    FairRootManager* mgr = FairRootManager::Instance();
    if (NULL == mgr)
        LOG(fatal) << "FairRootManager not found";

    header = dynamic_cast<R3BEventHeader*>(mgr->GetObject("EventHeader."));
    if (header)
        R3BLOG(info, "EventHeader. was found");
    else
        R3BLOG(info, "EventHeader. was not found");

    FairRunAna* run = FairRunAna::Instance();

    // Get objects for detectors on all levels
    fMCTrack = (TClonesArray*)mgr->GetObject("MCTrack");
    fTrack = (TClonesArray*)mgr->GetObject("TrackData");
    maxevent = mgr->CheckMaxEventNo();

    // get access to data of Califa
    fMappedItemsCalifa = (TClonesArray*)mgr->GetObject("CalifaMappedData");
    if (!fMappedItemsCalifa)
        LOG(warning) << "GlobalAnalysis: CalifaMappedData not found";
    fCalItemsCalifa = (TClonesArray*)mgr->GetObject("CalifaCrystalCalData");
    if (!fCalItemsCalifa)
        LOG(warning) << "GlobalAnalysis: CalifaCrystalCalData not found";
    fHitItemsCalifa = (TClonesArray*)mgr->GetObject("CalifaClusterData");
    if (!fHitItemsCalifa)
        LOG(warning) << "GlobalAnalysis: CalifaClusterData not found";

    if (fMappedItemsCalifa)
        mgr->Register("CalifaMappedData", "Land", fCalifaMappedItems, kTRUE);
    if (fCalItemsCalifa)
        mgr->Register("CalifaCrystalCalData", "Land", fCalifaCalItems, kTRUE);
    if (fHitItemsCalifa)
        mgr->Register("CalifaClusterData", "Land", fCalifaHitItems, kTRUE);

    //------------------------------------------------------------------------
    // graphical cuts
    //------------------------------------------------------------------------
    if (fGraphCuts)
    {
        LOG(info) << "READING GRAPHICAL CUTS " << endl;
        cut_EHe_EC = NULL;
        if (!fEvsE)
        {
            LOG(warn) << "No file with graphical cuts  fEvsE loaded!" << endl;
        }
        else
        {
            TFile* fevse = TFile::Open(fEvsE, "read");
            cut_EHe_EC = dynamic_cast<TCutG*>(fevse->Get("evse_new"));

            LOG(info) << "GRAPHICAL CUTS fEvsE ARE READ" << endl;
        }
    }
    if (!fSimu)
    {
        cut_CalifaTof = NULL;
        if (!fCalifaCut)
        {
            LOG(warn) << "No file with graphical cuts fCalifaCut loaded!" << endl;
        }
        else
        {
            TFile* fcalifacut = TFile::Open(fCalifaCut, "read");
            cut_CalifaTof = dynamic_cast<TCutG*>(fcalifacut->Get("cluster_sorting"));

            LOG(info) << "GRAPHICAL CUTS fCalifaCut ARE READ" << endl;
        }
    }
    //------------------------------------------------------------------------
    // create histograms of all detectors
    //------------------------------------------------------------------------

    //    char strNameC[255];
    //    sprintf(strNameC, "Cave_C_position");
    //    fh_Cave_position = new TH2F(strNameC, "", 2100, -100., 2000., 1000, -500., 500.);

    //-----------------------------------------------------------------------
    //-----------------------------------------------------------------------
    //

    R3BLOG(debug, "Starting to define histograms");
    fh_px_He = new TH1F("px_He", " px He", 500, -500., 500);
    fh_px_He->GetXaxis()->SetTitle("px / MeV/c");
    fh_px_He->GetYaxis()->SetTitle("counts");

    fh_py_He = new TH1F("py_He", " py He", 500, -500., 500);
    fh_py_He->GetXaxis()->SetTitle("py / MeV/c");
    fh_py_He->GetYaxis()->SetTitle("counts");

    fh_pz_He = new TH1F("pz_He", " pz He", 1000, 2000., 5000.);
    fh_pz_He->GetXaxis()->SetTitle("pz / MeV/c");
    fh_pz_He->GetYaxis()->SetTitle("counts");

    fh_p_He = new TH1F("p_He", " p He", 1000, 2000., 5000.);
    fh_p_He->GetXaxis()->SetTitle("pz / MeV/c");
    fh_p_He->GetYaxis()->SetTitle("counts");

    fh_px_C = new TH1F("px_C", " px C", 500, -500., 500);
    fh_px_C->GetXaxis()->SetTitle("px / MeV/c");
    fh_px_C->GetYaxis()->SetTitle("counts");

    fh_py_C = new TH1F("py_C", " py C", 500, -500., 500);
    fh_py_C->GetXaxis()->SetTitle("py / MeV/c");
    fh_py_C->GetYaxis()->SetTitle("counts");

    fh_pz_C = new TH1F("pz_C", " pz C", 2000, 10000., 16000.);
    fh_pz_C->GetXaxis()->SetTitle("pz / MeV/c");
    fh_pz_C->GetYaxis()->SetTitle("counts");

    fh_p_C = new TH1F("p_C", " p C", 200, 10000., 16000.);
    fh_p_C->GetXaxis()->SetTitle("pz / MeV/c");
    fh_p_C->GetYaxis()->SetTitle("counts");

    fh_px_O = new TH1F("px_O", " px O", 500, -500., 500);
    fh_px_O->GetXaxis()->SetTitle("px / MeV/c");
    fh_px_O->GetYaxis()->SetTitle("counts");

    fh_py_O = new TH1F("py_O", " py O", 500, -500., 500);
    fh_py_O->GetXaxis()->SetTitle("py / MeV/c");
    fh_py_O->GetYaxis()->SetTitle("counts");

    fh_pz_O = new TH1F("pz_O", " pz O", 1000, 17000., 18000.);
    fh_pz_O->GetXaxis()->SetTitle("pz / MeV/c");
    fh_pz_O->GetYaxis()->SetTitle("counts");

    fh_p_O = new TH1F("p_O", " p O", 1000, 17000., 18000.);
    fh_p_O->GetXaxis()->SetTitle("pz / MeV/c");
    fh_p_O->GetYaxis()->SetTitle("counts");

    fh_px_He_nc = new TH1F("px_He_nc", " px He all chi2", 500, -500., 500);
    fh_px_He_nc->GetXaxis()->SetTitle("px / MeV/c");
    fh_px_He_nc->GetYaxis()->SetTitle("counts");

    fh_py_He_nc = new TH1F("py_He_nc", " py He all chi2", 500, -500., 500);
    fh_py_He_nc->GetXaxis()->SetTitle("py / MeV/c");
    fh_py_He_nc->GetYaxis()->SetTitle("counts");

    fh_pz_He_nc = new TH1F("pz_He_nc", " pz He all chi2", 1000, 2000., 5000.);
    fh_pz_He_nc->GetXaxis()->SetTitle("pz / MeV/c");
    fh_pz_He_nc->GetYaxis()->SetTitle("counts");

    fh_p_He_nc = new TH1F("p_He_nc", " p He all chi2", 1000, 2000., 5000.);
    fh_p_He_nc->GetXaxis()->SetTitle("pz / MeV/c");
    fh_p_He_nc->GetYaxis()->SetTitle("counts");

    fh_px_C_nc = new TH1F("px_C_nc", " px C all chi2", 500, -500., 500);
    fh_px_C_nc->GetXaxis()->SetTitle("px / MeV/c");
    fh_px_C_nc->GetYaxis()->SetTitle("counts");

    fh_py_C_nc = new TH1F("py_C_nc", " py C all chi2", 500, -500., 500);
    fh_py_C_nc->GetXaxis()->SetTitle("py / MeV/c");
    fh_py_C_nc->GetYaxis()->SetTitle("counts");

    fh_pz_C_nc = new TH1F("pz_C_nc", " pz C all chi2", 2000, 10000., 16000.);
    fh_pz_C_nc->GetXaxis()->SetTitle("pz / MeV/c");
    fh_pz_C_nc->GetYaxis()->SetTitle("counts");

    fh_p_C_nc = new TH1F("p_C_nc", " p C all chi2", 2000, 10000., 16000.);
    fh_p_C_nc->GetXaxis()->SetTitle("pz / MeV/c");
    fh_p_C_nc->GetYaxis()->SetTitle("counts");

    fh_px_O_nc = new TH1F("px_O_nc", " px O all chi2", 500, -500., 500);
    fh_px_O_nc->GetXaxis()->SetTitle("px / MeV/c");
    fh_px_O_nc->GetYaxis()->SetTitle("counts");

    fh_py_O_nc = new TH1F("py_O_nc", " py O all chi2", 500, -500., 500);
    fh_py_O_nc->GetXaxis()->SetTitle("py / MeV/c");
    fh_py_O_nc->GetYaxis()->SetTitle("counts");

    fh_pz_O_nc = new TH1F("pz_O_nc", " pz O all chi2", 1000, 16000., 19000.);
    fh_pz_O_nc->GetXaxis()->SetTitle("pz / MeV/c");
    fh_pz_O_nc->GetYaxis()->SetTitle("counts");

    fh_p_O_nc = new TH1F("p_O_nc", " p O all chi2", 1000, 16000., 19000.);
    fh_p_O_nc->GetXaxis()->SetTitle("pz / MeV/c");
    fh_p_O_nc->GetYaxis()->SetTitle("counts");

    fh_target_xy_nc = new TH2F("target_xy_nc", "target xy all chi2", 200, -1., 1., 200, -1., 1.);
    fh_target_xy_nc->GetXaxis()->SetTitle("x / cm");
    fh_target_xy_nc->GetYaxis()->SetTitle("y / cm");

    fh_target_xx = new TH2F("targetX_fi23X", "target x vs fi23 x ", 300, -1.5, 1.5, 429, -6.006, 6.006);
    fh_target_xx->GetXaxis()->SetTitle("fi23X / cm");
    fh_target_xx->GetYaxis()->SetTitle("targetX / cm");

    fh_target_yy = new TH2F("targetY_fi23Y", "target y vs fi23 y ", 300, -1.5, 1.5, 429, -6.006, 6.006);
    fh_target_yy->GetXaxis()->SetTitle("fi23Y / cm");
    fh_target_yy->GetYaxis()->SetTitle("targetY / cm");

    fh_target_xy = new TH2F("target_xy", "target xy ", 200, -1., 1., 200, -1., 1.);
    fh_target_xy->GetXaxis()->SetTitle("x / cm");
    fh_target_xy->GetYaxis()->SetTitle("y / cm");

    fh_py_yfi23 = new TH2F("fh_py_yfi23", "tracker fh_py_yfi23 ", 429, -6.006, 6.006, 2000, -1000, 1000);
    fh_py_yfi23->GetYaxis()->SetTitle("pY / MeV/c");
    fh_py_yfi23->GetXaxis()->SetTitle("yfib23 / cm");

    fh_px_xfi23 = new TH2F("fh_px_xfi23", "tracker fh_px_xfi23 ", 429, -6.006, 6.006, 2000, -1000, 1000);
    fh_px_xfi23->GetYaxis()->SetTitle("px / MeV/c");
    fh_px_xfi23->GetXaxis()->SetTitle("xfib23 / cm");

    fh_p_vs_x0 = new TH2F("fh_p_x0", "p vs x0 ", 60, -1.5, 1.5, 800, 2000, 18000);
    fh_p_vs_x0->GetYaxis()->SetTitle("p / MeV/c");
    fh_p_vs_x0->GetXaxis()->SetTitle("xtarget / cm");

    fh_psum_vs_x0 = new TH2F("fh_psum_x0", "psum vs x0 ", 60, -1.5, 1.5, 200, 16000, 18000);
    fh_psum_vs_x0->GetYaxis()->SetTitle("psum / MeV/c");
    fh_psum_vs_x0->GetXaxis()->SetTitle("xtarget / cm");

    fh_p_vs_y0 = new TH2F("fh_p_y0", "p vs y0 ", 60, -1.5, 1.5, 800, 2000, 18000);
    fh_p_vs_y0->GetYaxis()->SetTitle("p / MeV/c");
    fh_p_vs_y0->GetXaxis()->SetTitle("ytarget / cm");

    fh_psum_vs_y0 = new TH2F("fh_psum_y0", "psum vs y0 ", 60, -1.5, 1.5, 200, 16000, 18000);
    fh_psum_vs_y0->GetYaxis()->SetTitle("psum / MeV/c");
    fh_psum_vs_y0->GetXaxis()->SetTitle("ytarget / cm");

    if (fSimu)
    {
        fh_target_xy_MC = new TH2F("target_xy_MC", "target xy MC", 200, -1., 1., 200, -1., 1.);
        fh_target_xy_MC->GetXaxis()->SetTitle("x / cm");
        fh_target_xy_MC->GetYaxis()->SetTitle("y / cm");

        fh_px_px_mc = new TH2F("px_px_MC", "px He vs. px C MC", 400, -500, 500, 400, -500., 500.);
        fh_px_px_mc->GetXaxis()->SetTitle("px / MeV/c");
        fh_px_px_mc->GetYaxis()->SetTitle("px / MeV/c");

        fh_py_py_mc = new TH2F("py_py_MC", "py He vs. py C MC", 400, -500, 500, 400, -500., 500.);
        fh_py_py_mc->GetXaxis()->SetTitle("py / MeV/c");
        fh_py_py_mc->GetYaxis()->SetTitle("py / MeV/c");

        fh_pz_pz_mc = new TH2F("pz_pz_MC", "pz He vs. pz C MC", 400, 2000, 6000, 1000, 10000., 20000.);
        fh_pz_pz_mc->GetXaxis()->SetTitle("pz / MeV/c");
        fh_pz_pz_mc->GetYaxis()->SetTitle("pz / MeV/c");

        fh_dx = new TH1F("tracker_dx", "tracker dx ", 200, -5., 5);
        fh_dx->GetXaxis()->SetTitle("dx / cm");
        fh_dx->GetYaxis()->SetTitle("counts");

        fh_dy = new TH1F("tracker_dy", "tracker dy ", 200, -5., 5.);
        fh_dy->GetXaxis()->SetTitle("dy / cm");
        fh_dy->GetYaxis()->SetTitle("counts");

        fh_dz = new TH1F("tracker_dz", "tracker dz ", 200, -5., 5.);
        fh_dz->GetXaxis()->SetTitle("dz / cm");
        fh_dz->GetYaxis()->SetTitle("counts");

        fh_dpxO = new TH1F("tracker_dpxO", "tracker dpx O", 10000, -1000, 1000);
        fh_dpxO->GetXaxis()->SetTitle("dpx / percent");
        fh_dpxO->GetYaxis()->SetTitle("counts");

        fh_dpyO = new TH1F("tracker_dpyO", "tracker dpy O", 10000, -1000, 1000);
        fh_dpyO->GetXaxis()->SetTitle("dpy / percent");
        fh_dpyO->GetYaxis()->SetTitle("counts");

        fh_dpzO = new TH1F("tracker_dpzO", "tracker dpz O", 400, -20, 20);
        fh_dpzO->GetXaxis()->SetTitle("dpz / percent");
        fh_dpzO->GetYaxis()->SetTitle("counts");

        fh_dpO = new TH1F("tracker_dpO", "tracker dp O", 10000, -10, 10);
        fh_dpO->GetXaxis()->SetTitle("dp / percent");
        fh_dpO->GetYaxis()->SetTitle("counts");

        fh_dpxC = new TH1F("tracker_dpxC", "tracker dpx C", 2000, -100, 100);
        fh_dpxC->GetXaxis()->SetTitle("dpx / percent");
        fh_dpxC->GetYaxis()->SetTitle("counts");

        fh_dpyC = new TH1F("tracker_dpyC", "tracker dpy C", 2000, -100, 100);
        fh_dpyC->GetXaxis()->SetTitle("dpy / percent");
        fh_dpyC->GetYaxis()->SetTitle("counts");

        fh_dpzC = new TH1F("tracker_dpzC", "tracker dpzC", 400, -20, 20);
        fh_dpzC->GetXaxis()->SetTitle("dpz / percent");
        fh_dpzC->GetYaxis()->SetTitle("counts");

        fh_dpC = new TH1F("tracker_dpC", "tracker dpC", 400, -20, 20);
        fh_dpC->GetXaxis()->SetTitle("dp / percent");
        fh_dpC->GetYaxis()->SetTitle("counts");

        fh_dp = new TH1F("tracker_dp", "tracker dp", 4000, -10, 10);
        fh_dp->GetXaxis()->SetTitle("dp / percent");
        fh_dp->GetYaxis()->SetTitle("counts");

        fh_dpxHe = new TH1F("tracker_dpxHe", "tracker Hedpx ", 2000, -100, 100);
        fh_dpxHe->GetXaxis()->SetTitle("dpx / percent");
        fh_dpxHe->GetYaxis()->SetTitle("counts");

        fh_dpyHe = new TH1F("tracker_dpyHe", "tracker Hedpy ", 2000, -100, 100);
        fh_dpyHe->GetXaxis()->SetTitle("dpy / percent");
        fh_dpyHe->GetYaxis()->SetTitle("counts");

        fh_dpzHe = new TH1F("tracker_dpzHe", "tracker dpz He", 400, -20, 20);
        fh_dpzHe->GetXaxis()->SetTitle("dpz / percent");
        fh_dpzHe->GetYaxis()->SetTitle("counts");

        fh_dpHe = new TH1F("tracker_dpHe", "tracker dp He", 400, -20, 20);
        fh_dpHe->GetXaxis()->SetTitle("dp / percent");
        fh_dpHe->GetYaxis()->SetTitle("counts");

        fh_dpy_dpx = new TH2F("tracker_dpy_dpx", "tracker dpy vs dpx ", 200, -100, 100, 200, -100., 100.);
        fh_dpy_dpx->GetXaxis()->SetTitle("dpx /  percent");
        fh_dpy_dpx->GetYaxis()->SetTitle("dpy / percent");

        // fh_dpy_dpx = new TH2F("tracker_dpy_dpx", "tracker py vs pz ", 2000, 16000., 18000, 800,-1000,1000);
        // fh_dpy_dpx->GetXaxis()->SetTitle("Pz / MeV/c");
        // fh_dpy_dpx->GetYaxis()->SetTitle("Py / MeV/c");

        fh_dp_2d = new TH2F("fh_dp_2d", "tracker dp vs y vs x ", 429, -6.006, 6.006, 429, -6.006, 6.006);
        fh_dp_2d->GetYaxis()->SetTitle("xfib23 / cm");
        fh_dp_2d->GetXaxis()->SetTitle("yfib23 / cm");

        fh_dpx_2d = new TH2F("fh_dpx_2d", "tracker dpx vs y vs x ", 429, -6.006, 6.006, 429, -6.006, 6.006);
        fh_dpx_2d->GetYaxis()->SetTitle("xfib23 / cm");
        fh_dpx_2d->GetXaxis()->SetTitle("yfib23 / cm");

        fh_dpy_2d = new TH2F("fh_dpy_2d", "tracker dpy vs y vs x ", 429, -6.006, 6.006, 429, -6.006, 6.006);
        fh_dpy_2d->GetYaxis()->SetTitle("xfib23 / cm");
        fh_dpy_2d->GetXaxis()->SetTitle("yfib23 / cm");

        fh_dpz_2d = new TH2F("fh_dpz_2d", "tracker dpz vs y vs x ", 429, -6.006, 6.006, 429, -6.006, 6.006);
        fh_dpz_2d->GetYaxis()->SetTitle("xfib23 / cm");
        fh_dpz_2d->GetXaxis()->SetTitle("yfib23 / cm");

        fh_thetax_dpx_C =
            new TH2F("tracker_thethax_dpx_C", "tracker theta_x vs dpx C", 2000, -200, 200, 200, -50., 50.);
        fh_thetax_dpx_C->GetXaxis()->SetTitle("dpx / MeV/c");
        fh_thetax_dpx_C->GetYaxis()->SetTitle("theta_x / mrad ");

        fh_thetay_dpy_C =
            new TH2F("tracker_thethay_dpy_C", "tracker theta_y vs dpy_C ", 2000, -200, 200, 200, -50., 50.);
        fh_thetay_dpy_C->GetXaxis()->SetTitle("dpy / MeV/c");
        fh_thetay_dpy_C->GetYaxis()->SetTitle("theta_y / mrad ");

        fh_thetax_dpx_He =
            new TH2F("tracker_thethax_dpx_He", "tracker theta_x vs dpx He", 2000, -200, 200, 200, -50., 50.);
        fh_thetax_dpx_He->GetXaxis()->SetTitle("dpx / MeV/c");
        fh_thetax_dpx_He->GetYaxis()->SetTitle("theta_x / mrad ");

        fh_thetay_dpy_He =
            new TH2F("tracker_thethay_dpy_He", "tracker theta_y vs dpy_He ", 2000, -200, 200, 200, -50., 50.);
        fh_thetay_dpy_He->GetXaxis()->SetTitle("dpy / MeV/c");
        fh_thetay_dpy_He->GetYaxis()->SetTitle("theta_y / mrad ");

        fh_thetax_dpy = new TH2F("tracker_thethax_dpy", "tracker theta_x vs dpy ", 2000, -200, 200, 200, -50., 50.);
        fh_thetax_dpy->GetXaxis()->SetTitle("dpy / MeV/c");
        fh_thetax_dpy->GetYaxis()->SetTitle("theta_x / mrad ");

        fh_thetay_dpx = new TH2F("tracker_thethay_dpx", "tracker theta_y vs dpx ", 2000, -200, 200, 200, -50., 50.);
        fh_thetay_dpx->GetXaxis()->SetTitle("dpx / MeV/c");
        fh_thetay_dpx->GetYaxis()->SetTitle("theta_y / mrad ");

        fh_x_dpx = new TH2F("tracker_px_dpx", "tracker px vs dpx ", 1000, -50, 50, 400, -500., 500.);
        fh_x_dpx->GetXaxis()->SetTitle("dpx / percent");
        fh_x_dpx->GetYaxis()->SetTitle("px / MeV/c");

        // fh_x_dpx = new TH2F("tracker_xfib23_dpx", "tracker xfi23 vs dpx ",  429, -6.006, 6.006, 2000, 0, 2);
        // fh_x_dpx->GetXaxis()->SetTitle("dpx / percent");
        // fh_x_dpx->GetYaxis()->SetTitle("xfib23 / cm");

        fh_y_dpy = new TH2F("tracker_py_dpy", "tracker py vs dpy ", 1000, -50, 50, 400, -500., 500.);
        fh_y_dpy->GetXaxis()->SetTitle("dpy / percent");
        fh_y_dpy->GetYaxis()->SetTitle("py / MeV/c");

        // fh_y_dpy = new TH2F("tracker_yfib23_dpy", "tracker yfib23 vs Py_MC/Py ", 429, -6.006, 6.006, 2000, 0, 2);
        // fh_y_dpy->GetYaxis()->SetTitle("Py_MC/Py");
        // fh_y_dpy->GetXaxis()->SetTitle("yfib23 / cm");

        fh_py_yfi23_MC = new TH2F("fh_py_yfi23_MC", "tracker fh_py_yfi23_MC ", 429, -6.006, 6.006, 2000, -1000, 1000);
        fh_py_yfi23_MC->GetYaxis()->SetTitle("pY / MeV/c");
        fh_py_yfi23_MC->GetXaxis()->SetTitle("yfib23 / cm");

        fh_px_xfi23_MC = new TH2F("fh_px_xfi23_MC", "tracker fh_px_xfi23_MC ", 429, -6.006, 6.006, 2000, -1000, 1000);
        fh_px_xfi23_MC->GetYaxis()->SetTitle("px / MeV/c");
        fh_px_xfi23_MC->GetXaxis()->SetTitle("xfib23 / cm");

        fh_py_pz_MC = new TH2F("py_pz_MC", "py vs. ptot MC", 2000, 16000., 18000., 800, -1000., 1000.);
        fh_py_pz_MC->GetXaxis()->SetTitle("ptot / MeV/c");
        fh_py_pz_MC->GetYaxis()->SetTitle("py / MeV/c");

        fh_phiMC_bc_cm = new TH1F("phiMC_bc_cm", "phiMC bc in cm-system", 360, 0., 360.);
        fh_phiMC_bc_cm->GetXaxis()->SetTitle("angle / degree");
        fh_phiMC_bc_cm->GetYaxis()->SetTitle("counts");

        fh_phiMC_bc_cm_polar = new TH2F("phiMC_bc_cm_polar", "phiMC_bc_cm_polar", 360, 0., 360., 100, 0., 10.);
        fh_phiMC_bc_cm_polar->GetXaxis()->SetTitle("angle / degree");
        fh_phiMC_bc_cm_polar->GetYaxis()->SetTitle("counts");

        fh_theta26_simu = new TH1F("theta26_simu", "theta 26 simulation", 500, 0., 5);
        fh_theta26_simu->GetXaxis()->SetTitle("angle / degree");
        fh_theta26_simu->GetYaxis()->SetTitle("counts");

        fh_minv_simu = new TH1F("minv_MC", "minv_MC ", 500, 14900, 14920);
        fh_minv_simu->GetXaxis()->SetTitle("minv / MeV");
        fh_minv_simu->GetYaxis()->SetTitle("counts");

        fh_psum_MC = new TH1F("sum_p_MC", "sum_p_MC", 2000, 16000., 18000.);
        fh_psum_MC->GetXaxis()->SetTitle("p / MeV/c");
        fh_psum_MC->GetYaxis()->SetTitle("count");

        fh_pzsum_MC = new TH1F("sum_pz_MC", "sum_pz_MC", 2000, 16000., 18000.);
        fh_pzsum_MC->GetXaxis()->SetTitle("pz / MeV/c");
        fh_pzsum_MC->GetYaxis()->SetTitle("count");

        fh_pxsum_MC = new TH1F("sum_px_MC", "sum_px_MC", 500, -500., 500.);
        fh_pxsum_MC->GetXaxis()->SetTitle("px / MeV/c");
        fh_pxsum_MC->GetYaxis()->SetTitle("count");

        fh_pysum_MC = new TH1F("sum_py_MC", "sum_py_MC", 500, -500., 500.);
        fh_pysum_MC->GetXaxis()->SetTitle("py / MeV/c");
        fh_pysum_MC->GetYaxis()->SetTitle("count");

        fh_dtheta = new TH1F("dtheta26_cm", "delta theta ", 1000, -5., 5.);
        fh_dtheta->GetYaxis()->SetTitle("counts");
        fh_dtheta->GetXaxis()->SetTitle("Delta theta / deg");

        fh_Erel_vs_thetaMC = new TH2F("Erel_vs_thetaMC", "Erel vs. theta MC", 125, 0., 5, 200, 0, 20.);
        fh_Erel_vs_thetaMC->GetXaxis()->SetTitle("angle / deg");
        fh_Erel_vs_thetaMC->GetYaxis()->SetTitle("Erel / MeV");

        fh_Erel_vs_phibcMC = new TH2F("Erel_vs_phibccmMC", "Erel vs. phi_bc_cm MCtrack", 200, 0, 20., 400, 0., 400.);
        fh_Erel_vs_phibcMC->GetYaxis()->SetTitle("phi angle bc / deg");
        fh_Erel_vs_phibcMC->GetXaxis()->SetTitle("Erel / MeV");

        fh_Erel_simu = new TH1F("Erel_simu", "Erel simulation", 2000, 0., 20);
        fh_Erel_simu->GetXaxis()->SetTitle("Erel / MeV");
        fh_Erel_simu->GetYaxis()->SetTitle("counts");

        fh_ErelB_vs_phibc_MC = new TH2F("ErelB_vs_phibccm_MC", "ErelB vs. phi_bc_cm MC", 2000, 0, 20., 400, 0., 400.);
        fh_ErelB_vs_phibc_MC->GetYaxis()->SetTitle("phi angle bc / deg");
        fh_ErelB_vs_phibc_MC->GetXaxis()->SetTitle("Erel / MeV");

        fh_erelB_vs_ptransHe_MC =
            new TH2F("erelB_vs_ptransHe_MC", "ErelB vs. ptransHe MC", 2000, -1000, 1000, 150, 0., 15);
        fh_erelB_vs_ptransHe_MC->GetXaxis()->SetTitle("ptrans / MeV/c");
        fh_erelB_vs_ptransHe_MC->GetYaxis()->SetTitle("Erel / MeV");

        fh_erelB_vs_ptransC_MC =
            new TH2F("erelB_vs_ptransC_MC", "ErelB vs. ptransC MC", 2000, -1000, 1000, 150, 0., 15);
        fh_erelB_vs_ptransC_MC->GetXaxis()->SetTitle("ptrans / MeV/c");
        fh_erelB_vs_ptransC_MC->GetYaxis()->SetTitle("Erel / MeV");

        fh_thetaB_4He_cm_MC = new TH2F("thetaB_4He_cm_MC", "ThetaB 4He in cm-system MC", 200, -1., 1., 2000, 0, 20);
        fh_thetaB_4He_cm_MC->GetXaxis()->SetTitle("cos(theta_cm)");
        fh_thetaB_4He_cm_MC->GetYaxis()->SetTitle("ErelB / MeV");

        fh_phiB_4He_cm_MC = new TH2F("phiB_4He_cm_MC", "phiB 4He in cm-system MC", 360, 0., 360., 2000, 0, 20);
        fh_phiB_4He_cm_MC->GetXaxis()->SetTitle("angle / degree");
        fh_phiB_4He_cm_MC->GetYaxis()->SetTitle("Erel / MeV");

        fh_thetaB_12C_cm_MC = new TH2F("thetaB_12C_cm_MC", "ThetaB 12C in cm-system MC", 200, -1, 1., 2000, 0, 20);
        fh_thetaB_12C_cm_MC->GetXaxis()->SetTitle("cos(theta_cm)");
        fh_thetaB_12C_cm_MC->GetYaxis()->SetTitle("Erel / MeV");

        fh_phiB_12C_cm_MC = new TH2F("phiB_12C_cm_MC", "phiB 12C in cm-system MC", 360, 0., 360., 2000, 0, 20);
        fh_phiB_12C_cm_MC->GetXaxis()->SetTitle("angle / degree");
        fh_phiB_12C_cm_MC->GetYaxis()->SetTitle("Erel / MeV");

        fh_ErelB_vs_theta16O_MC = new TH2F("ErelB_vs_theta16O_MC", "ErelB vs. theta 16O* MC", 125, 0., 5, 2000, 0, 20.);
        fh_ErelB_vs_theta16O_MC->GetXaxis()->SetTitle("angle / deg");
        fh_ErelB_vs_theta16O_MC->GetYaxis()->SetTitle("ErelB / MeV");

        fh_psum_vs_dx0 = new TH2F("fh_psum_vs_dx0", "psum vs x0track-x0MC", 600, -3., 3., 1500, 16000., 19000);

        fh_pC_vs_dx0 = new TH2F("fh_pC_vs_dx0", "pzC vs x0track-x0MC", 600, -3., 3., 200, 12000., 14000);

        fh_x0_vs_dx0 = new TH2F("fh_x0_vs_x0", "x0track vs x0MC", 600, -1.5, 1.5, 600, -1.5, 1.5);
        fh_pC_vs_dx0_mc = new TH2F("fh_pC_vs_dx0_mc", "pzC_MC vs x0track-x0MC", 600, -3., 3., 200, 12000., 14000);

        fh_pC_vs_x0 = new TH2F("fh_pC_vs_x0", "pzC vs x0track", 600, -1.5, 1.5, 200, 12000., 14000);
        fh_pC_vs_x0_mc = new TH2F("fh_pC_vs_x0_mc", "pzC_MC vs x0MC", 600, -1.5, 1.5, 200, 12000., 14000);

        fh_dpC_vs_dx0 = new TH2F("fh_dpC_vs_dx0", "pzC-pzC_MC vs x0track-x0MC", 600, -3., 3., 2000, -1000, 1000);

        fh_pHe_vs_dx0 = new TH2F("fh_pHe_vs_dx0", "pzHe vs x0track-x0MC", 600, -3., 3., 200, 3000., 5000);
        fh_pHe_vs_dx0_mc = new TH2F("fh_pHe_vs_dx0_mc", "pzHe_MC vs x0track-x0MC", 600, -3., 3., 200, 3000., 5000);

        fh_pHe_vs_x0 = new TH2F("fh_pHe_vs_x0", "pzHe vs x0track", 600, -3., 3., 200, 3000., 5000);
        fh_pHe_vs_x0_mc = new TH2F("fh_pHe_vs_x0_mc", "pzHe_MC vs x0MC", 600, -1.5, 1.5, 200, 3000., 5000);

        fh_dpHe_vs_dx0 = new TH2F("fh_dpHe_vs_dx0", "pzHe-pzHe_MC vs x0track-x0MC", 600, -3., 3., 2000, -1000, 1000);
    }
    fh_dErel = new TH1F("dErel", "delta Erel ", 4000, -10, 10); // 400, -100, 100.);
    fh_dErel->GetYaxis()->SetTitle("counts");
    fh_dErel->GetXaxis()->SetTitle("Delta Erel / %");

    fh_thetax_px = new TH2F("thethax_px", "theta_x vs px ", 500, -500, 500, 200, -100., 100.);
    fh_thetax_px->GetXaxis()->SetTitle("px / MeV/c");
    fh_thetax_px->GetYaxis()->SetTitle("theta_x ");

    fh_thetay_py = new TH2F("thethay_py", "theta_y vs py ", 500, -500, 500, 200, -100., 100.);
    fh_thetay_py->GetXaxis()->SetTitle("py / MeV/c");
    fh_thetay_py->GetYaxis()->SetTitle("theta_y ");

    fh_thetax_py = new TH2F("thethax_py", "theta_x vs py ", 200, -100, 100, 200, -100., 100.);
    fh_thetax_py->GetXaxis()->SetTitle("py / MeV/c");
    fh_thetax_py->GetYaxis()->SetTitle("theta_x ");

    fh_thetay_px = new TH2F("thethay_px", "theta_y vs px ", 200, -100, 100, 200, -100., 100.);
    fh_thetay_px->GetXaxis()->SetTitle("px / MeV/c");
    fh_thetay_px->GetYaxis()->SetTitle("theta_y ");

    fh_px_px = new TH2F("px_px", "px He vs. px C", 600, -600, 600, 600, -600., 600.);
    fh_px_px->GetXaxis()->SetTitle("px / MeV/c");
    fh_px_px->GetYaxis()->SetTitle("px / MeV/c");

    fh_py_py = new TH2F("py_py", "py He vs. py C", 600, -600, 600, 600, -600., 600.);
    fh_py_py->GetXaxis()->SetTitle("py / MeV/c");
    fh_py_py->GetYaxis()->SetTitle("py / MeV/c");

    fh_pz_pz = new TH2F("pz_pz", "pz He vs. pz C", 400, 0, 10000, 400, 10000., 20000.);
    fh_pz_pz->GetXaxis()->SetTitle("pz / MeV/c");
    fh_pz_pz->GetYaxis()->SetTitle("pz / MeV/c");

    fh_p_p = new TH2F("p_p", "p He vs. p C", 400, 0, 10000, 400, 10000., 20000.);
    fh_p_p->GetXaxis()->SetTitle("p / MeV/c");
    fh_p_p->GetYaxis()->SetTitle("p / MeV/c");

    fh_py_pz = new TH2F("py_pz", "py vs. ptot ", 2000, 16000., 18000., 800, -1000., 1000.);
    fh_py_pz->GetXaxis()->SetTitle("ptot / MeV/c");
    fh_py_pz->GetYaxis()->SetTitle("py / MeV/c");

    fh_px_px_nc = new TH2F("px_px_nc", "px He vs. px C all chi2", 600, -600, 600, 600, -600., 600.);
    fh_px_px_nc->GetXaxis()->SetTitle("px / MeV/c");
    fh_px_px_nc->GetYaxis()->SetTitle("px / MeV/c");

    fh_py_py_nc = new TH2F("py_py_nc", "py He vs. py C all chi2", 600, -600, 600, 600, -600., 600.);
    fh_py_py_nc->GetXaxis()->SetTitle("py / MeV/c");
    fh_py_py_nc->GetYaxis()->SetTitle("py / MeV/c");

    fh_pz_pz_nc = new TH2F("pz_pz_nc", "pz He vs. pz C all chi2", 400, 2000, 6000, 500, 10000., 20000.);
    fh_pz_pz_nc->GetXaxis()->SetTitle("pz / MeV/c");
    fh_pz_pz_nc->GetYaxis()->SetTitle("pz / MeV/c");

    // Analysis
    fh_chi2 = new TH1F("chi2", "chi2 ", 1000, 0., 1000);
    fh_chi2->GetXaxis()->SetTitle("Chi2");
    fh_chi2->GetYaxis()->SetTitle("counts");

    fh_chiy_vs_chix_nc = new TH2F("chiC_vs_chiHe", "chi2 C vs chi2 He", 100, 0., 50000., 100, 0., 50000.);
    fh_chiy_vs_chix_nc->GetXaxis()->SetTitle("chi2 He");
    fh_chiy_vs_chix_nc->GetYaxis()->SetTitle("chi2 C");

    fh_chiy_vs_chix = new TH2F("chiC_vs_chiHe_after_cut", "chi2 C vs chi2 He after cut", 50, 0., 500., 50, 0., 500.);
    fh_chiy_vs_chix->GetXaxis()->SetTitle("chi2 He");
    fh_chiy_vs_chix->GetYaxis()->SetTitle("chi2 C");

    fh_theta26 = new TH1F("theta26", "theta between alpha and carbon ", 125, 0., 5);
    fh_theta26->GetXaxis()->SetTitle("angle / degree");
    fh_theta26->GetYaxis()->SetTitle("counts");

    fh_theta26_nc = new TH1F("theta26_nc", "theta between alpha and carbon all chi2 ", 125, 0., 5);
    fh_theta26_nc->GetXaxis()->SetTitle("angle / degree");
    fh_theta26_nc->GetYaxis()->SetTitle("counts");

    fh_theta26_vs_chi = new TH2F("theta_vs_chi", "theta vs chi", 500, 0., 500., 125, 0., 5);
    fh_theta26_vs_chi->GetXaxis()->SetTitle("chi ");
    fh_theta26_vs_chi->GetYaxis()->SetTitle("angle / degree");

    fh_theta_16O = new TH1F("theta_16O", "theta of excited 16O ", 500, 0., 5);
    fh_theta_16O->GetXaxis()->SetTitle("angle / degree");
    fh_theta_16O->GetYaxis()->SetTitle("counts");

    fh_phi_16O = new TH1F("phi_16O", "phi of excited 16O ", 360, 0., 360);
    fh_phi_16O->GetXaxis()->SetTitle("angle / degree");
    fh_phi_16O->GetYaxis()->SetTitle("counts");

    fh_theta_bc_cm = new TH2F("theta_bc_cm", "Theta bc in cm-system", 180, 0., 180., 2000, 0, 20);
    fh_theta_bc_cm->GetXaxis()->SetTitle("angle / degree");
    fh_theta_bc_cm->GetYaxis()->SetTitle("Erel / MeV");

    fh_theta_bc_cm_nc = new TH1F("theta_bc_cm_nc", "Theta bc in cm-system all chi2", 180, 0., 180.);
    fh_theta_bc_cm_nc->GetXaxis()->SetTitle("angle / degree");
    fh_theta_bc_cm_nc->GetYaxis()->SetTitle("counts");

    fh_phi_bc_cm_nc = new TH1F("phi_bc_cm_nc", "phi bc in cm-system all chi2", 360, 0., 360.);
    fh_phi_bc_cm_nc->GetXaxis()->SetTitle("angle / degree");
    fh_phi_bc_cm_nc->GetYaxis()->SetTitle("counts");

    fh_phi_bc_cm = new TH1F("phi_bc_cm", "phi bc in cm-system", 360, 0., 360.);
    fh_phi_bc_cm->GetXaxis()->SetTitle("angle / degree");
    fh_phi_bc_cm->GetYaxis()->SetTitle("counts");

    fh_phi_bc_cm_polar = new TH2F("phi_bc_cm_polar", "phi_bc_cm_polar", 360, 0., 360., 100, 0., 10.);
    fh_phi_bc_cm_polar->GetXaxis()->SetTitle("angle / degree");
    fh_phi_bc_cm_polar->GetYaxis()->SetTitle("counts");

    fh_minv = new TH1F("Target_exc", "target_exc ", 20000, -200, 200);
    fh_minv->GetXaxis()->SetTitle("target exc / MeV");
    fh_minv->GetYaxis()->SetTitle("counts");

    fh_Erel = new TH1F("Erel", "Erel ", 2000, 0., 20);
    fh_Erel->GetXaxis()->SetTitle("Erel / MeV");
    fh_Erel->GetYaxis()->SetTitle("counts");

    fh_Erel_vs_chi = new TH2F("Erel_vs_chi", "Erel vs chi", 1000, 0, 1000, 150, 0., 15);
    fh_Erel_vs_chi->GetYaxis()->SetTitle("Erel / MeV");
    fh_Erel_vs_chi->GetXaxis()->SetTitle("chi2");

    fh_phi26_vs_chi = new TH2F("phi16O_vs_phi_bc_cm", "phi_bc_cm vs phi_16O", 360, 0, 360, 360, 0., 360);
    fh_phi26_vs_chi->GetYaxis()->SetTitle("phi_bc_cm / deg");
    fh_phi26_vs_chi->GetXaxis()->SetTitle("phi16O / deg");

    fh_psum_vs_theta26 = new TH2F("psum_vs_theta26", "psum_vs_theta26", 125, 0., 5, 600, 16000., 19000);
    fh_psum_vs_theta26->GetYaxis()->SetTitle("theta26 / deg");
    fh_psum_vs_theta26->GetXaxis()->SetTitle("psum / MeV/c");

    fh_Erel_nc = new TH1F("Erel_nc", "Erel all chi2", 2000, 0., 20);
    fh_Erel_nc->GetXaxis()->SetTitle("Erel / MeV");
    fh_Erel_nc->GetYaxis()->SetTitle("counts");

    fh_ErelL = new TH1F("ErelL", "Erel links", 2000, 0., 20);
    fh_ErelL->GetXaxis()->SetTitle("Erel / MeV");
    fh_ErelL->GetYaxis()->SetTitle("counts");

    fh_ErelR = new TH1F("ErelR", "Erel rechts ", 2000, 0., 20);
    fh_ErelR->GetXaxis()->SetTitle("Erel / MeV");
    fh_ErelR->GetYaxis()->SetTitle("counts");

    fh_ErelB = new TH1F("ErelB", "Erel both sides ", 2000, 0., 20);
    fh_ErelB->GetXaxis()->SetTitle("Erel / MeV");
    fh_ErelB->GetYaxis()->SetTitle("counts");

    fh_ErelB_nc = new TH1F("ErelB_nc", "Erel both sides all chi2", 2000, 0., 20);
    fh_ErelB_nc->GetXaxis()->SetTitle("Erel / MeV");
    fh_ErelB_nc->GetYaxis()->SetTitle("counts");

    fh_Erel_vs_psum = new TH2F("Erel_vs_psum", "Erel vs. psum", 500, 16000., 18000., 200, 0, 20.);
    fh_Erel_vs_psum->GetXaxis()->SetTitle("psum / MeV/c");
    fh_Erel_vs_psum->GetYaxis()->SetTitle(" Erel / MeV");

    fh_Erel_vs_px = new TH2F("fh_Erel_vs_px", "Erel vs. px", 120, -600., 600., 200, 0, 20.);
    fh_Erel_vs_px->GetXaxis()->SetTitle("px / MeV/c");
    fh_Erel_vs_px->GetYaxis()->SetTitle(" Erel / MeV");

    fh_Erel_vs_py = new TH2F("fh_Erel_vs_py", "Erel vs. py", 120, -600., 600., 200, 0, 20.);
    fh_Erel_vs_py->GetXaxis()->SetTitle("py / MeV/c");
    fh_Erel_vs_py->GetYaxis()->SetTitle(" Erel / MeV");

    fh_Erel_vs_x0 = new TH2F("Erel_vs_x0", "Erel vs. x0", 300, -1.5, 1.5, 200, 0, 20.);
    fh_Erel_vs_x0->GetXaxis()->SetTitle("x0 / cm");
    fh_Erel_vs_x0->GetYaxis()->SetTitle(" Erel / MeV");

    fh_Erel_vs_y0 = new TH2F("Erel_vs_y0", "Erel vs. y0", 300, -1.5, 1.5, 200, 0, 20.);
    fh_Erel_vs_y0->GetXaxis()->SetTitle("y0 / cm");
    fh_Erel_vs_y0->GetYaxis()->SetTitle(" Erel / MeV");

    fh_Erel_vs_xfi23 = new TH2F("Erel_vs_xfi23", "Erel vs. xfi23C-xfi23He", 429, -6.006, 6.006, 200, 0, 20.);
    fh_Erel_vs_xfi23->GetXaxis()->SetTitle("xfi23C-xfi23He / cm");
    fh_Erel_vs_xfi23->GetYaxis()->SetTitle(" Erel / MeV");

    fh_Erel_vs_yfi23 = new TH2F("Erel_vs_yfi23", "Erel vs. yfi23C-yfi23He", 429, -6.006, 6.006, 200, 0, 20.);
    fh_Erel_vs_yfi23->GetXaxis()->SetTitle("yfi23C-yfi23He / cm");
    fh_Erel_vs_yfi23->GetYaxis()->SetTitle(" Erel / MeV");

    fh_Erel_vs_dr_tofd = new TH2F("Erel_vs_dr_tofd", "Erel vs. dr_tofd", 200, 0, 100, 200, 0, 20.);
    fh_Erel_vs_dr_tofd->GetXaxis()->SetTitle("sqrt((xC-xHe)^2+(yC-yHe)^2 / cm");
    fh_Erel_vs_dr_tofd->GetYaxis()->SetTitle(" Erel / MeV");

    fh_Erel_vs_dr_limit_tofd = new TH2F("Erel_vs_dr_tofd_limit", "Erel vs. dr limit_tofd", 200, 0, 100, 200, 0, 20.);
    fh_Erel_vs_dr_limit_tofd->GetXaxis()->SetTitle("sqrt((xC-xHe)^2+(yC-yHe)^2 / cm");
    fh_Erel_vs_dr_limit_tofd->GetYaxis()->SetTitle(" Erel / MeV");

    fh_Erel_limit_tofd = new TH1F("Erel_limit_tofd", "Erel limit_tofd", 2000, 0., 20);
    fh_Erel_limit_tofd->GetXaxis()->SetTitle("Erel / MeV");
    fh_Erel_limit_tofd->GetYaxis()->SetTitle("counts");

    fh_Erel_vs_dr_fi23 = new TH2F("Erel_vs_dr_fi23", "Erel vs. dr_fi23", 400, 0, 20, 200, 0, 20.);
    fh_Erel_vs_dr_fi23->GetXaxis()->SetTitle("sqrt((xC-xHe)^2+(yC-yHe)^2 / cm");
    fh_Erel_vs_dr_fi23->GetYaxis()->SetTitle(" Erel / MeV");

    fh_Erel_vs_dr_limit_fi23 = new TH2F("Erel_vs_dr_fi23_limit", "Erel vs. dr limit_fi23", 400, 0, 20, 200, 0, 20.);
    fh_Erel_vs_dr_limit_fi23->GetXaxis()->SetTitle("sqrt((xC-xHe)^2+(yC-yHe)^2 / cm");
    fh_Erel_vs_dr_limit_fi23->GetYaxis()->SetTitle(" Erel / MeV");

    fh_Erel_limit_fi23 = new TH1F("Erel_limit_fi23", "Erel limit_fi23", 2000, 0., 20);
    fh_Erel_limit_fi23->GetXaxis()->SetTitle("Erel / MeV");
    fh_Erel_limit_fi23->GetYaxis()->SetTitle("counts");

    fh_dr_vs_dr = new TH2F("dr_vs_dr", "drtofd vs drfi23", 400, 0, 10, 200, 0, 100.);
    fh_dr_vs_dr->GetXaxis()->SetTitle("sqrt((xC-xHe)^2+(yC-yHe)^2 fi23 / cm");
    fh_dr_vs_dr->GetYaxis()->SetTitle(" sqrt((xC-xHe)^2+(yC-yHe)^2 tofd / cm");

    fh_dr_vs_dr_limit = new TH2F("dr_vs_dr_limit", "drtofd vs drfi23_limit", 400, 0, 10, 200, 0, 100.);
    fh_dr_vs_dr_limit->GetXaxis()->SetTitle("sqrt((xC-xHe)^2+(yC-yHe)^2 fi23 / cm");
    fh_dr_vs_dr_limit->GetYaxis()->SetTitle(" sqrt((xC-xHe)^2+(yC-yHe)^2 tofd / cm");

    fh_mass_nc = new TH2F("mass_nc", "mass He vs mass C all chi2", 200, 10, 20, 500, 0, 5);
    fh_mass_nc->GetXaxis()->SetTitle("mass C");
    fh_mass_nc->GetYaxis()->SetTitle(" mass H");

    fh_mass = new TH2F("mass", "mass He vs mass C selected chi2", 200, 10, 20, 500, 0, 5);
    fh_mass->GetXaxis()->SetTitle("mass C");
    fh_mass->GetYaxis()->SetTitle(" mass H");

    fh_energy_nc = new TH2F("energy_nc", "energy He vs energy C all chi2", 200, 10, 20, 200, 0, 10);
    fh_energy_nc->GetXaxis()->SetTitle("energy C / GeV");
    fh_energy_nc->GetYaxis()->SetTitle(" energy H / GeV");

    fh_energy = new TH2F("energy", "energy He vs energy C selected chi2", 200, 10, 20, 200, 0, 10);
    fh_energy->GetXaxis()->SetTitle("energy C / GeV");
    fh_energy->GetYaxis()->SetTitle(" energy H / GeV");

    fh_psum = new TH1F("sum_p", "ptot", 2000, 16000., 18000.);
    fh_psum->GetXaxis()->SetTitle("p / MeV/c");
    fh_psum->GetYaxis()->SetTitle("count");

    fh_pzsum = new TH1F("sum_pz", "pztot", 2000, 16000., 18000.);
    fh_pzsum->GetXaxis()->SetTitle("pz / MeV/c");
    fh_pzsum->GetYaxis()->SetTitle("count");

    fh_psum_nc = new TH1F("sum_p_nc", "ptot all chi2", 2000, 16000., 18000.);
    fh_psum_nc->GetXaxis()->SetTitle("p / MeV/c");
    fh_psum_nc->GetYaxis()->SetTitle("count");

    fh_psum_vs_chi = new TH2F("sum_p_vs_chi", "ptot vs chi2", 1000, 0, 1000, 400, 16800., 17800.);
    fh_psum_vs_chi->GetYaxis()->SetTitle("p / MeV/c");
    fh_psum_vs_chi->GetXaxis()->SetTitle("chi2");

    fh_Erel_vs_event = new TH2F("Erel_vs_event", "Erel vs. event number", 10000, 0, 1e6, 250, 0., 10.);
    fh_Erel_vs_event->GetXaxis()->SetTitle("Event number");
    fh_Erel_vs_event->GetYaxis()->SetTitle("Erel / MeV");

    fh_psum_vs_event_nc =
        new TH2F("psum_vs_event_nc", "psum vs. event number all chi2", 100, 0, 1e6, 80, 17000., 17800.);
    fh_psum_vs_event_nc->GetXaxis()->SetTitle("Event number");
    fh_psum_vs_event_nc->GetYaxis()->SetTitle("p / MeV/c");

    fh_Erel_vs_theta16O = new TH2F("Erel_vs_theta16O", "Erel vs. theta 16O*", 125, 0., 5, 2000, 0, 20.);
    fh_Erel_vs_theta16O->GetXaxis()->SetTitle("angle / deg");
    fh_Erel_vs_theta16O->GetYaxis()->SetTitle("Erel / MeV");
    /*
        fh_Erel_vs_theta26 = new TH2F("Erel_vs_theta26", "Erel vs. opening angle", 125, 0., 5, 2000, 0, 20.);
        fh_Erel_vs_theta26->GetXaxis()->SetTitle("angle / deg");
        fh_Erel_vs_theta26->GetYaxis()->SetTitle("Erel / MeV");

        fh_Erel_vs_theta26_max = new TH2F("Erel_vs_theta26_max", "Erel vs. max opening angle", 125, 0., 5, 2000,
       0, 20.); fh_Erel_vs_theta26_max->GetXaxis()->SetTitle("angle / deg");
        fh_Erel_vs_theta26_max->GetYaxis()->SetTitle("Erel / MeV");
    */
    fh_ErelB_vs_theta16O = new TH2F("ErelB_vs_theta16O", "ErelB vs. theta 16O*", 125, 0., 5, 2000, 0, 20.);
    fh_ErelB_vs_theta16O->GetXaxis()->SetTitle("angle / deg");
    fh_ErelB_vs_theta16O->GetYaxis()->SetTitle("ErelB / MeV");

    fh_Erel_vs_theta16O_3He12C =
        new TH2F("Erel_vs_theta16O_3He12C", "Erel vs. theta 16O* for 3He+12C", 125, 0., 5, 500, 0, 20.);
    fh_Erel_vs_theta16O_3He12C->GetXaxis()->SetTitle("angle / deg");
    fh_Erel_vs_theta16O_3He12C->GetYaxis()->SetTitle("Erel / MeV");

    fh_pHe_vs_theta26 = new TH2F("pHe_vs_theta", "theta vs pHe", 400, 2000, 6000, 200, 0, 100.);
    fh_pHe_vs_theta26->GetYaxis()->SetTitle("angle / mrad");
    fh_pHe_vs_theta26->GetXaxis()->SetTitle("p / MeV/c");

    fh_ErelB_vs_phibc_bg =
        new TH2F("ErelB_vs_phibccm_bg", "ErelB vs. phi_bc_cm below grazing", 2000, 0, 20., 400, 0., 400.);
    fh_ErelB_vs_phibc_bg->GetYaxis()->SetTitle("phi angle bc / deg");
    fh_ErelB_vs_phibc_bg->GetXaxis()->SetTitle("Erel / MeV");

    fh_ErelB_vs_theta26_bg = new TH2F("ErelB_vs_theta_bg", "ErelB vs. theta below grazing", 125, 0., 5, 2000, 0, 20.);
    fh_ErelB_vs_theta26_bg->GetXaxis()->SetTitle("angle / deg");
    fh_ErelB_vs_theta26_bg->GetYaxis()->SetTitle("Erel / MeV");

    fh_erelB_vs_ptransHe_bg =
        new TH2F("erelB_vs_ptransHe_bg", "ErelB vs. ptransHe below grazing", 2000, -1000, 1000, 150, 0., 15);
    fh_erelB_vs_ptransHe_bg->GetXaxis()->SetTitle("ptrans / MeV/c");
    fh_erelB_vs_ptransHe_bg->GetYaxis()->SetTitle("Erel / MeV");

    fh_erelB_vs_ptransC_bg =
        new TH2F("erelB_vs_ptransC_bg", "ErelB vs. ptransC below grazing", 2000, -1000, 1000, 150, 0., 15);
    fh_erelB_vs_ptransC_bg->GetXaxis()->SetTitle("ptrans / MeV/c");
    fh_erelB_vs_ptransC_bg->GetYaxis()->SetTitle("Erel / MeV");

    fh_thetaB_4He_cm_bg =
        new TH2F("thetaB_4He_cm_bg", "ThetaB 4He in cm-system below grazing", 200, -1., 1., 2000, 0, 20);
    fh_thetaB_4He_cm_bg->GetXaxis()->SetTitle("cos(theta_cm)");
    fh_thetaB_4He_cm_bg->GetYaxis()->SetTitle("ErelB / MeV");

    fh_phiB_4He_cm_bg = new TH2F("phiB_4He_cm_bg", "phiB 4He in cm-system below grazing", 360, 0., 360., 2000, 0, 20);
    fh_phiB_4He_cm_bg->GetXaxis()->SetTitle("angle / degree");
    fh_phiB_4He_cm_bg->GetYaxis()->SetTitle("Erel / MeV");

    fh_thetaB_12C_cm_bg =
        new TH2F("thetaB_12C_cm_bg", "ThetaB 12C in cm-system below grazing", 200, -1, 1., 2000, 0, 20);
    fh_thetaB_12C_cm_bg->GetXaxis()->SetTitle("cos(theta_cm)");
    fh_thetaB_12C_cm_bg->GetYaxis()->SetTitle("Erel / MeV");

    fh_phiB_12C_cm_bg = new TH2F("phiB_12C_cm_bg", "phiB 12C in cm-system below grazing", 360, 0., 360., 2000, 0, 20);
    fh_phiB_12C_cm_bg->GetXaxis()->SetTitle("angle / degree");
    fh_phiB_12C_cm_bg->GetYaxis()->SetTitle("Erel / MeV");

    fh_ErelB_vs_phibc_ag =
        new TH2F("ErelB_vs_phibccm_ag", "ErelB vs. phi_bc_cm above grazing", 2000, 0, 20., 400, 0., 400.);
    fh_ErelB_vs_phibc_ag->GetYaxis()->SetTitle("phi angle bc / deg");
    fh_ErelB_vs_phibc_ag->GetXaxis()->SetTitle("Erel / MeV");

    fh_ErelB_vs_theta26_ag = new TH2F("ErelB_vs_theta_ag", "ErelB vs. theta above grazing", 125, 0., 5, 2000, 0, 20.);
    fh_ErelB_vs_theta26_ag->GetXaxis()->SetTitle("angle / deg");
    fh_ErelB_vs_theta26_ag->GetYaxis()->SetTitle("Erel / MeV");

    fh_erelB_vs_ptransHe_ag =
        new TH2F("erelB_vs_ptransHe_ag", "ErelB vs. ptransHe above grazing", 2000, -1000, 1000, 150, 0., 15);
    fh_erelB_vs_ptransHe_ag->GetXaxis()->SetTitle("ptrans / MeV/c");
    fh_erelB_vs_ptransHe_ag->GetYaxis()->SetTitle("Erel / MeV");

    fh_erelB_vs_ptransC_ag =
        new TH2F("erelB_vs_ptransC_ag", "ErelB vs. ptransC above grazing", 2000, -1000, 1000, 150, 0., 15);
    fh_erelB_vs_ptransC_ag->GetXaxis()->SetTitle("ptrans / MeV/c");
    fh_erelB_vs_ptransC_ag->GetYaxis()->SetTitle("Erel / MeV");

    fh_thetaB_4He_cm_ag =
        new TH2F("thetaB_4He_cm_ag", "ThetaB 4He in cm-system above grazing", 200, -1., 1., 2000, 0, 20);
    fh_thetaB_4He_cm_ag->GetXaxis()->SetTitle("cos(theta_cm)");
    fh_thetaB_4He_cm_ag->GetYaxis()->SetTitle("ErelB / MeV");

    fh_phiB_4He_cm_ag = new TH2F("phiB_4He_cm_ag", "phiB 4He in cm-system above grazing", 360, 0., 360., 2000, 0, 20);
    fh_phiB_4He_cm_ag->GetXaxis()->SetTitle("angle / degree");
    fh_phiB_4He_cm_ag->GetYaxis()->SetTitle("Erel / MeV");

    fh_thetaB_12C_cm_ag =
        new TH2F("thetaB_12C_cm_ag", "ThetaB 12C in cm-system above grazing", 200, -1, 1., 2000, 0, 20);
    fh_thetaB_12C_cm_ag->GetXaxis()->SetTitle("cos(theta_cm)");
    fh_thetaB_12C_cm_ag->GetYaxis()->SetTitle("Erel / MeV");

    fh_phiB_12C_cm_ag = new TH2F("phiB_12C_cm_ag", "phiB 12C in cm-system above grazing", 360, 0., 360., 2000, 0, 20);
    fh_phiB_12C_cm_ag->GetXaxis()->SetTitle("angle / degree");
    fh_phiB_12C_cm_ag->GetYaxis()->SetTitle("Erel / MeV");

    fh_Erel_vs_phibc_bg =
        new TH2F("Erel_vs_phibccm_bg", "Erel vs. phi_bc_cm below grazing", 2000, 0, 20., 400, 0., 400.);
    fh_Erel_vs_phibc_bg->GetYaxis()->SetTitle("phi angle bc / deg");
    fh_Erel_vs_phibc_bg->GetXaxis()->SetTitle("Erel / MeV");

    fh_Erel_vs_theta26_bg = new TH2F("Erel_vs_theta_bg", "Erel vs. theta below grazing", 125, 0., 5, 2000, 0, 20.);
    fh_Erel_vs_theta26_bg->GetXaxis()->SetTitle("angle / deg");
    fh_Erel_vs_theta26_bg->GetYaxis()->SetTitle("Erel / MeV");

    fh_erel_vs_ptransHe_bg =
        new TH2F("erel_vs_ptransHe_bg", "ErelB vs. ptransHe below grazing", 2000, -1000, 1000, 150, 0., 15);
    fh_erel_vs_ptransHe_bg->GetXaxis()->SetTitle("ptrans / MeV/c");
    fh_erel_vs_ptransHe_bg->GetYaxis()->SetTitle("Erel / MeV");

    fh_erel_vs_ptransC_bg =
        new TH2F("erel_vs_ptransC_bg", "Erel vs. ptransC below grazing", 2000, -1000, 1000, 150, 0., 15);
    fh_erel_vs_ptransC_bg->GetXaxis()->SetTitle("ptrans / MeV/c");
    fh_erel_vs_ptransC_bg->GetYaxis()->SetTitle("Erel / MeV");

    fh_theta_4He_cm_bg = new TH2F("theta_4He_cm_bg", "Theta 4He in cm-system below grazing", 200, -1., 1., 2000, 0, 20);
    fh_theta_4He_cm_bg->GetXaxis()->SetTitle("cos(theta_cm)");
    fh_theta_4He_cm_bg->GetYaxis()->SetTitle("Erel / MeV");

    fh_phi_4He_cm_bg = new TH2F("phi_4He_cm_bg", "phi 4He in cm-system below grazing", 360, 0., 360., 2000, 0, 20);
    fh_phi_4He_cm_bg->GetXaxis()->SetTitle("angle / degree");
    fh_phi_4He_cm_bg->GetYaxis()->SetTitle("Erel / MeV");

    fh_theta_12C_cm_bg = new TH2F("theta_12C_cm_bg", "Theta 12C in cm-system below grazing", 200, -1, 1., 2000, 0, 20);
    fh_theta_12C_cm_bg->GetXaxis()->SetTitle("cos(theta_cm)");
    fh_theta_12C_cm_bg->GetYaxis()->SetTitle("Erel / MeV");

    fh_phi_12C_cm_bg = new TH2F("phi_12C_cm_bg", "phi 12C in cm-system below grazing", 360, 0., 360., 2000, 0, 20);
    fh_phi_12C_cm_bg->GetXaxis()->SetTitle("angle / degree");
    fh_phi_12C_cm_bg->GetYaxis()->SetTitle("Erel / MeV");

    fh_Erel_vs_phibc_ag =
        new TH2F("Erel_vs_phibccm_ag", "Erel vs. phi_bc_cm above grazing", 2000, 0, 20., 400, 0., 400.);
    fh_Erel_vs_phibc_ag->GetYaxis()->SetTitle("phi angle bc / deg");
    fh_Erel_vs_phibc_ag->GetXaxis()->SetTitle("Erel / MeV");

    fh_Erel_vs_theta26_ag = new TH2F("Erel_vs_theta_ag", "Erel vs. theta above grazing", 125, 0., 5, 2000, 0, 20.);
    fh_Erel_vs_theta26_ag->GetXaxis()->SetTitle("angle / deg");
    fh_Erel_vs_theta26_ag->GetYaxis()->SetTitle("Erel / MeV");

    fh_erel_vs_ptransHe_ag =
        new TH2F("erel_vs_ptransHe_ag", "Erel vs. ptransHe above grazing", 2000, -1000, 1000, 150, 0., 15);
    fh_erel_vs_ptransHe_ag->GetXaxis()->SetTitle("ptrans / MeV/c");
    fh_erel_vs_ptransHe_ag->GetYaxis()->SetTitle("Erel / MeV");

    fh_erel_vs_ptransC_ag =
        new TH2F("erel_vs_ptransC_ag", "Erel vs. ptransC above grazing", 2000, -1000, 1000, 150, 0., 15);
    fh_erel_vs_ptransC_ag->GetXaxis()->SetTitle("ptrans / MeV/c");
    fh_erel_vs_ptransC_ag->GetYaxis()->SetTitle("Erel / MeV");

    fh_theta_4He_cm_ag = new TH2F("theta_4He_cm_ag", "Theta 4He in cm-system above grazing", 200, -1., 1., 2000, 0, 20);
    fh_theta_4He_cm_ag->GetXaxis()->SetTitle("cos(theta_cm)");
    fh_theta_4He_cm_ag->GetYaxis()->SetTitle("Erel / MeV");

    fh_phi_4He_cm_ag = new TH2F("phi_4He_cm_ag", "phi 4He in cm-system above grazing", 360, 0., 360., 2000, 0, 20);
    fh_phi_4He_cm_ag->GetXaxis()->SetTitle("angle / degree");
    fh_phi_4He_cm_ag->GetYaxis()->SetTitle("Erel / MeV");

    fh_theta_12C_cm_ag = new TH2F("theta_12C_cm_ag", "Theta 12C in cm-system above grazing", 200, -1, 1., 2000, 0, 20);
    fh_theta_12C_cm_ag->GetXaxis()->SetTitle("cos(theta_cm)");
    fh_theta_12C_cm_ag->GetYaxis()->SetTitle("Erel / MeV");

    fh_phi_12C_cm_ag = new TH2F("phi_12C_cm_ag", "phi 12C in cm-system above grazing", 360, 0., 360., 2000, 0, 20);
    fh_phi_12C_cm_ag->GetXaxis()->SetTitle("angle / degree");
    fh_phi_12C_cm_ag->GetYaxis()->SetTitle("Erel / MeV");

    fh_psum_vs_theta26_nc = new TH2F("psum_vs_theta_nc", "theta vs psum nc", 125, 0., 5, 600, 16000., 19000);
    fh_psum_vs_theta26_nc->GetYaxis()->SetTitle("angle / deg");
    fh_psum_vs_theta26_nc->GetXaxis()->SetTitle("p / MeV/c");

    Int_t nbinsx = 100, nbinsy = 100;
    Double_t binminx = 0, binmaxx = 10, binminy = 0, binmaxy = 10;
    for (Int_t det = 0; det < N_DET_MAX; det++)
    {
        if (det < 2)
        {
            nbinsx = 429;
            binminx = -6.006;
            binmaxx = 6.006;
            nbinsy = 429;
            binminy = -6.006;
            binmaxy = 6.006;
        }
        else if (det > 1 && det < 6)
        {
            nbinsx = 120;
            binminx = -30;
            binmaxx = 30;
            nbinsy = 120;
            binminy = -30;
            binmaxy = 30;
        }
        else
        {
            nbinsx = 89;
            binminx = -60.075;
            binmaxx = 60.075;
            nbinsy = 120;   // 89;
            binminy = -60.; // 075;
            binmaxy = 60.;  // 075;
        }
        fh_psum_vs_xdet[det] = new TH2F(Form("psum_vs_x%s", fDetectorNames[det]),
                                        Form("psum vs x%s", fDetectorNames[det]),
                                        nbinsx,
                                        binminx,
                                        binmaxx,
                                        1500,
                                        16000.,
                                        19000);

        fh_psum_vs_ydet[det] = new TH2F(Form("psum_vs_y%s", fDetectorNames[det]),
                                        Form("psum vs y%s", fDetectorNames[det]),
                                        nbinsy,
                                        binminy,
                                        binmaxy,
                                        1500,
                                        16000.,
                                        19000);

        fh_p_vs_xdet[det] = new TH2F(Form("p_vs_x%s", fDetectorNames[det]),
                                     Form("p vs x%s", fDetectorNames[det]),
                                     nbinsx,
                                     binminx,
                                     binmaxx,
                                     1400,
                                     2000.,
                                     16000);
    }
    fh_yfi23_vs_ytofd_bc = new TH2F("yfi23_vs_ytofd_bc", " yFi23b vs yTofd bc", 120, -60, 60, 429, -6.006, 6.006);
    fh_yfi23_vs_ytofd = new TH2F("yfi23_vs_ytofd", " yFi23b vs yTofd", 120, -60, 60, 429, -6.006, 6.006);
    fh_xfi23_vs_xtofd_bc =
        new TH2F("xfi23_vs_xtofd_bc", " xFi23b vs xTofd bc", 89, -60.075, 60.075, 429, -6.006, 6.006);
    fh_xfi23_vs_xtofd = new TH2F("xfi23_vs_xtofd", " xFi23b vs xTofd", 89, -60.075, 60.075, 429, -6.006, 6.006);
    fh_xy_fib23_bc = new TH2F("xy_fi23_bc", " xFi23b vs yfi23 bc", 429, -6.006, 6.006, 429, -6.006, 6.006);
    fh_xy_fib23 = new TH2F("xy_fi23", " xFi23b vs yfi23b", 429, -6.006, 6.006, 429, -6.006, 6.006);
    fh_xx_fib23 = new TH2F("xx_fi23", " xFi23a C vs xfi23a He", 429, -6.006, 6.006, 429, -6.006, 6.006);
    fh_yy_fib23 = new TH2F("yy_fi23", " yFi23a C vs yfi23a He", 429, -6.006, 6.006, 429, -6.006, 6.006);
    fh_dt_fib23_bc = new TH1F("dt_fib23_bc", " tFi23b - tfi23a bc", 120, -30, 30);

    if (fHitItemsCalifa)
    {
        fh_Erel_withCalifa = new TH1F("Erel_withCalifa", "Erel with is_califa", 2000, 0., 20);
        fh_Erel_withCalifa->GetXaxis()->SetTitle("Erel / MeV");
        fh_Erel_withCalifa->GetYaxis()->SetTitle("counts");

        fh_Erel_withCalifa_2d = new TH2F("Erel_vs_numCrystalCalifa", "Erel vs Numb of Crystal", 40, 0, 40, 200, 0., 20);
        fh_Erel_withCalifa_2d->GetYaxis()->SetTitle("Erel / MeV");
        fh_Erel_withCalifa_2d->GetXaxis()->SetTitle("No crystals in cluster");

        fh_califaenergy_2d =
            new TH2F("Ecalifa_vs_numCrystalCalifa", "Ecalifa vs Numb of Crystal", 40, 0, 40, 1000, 0., 10000.);
        fh_califaenergy_2d->GetYaxis()->SetTitle("Ecalifa / keV");
        fh_califaenergy_2d->GetXaxis()->SetTitle("No crystals in cluster");

        fh_califaenergy_barrel_2d = new TH2F(
            "Ecalifa_barrel_vs_numCrystalCalifa", "Ecalifa_barrel vs Numb of Crystal", 40, 0, 40, 1000, 0., 10000.);
        fh_califaenergy_barrel_2d->GetYaxis()->SetTitle("Ecalifa / keV");
        fh_califaenergy_barrel_2d->GetXaxis()->SetTitle("No crystals in cluster");

        fh_califaenergy_iphos_2d = new TH2F(
            "Ecalifa_iphos_vs_numCrystalCalifa", "Ecalifa_iphos vs Numb of Crystal", 40, 0, 40, 1000, 0., 10000.);
        fh_califaenergy_iphos_2d->GetYaxis()->SetTitle("Ecalifa / keV");
        fh_califaenergy_iphos_2d->GetXaxis()->SetTitle("No crystals in cluster");

        fh_esumcluster_dc_2d =
            new TH2F("fh_esumcluster_dc_2d", "Ecalifa vs Numb of Crystal", 40, 0, 40, 1000, 0., 10000.);
        fh_esumcluster_dc_2d->GetYaxis()->SetTitle("Ecalifa / keV");
        fh_esumcluster_dc_2d->GetXaxis()->SetTitle("No crystals in cluster");

        fh_Erel_vs_theta16O_withcalifa = new TH2F("Erel_vs_theta16O_withcalifa",
                                                  "Erel vs. theta 16O* with califa cut on Ecalifa&tof",
                                                  125,
                                                  0.,
                                                  5,
                                                  200,
                                                  0,
                                                  20.);
        fh_Erel_vs_theta16O_withcalifa->GetXaxis()->SetTitle("angle / deg");
        fh_Erel_vs_theta16O_withcalifa->GetYaxis()->SetTitle("Erel / MeV");

        fh_califa_hitenergy_ag = new TH2F(
            "califa_hit_energy_vs_erel_ag", "Califa hitE vs Erel above grazing", 2000, 0, 20, 1000, 0., 10000.);
        fh_califa_hitenergy_ag->GetYaxis()->SetTitle("Califa energy / keV");
        fh_califa_hitenergy_ag->GetXaxis()->SetTitle("Erel / MeV");

        fh_califa_hitenergy_bg = new TH2F(
            "califa_hit_energy_vs_erel_bg", "Califa hitE vs Erel below grazing", 2000, 0, 20, 1000, 0., 10000.);
        fh_califa_hitenergy_bg->GetYaxis()->SetTitle("Califa energy / keV");
        fh_califa_hitenergy_bg->GetXaxis()->SetTitle("Erel / MeV");

        fh_califa_hitenergy_rand_ag = new TH2F("califa_hit_energy_vs_erel_rand_ag",
                                               "Califa hitE vs Erel above grazing rand coinc",
                                               2000,
                                               0,
                                               20,
                                               1000,
                                               0.,
                                               10000.);
        fh_califa_hitenergy_rand_ag->GetYaxis()->SetTitle("Califa energy / keV");
        fh_califa_hitenergy_rand_ag->GetXaxis()->SetTitle("Erel / MeV");

        fh_califa_hitenergy_rand_bg = new TH2F("califa_hit_energy_vs_erel_rand_bg",
                                               "Califa hitE vs Erel below grazing rand coinc",
                                               2000,
                                               0,
                                               20,
                                               1000,
                                               0.,
                                               10000.);
        fh_califa_hitenergy_rand_bg->GetYaxis()->SetTitle("Califa energy / keV");
        fh_califa_hitenergy_rand_bg->GetXaxis()->SetTitle("Erel / MeV");

        fh_califa_hitenergy_boost =
            new TH1F("hit_califa_energy_DC", "Califa hitE is_track true DopCorr", 1000, 0., 10000.);
        fh_califa_hitenergy_boost->GetYaxis()->SetTitle("Counts");
        fh_califa_hitenergy_boost->GetXaxis()->SetTitle("Energy / keV");

        fh_califa_barrel_hitenergy_boost =
            new TH1F("hit_califa_barrel_energy_DC", "Califa_barrel hitE is_track true DopCorr", 1000, 0., 10000.);
        fh_califa_barrel_hitenergy_boost->GetYaxis()->SetTitle("Counts");
        fh_califa_barrel_hitenergy_boost->GetXaxis()->SetTitle("Energy / keV");

        fh_califa_iphos_hitenergy_boost =
            new TH1F("hit_califa_iphos_energy_DC", "Califa_iphos hitE is_track true DopCorr", 1000, 0., 10000.);
        fh_califa_iphos_hitenergy_boost->GetYaxis()->SetTitle("Counts");
        fh_califa_iphos_hitenergy_boost->GetXaxis()->SetTitle("Energy / keV");

        fh_califa_hitenergycorr_boost =
            new TH1F("hit_califa_energycorr_DC", "Califa hitE reodered is_track true DopCorr", 1000, 0., 10000.);
        fh_califa_hitenergycorr_boost->GetYaxis()->SetTitle("Counts");
        fh_califa_hitenergycorr_boost->GetXaxis()->SetTitle("Energy / keV");

        fh_califa_hitenergy_boost_rand =
            new TH1F("hit_califa_energy_rand", "Califa hitE is_track true DopCorr", 1000, 0., 10000.);
        fh_califa_hitenergy_boost_rand->GetYaxis()->SetTitle("Counts");
        fh_califa_hitenergy_boost_rand->GetXaxis()->SetTitle("Energy / keV");

        fh_Erel_vs_theta16O_withcalifa_rand = new TH2F("Erel_vs_theta16O_withcalifa_rand",
                                                       "Erel vs. theta 16O* with califa cut on Ecalifa&tof",
                                                       125,
                                                       0.,
                                                       5,
                                                       200,
                                                       0,
                                                       20.);
        fh_Erel_vs_theta16O_withcalifa_rand->GetXaxis()->SetTitle("angle / deg");
        fh_Erel_vs_theta16O_withcalifa_rand->GetYaxis()->SetTitle("Erel / MeV");

        fh_Erel_withCalifa_tof = new TH2F("Erel_vs_tofCalifa", "Erei vs tof Califa", 1500, 1000., 4000, 200, 0., 20.);
        fh_Erel_withCalifa_tof->GetYaxis()->SetTitle("Erel / MeV");
        fh_Erel_withCalifa_tof->GetXaxis()->SetTitle("tof / ns");

        fh_Erel_withCalifa_motherId =
            new TH2F("Ecalifa_vs_motherId", "Ecalifa vs motherId Califa", 1600, 900., 2500, 100, 0., 10000);
        fh_Erel_withCalifa_motherId->GetYaxis()->SetTitle("Ecalifa / keV");
        fh_Erel_withCalifa_motherId->GetXaxis()->SetTitle("motherId");

        fh_phi_vs_ecalifa_barrel =
            new TH2F("Phi_vs_Ecalifa_barrel", "Califa phi vs energy barrel", 100, 0., 10000, 360, -180, 180);
        fh_phi_vs_ecalifa_barrel->GetYaxis()->SetTitle("phi Califa / deg");
        fh_phi_vs_ecalifa_barrel->GetXaxis()->SetTitle("Energy Califa / keV");

        fh_phi_vs_ecalifa_iphos =
            new TH2F("Phi_vs_Ecalifa_iphos", "Califa phi vs energy iphos", 100, 0., 10000, 360, -180, 180);
        fh_phi_vs_ecalifa_iphos->GetYaxis()->SetTitle("phi Califa / deg");
        fh_phi_vs_ecalifa_iphos->GetXaxis()->SetTitle("Energy Califa / keV");

        fh_califa_hitenergy_select = new TH1F("ErelEcalifaSum", "Erel+Ecalifa", 2000, 0., 20.);
        fh_califa_hitenergy_select->GetYaxis()->SetTitle("Counts");
        fh_califa_hitenergy_select->GetXaxis()->SetTitle("Energy / keV");

        fh_califa_tofd =
            new TH2F("fh_califa_tof_vs_motherId", "Califa tof vs motherId", 1500, 1000., 4000., 1600, 900., 2500);
        fh_califa_tofd->GetYaxis()->SetTitle("motherId");
        fh_califa_tofd->GetXaxis()->SetTitle("tof / ns");

        fh_califa_tofd_cut =
            new TH2F("fh_califa_tofd_cut", "Califa tof vs Califa energy", 100, 0., 10000, 1500, 1000., 4000);
        fh_califa_tofd_cut->GetYaxis()->SetTitle("tof / ns ");
        fh_califa_tofd_cut->GetXaxis()->SetTitle("Ecalifa / keV");

        fh_califa_tofd_ecr = new TH2F(
            "fh_califa_tofdr_ecr", "Califa crystal tof vs Califa crystal energy", 100, 0, 10000, 1500, 1000., 4000);
        fh_califa_tofd_ecr->GetYaxis()->SetTitle("tof / ns ");
        fh_califa_tofd_ecr->GetXaxis()->SetTitle("ecr / keV");

        fh_califa_ecl = new TH2F(
            "fh_califa_ecl", "Califa cluster energy vs Califa sorted cluster energy", 100, 0., 10000, 100, 0, 10000);
        fh_califa_ecl->GetYaxis()->SetTitle("Cluster energy / keV ");
        fh_califa_ecl->GetXaxis()->SetTitle("Sorted cluster energy / keV");

        fh_crystalNb = new TH2F("fh_crystalNb", "Califa crystalNbnew vs Califa crystalNb", 40, 0., 40, 40, 0, 40);
        fh_crystalNb->GetYaxis()->SetTitle("CrystalNb sorted ");
        fh_crystalNb->GetXaxis()->SetTitle("CrystalNb");

        fh_califa_dtof_ecr = new TH2F(
            "fh_califa_dtof_ecr", "Califa crystal tof vs Califa crystal energy", 100, 0., 10000, 1500, 1000., 4000);
        fh_califa_dtof_ecr->GetYaxis()->SetTitle("tof / ns ");
        fh_califa_dtof_ecr->GetXaxis()->SetTitle("Ecr / keV");

        fh_califa_time_cId = new TH2F(
            "fh_califa_time_cId", "Califa crystal time vs Califa crystal Id", 1600, 900., 2500, 1500, 1000., 4000);
        fh_califa_time_cId->GetYaxis()->SetTitle("time / ns ");
        fh_califa_time_cId->GetXaxis()->SetTitle("crystal Id");

        fh_califa_dtime_ecr = new TH2F(
            "fh_califa_dtime_ecr", "Califa crystal dt vs Califa crystal energy", 100, 0., 10000, 2000, -2000., 2000);
        fh_califa_dtime_ecr->GetYaxis()->SetTitle("dtime / ns ");
        fh_califa_dtime_ecr->GetXaxis()->SetTitle("Ecalifa / keV");
    }
    R3BLOG(debug, "Finished with Init");

    return kSUCCESS;
}

void R3BGlobalAnalysisS494::Exec(Option_t* option)
{
    R3BLOG(debug, "Entering Exec");

    if (fNEvents / 100000. == (int)fNEvents / 100000)
        std::cout << "\rEvents: " << fNEvents << " / " << maxevent << " (" << (int)(fNEvents * 100. / maxevent)
                  << " %) " << std::flush;

    // cout << "New event " << fNEvents << endl;

    FairRootManager* mgr = FairRootManager::Instance();
    if (NULL == mgr)
        LOG(fatal) << "FairRootManager not found";

    if (header)
    {
        timeTS = header->GetTimeStamp();
        eventTS = header->GetEventno();
    }
    //  cout<<"header time eventnum and time: "<<eventTS<<", "<<timeTS<<endl;

    pHex = -100.;
    pHey = -100.;
    pHez = -100.;
    pCx = -100.;
    pCy = -100.;
    pCz = -100.;
    XHe = -100;
    YHe = -100;
    ZHe = -100;
    XC = -100;
    YC = -100;
    ZC = -100;
    mtrackHe = -100;
    mtrackC = -100;
    pHex_mc = -100.;
    pHey_mc = -100.;
    pHez_mc = -100.;
    pCx_mc = -100.;
    pCy_mc = -100.;
    pCz_mc = -100.;
    XHe_mc = -100;
    YHe_mc = -100;
    ZHe_mc = -100;
    XC_mc = -100;
    YC_mc = -100;
    ZC_mc = -100;
    Double_t px, py, pz;
    Double_t theta_16O = 0., theta_26 = 0.;
    Double_t costh26 = 0.;
    Double_t chiHex = 1.e+36, chiHey = 1.e+36, chiCx = 1.e+36, chiCy = 1.e+36;
    Double_t chix = 1.e+36, chiy = 1.e+36;
    pbeam.SetXYZ(0, 0, 17391.5);

    Double_t beta_beam = 0.7593209;
    Double_t gamma_beam = 1.536768;
    Double_t m3He = 2809.41328, mHe, mC, mO, amu;
    if (!fSimu) // these are atima values in MeV/c2
    {
        mHe = 3.7273791e3;
        mC = 11.174862e3;
        mO = 14.89507961548e3;
        amu = 0.93149403e3;
    }
    else // these are geant3/4 values in MeV/c2
    {
        mHe = 3.7273791e3;
        mC = 11.17486e3;
        mO = 14.89508e3;
        amu = 0.931494028e3;
    }

    Bool_t is_alpha = false;
    Bool_t is_carbon = false;
    Bool_t is_oxygen = false;
    Bool_t is_tracked = false;

    Int_t nTracksFound = 0;
    Int_t nHitsTrack = 0;
    Int_t nHitsMCTrack = 0;

    vector<Double_t> posdetHe;
    vector<Double_t> posdetC;
    vector<Double_t> posdetO;

    Double_t tHetofd = -10000., tCtofd = -10000., tOtofd = -10000.;

    if (fTrack)
    {

        // read in Track parameter
        nHitsTrack = fTrack->GetEntriesFast();

        //  if(nHitsTrack < 1) return;

        // cout << "Track hits: " << nHitsTrack << endl;
        for (Int_t l = 0; l < nHitsTrack; l++)
        {
            R3BTrackData* aTrack = (R3BTrackData*)fTrack->At(l);

            if (aTrack->GetQ() == 2)
            {
                // 4He
                XHe = aTrack->GetX(); // cm
                YHe = aTrack->GetY(); // cm
                ZHe = aTrack->GetZ(); // cm

                pHex = aTrack->GetPx(); //  + ((std::rand() / (float)RAND_MAX) - 0.5)/100.* aTrack->GetPx(); // MeV/c
                pHey = aTrack->GetPy(); // + ((std::rand() / (float)RAND_MAX) - 0.5)/100.* aTrack->GetPx(); // MeV/c
                pHez = aTrack->GetPz(); // MeV/c
                mtrackHe = aTrack->GetAoZ() * aTrack->GetQ();

                chiHex = aTrack->GetChix();
                chiHey = aTrack->GetChiy();
                posdetHe = aTrack->GetDetPos();
                tHetofd = aTrack->GetTime();

                is_alpha = true;
                alpha.SetPxPyPzE(pHex, pHey, pHez, sqrt(pow(pHex, 2) + pow(pHey, 2) + pow(pHez, 2) + pow(mHe, 2)));
                if (mtrackHe > 2.4 && mtrackHe < 3.6)
                    helium3.SetPxPyPzE(
                        pHex, pHey, pHez, sqrt(pow(pHex, 2) + pow(pHey, 2) + pow(pHez, 2) + pow(m3He, 2)));

                LOG(debug) << "******************************************" << endl;
                LOG(debug) << "Track In 4He"
                           << "x " << XHe << " y " << YHe << " z " << ZHe << endl;
                LOG(debug) << "px " << pHex << " py " << pHey << " z " << pHez << endl;
                LOG(debug) << "chiHex " << chiHex << ", " << chiHey << endl;
            }
            if (aTrack->GetQ() == 6)
            {
                // 12C
                XC = aTrack->GetX();
                YC = aTrack->GetY();
                ZC = aTrack->GetZ();
                mtrackC = aTrack->GetAoZ() * aTrack->GetQ();

                pCx = aTrack->GetPx(); // + ((std::rand() / (float)RAND_MAX) - 0.5)/100.* aTrack->GetPx();
                pCy = aTrack->GetPy(); // + ((std::rand() / (float)RAND_MAX) - 0.5)/100.* aTrack->GetPx();
                pCz = aTrack->GetPz();

                chiCx = aTrack->GetChix();
                chiCy = aTrack->GetChiy();

                posdetC = aTrack->GetDetPos();
                tCtofd = aTrack->GetTime();

                is_carbon = true;
                carbon.SetPxPyPzE(pCx, pCy, pCz, sqrt(pow(pCx, 2) + pow(pCy, 2) + pow(pCz, 2) + pow(mC, 2)));

                LOG(debug) << "******************************************" << endl;
                LOG(debug) << "Track In 12C"
                           << "x " << XC << " y " << YC << " z " << ZC << endl;
                LOG(debug) << "px " << pCx << " py " << pCy << " z " << pCz << endl;
                LOG(debug) << "chiCx " << chiCx << ", " << chiCy << endl;
            }
            if (aTrack->GetQ() == 8)
            {
                // 160
                Xf = aTrack->GetX();
                Yf = aTrack->GetY();
                Zf = aTrack->GetZ();

                Pxf = aTrack->GetPx();
                Pyf = aTrack->GetPy();
                Pzf = aTrack->GetPz();
                is_oxygen = true;
                Pf_tot = sqrt((Pxf * Pxf) + (Pyf * Pyf) + (Pzf * Pzf));

                chix = aTrack->GetChix();
                chiy = aTrack->GetChiy();

                posdetO = aTrack->GetDetPos();
                tOtofd = aTrack->GetTime();

                LOG(debug) << "******************************************" << endl;
                LOG(debug) << "Track In 16O"
                           << "x " << Xf << " y " << Yf << " z " << Zf << endl;
                LOG(debug) << "px " << Pxf << " py " << Pyf << " z " << Pzf << endl;
            }
        }
    }
    if (fMCTrack && fSimu)
    {
        // read in Monte Carlo Track parameter

        nHitsMCTrack = fMCTrack->GetEntriesFast();
        //  cout << "MCTrack hits: " << nHitsMCTrack << endl;

        for (Int_t l = 0; l < nHitsMCTrack; l++)
        {
            //  cout << "Original MC Data is analyzed" << endl;
            R3BMCTrack* aTrack = (R3BMCTrack*)fMCTrack->At(l);

            Int_t PID = aTrack->GetPdgCode();
            Int_t mother = aTrack->GetMotherId();
            LOG(debug) << "PID " << PID << endl;
            if (mother < 0)
            {
                if (PID == 1000020040)
                {
                    // 4He
                    XHe_mc = aTrack->GetStartX(); // cm
                    YHe_mc = aTrack->GetStartY(); // cm
                    ZHe_mc = aTrack->GetStartZ(); // cm
                    THe_mc = aTrack->GetStartT(); // cm

                    pHex_mc = aTrack->GetPx() * 1000.; // MeV/c
                    pHey_mc = aTrack->GetPy() * 1000.;
                    pHez_mc = aTrack->GetPz() * 1000.;

                    is_alpha = true;
                    alphaMC.SetPxPyPzE(pHex_mc,
                                       pHey_mc,
                                       pHez_mc,
                                       sqrt(pow(pHex_mc, 2) + pow(pHey_mc, 2) + pow(pHez_mc, 2) + pow(mHe, 2)));

                    paMC = alphaMC.Vect();

                    LOG(debug) << "******************************************" << endl;
                    LOG(debug) << "MC Track In 4He"
                               << "x " << XHe_mc << " y " << YHe_mc << " z " << ZHe_mc << endl;
                    LOG(debug) << "px " << pHex_mc << " py " << pHey_mc << " pz " << pHez_mc << endl;
                }
                if (PID == 1000060120)
                {
                    // 12C
                    XC_mc = aTrack->GetStartX();
                    YC_mc = aTrack->GetStartY();
                    ZC_mc = aTrack->GetStartZ();
                    TC_mc = aTrack->GetStartT();

                    pCx_mc = aTrack->GetPx() * 1000.;
                    pCy_mc = aTrack->GetPy() * 1000.;
                    pCz_mc = aTrack->GetPz() * 1000.;

                    is_carbon = true;
                    carbonMC.SetPxPyPzE(
                        pCx_mc, pCy_mc, pCz_mc, sqrt(pow(pCx_mc, 2) + pow(pCy_mc, 2) + pow(pCz_mc, 2) + pow(mC, 2)));

                    pcMC = carbonMC.Vect();

                    LOG(debug) << "******************************************" << endl;
                    LOG(debug) << "MC Track In 12C"
                               << "x " << XC_mc << " y " << YC_mc << " z " << ZC_mc << endl;
                    LOG(debug) << "px " << pCx_mc << " py " << pCy_mc << " pz " << pCz_mc << endl;
                }
                if (PID == 1000080160)
                {
                    // 16O
                    Xf_mc = aTrack->GetStartX();
                    Yf_mc = aTrack->GetStartY();
                    Zf_mc = aTrack->GetStartZ();

                    Pxf_mc = aTrack->GetPx() * 1000.;
                    Pyf_mc = aTrack->GetPy() * 1000.;
                    Pzf_mc = aTrack->GetPz() * 1000.;

                    is_oxygen = true;
                    Pf_tot_mc = sqrt((Pxf_mc * Pxf_mc) + (Pyf_mc * Pyf_mc) + (Pzf_mc * Pzf_mc));

                    LOG(debug) << "******************************************" << endl;
                    LOG(debug) << "MC Track In 16O"
                               << "x " << Xf_mc << " y " << Yf_mc << " z " << Zf_mc << endl;
                    LOG(debug) << "px " << Pxf_mc << " py " << Pyf_mc << " z " << Pzf_mc << endl;
                }
            }
        }
    }

    // new version of analysis
    if (fTrack && is_carbon && is_alpha && fPairs && chiHex < 1.e+36) // 17364.240
    {
        fNeventstrack += 1;

        // (0,1) = (xfi23a,yfi23a), .. (12,13) = (xtofd,ytofd)
        fh_xy_fib23_bc->Fill(posdetHe.at(0), posdetHe.at(3));
        fh_xy_fib23_bc->Fill(posdetC.at(0), posdetC.at(3));
        fh_xfi23_vs_xtofd_bc->Fill(posdetHe.at(12), posdetHe.at(0));
        fh_xfi23_vs_xtofd_bc->Fill(posdetC.at(12), posdetC.at(0));
        fh_yfi23_vs_ytofd_bc->Fill(posdetHe.at(13), posdetHe.at(3));
        fh_yfi23_vs_ytofd_bc->Fill(posdetC.at(13), posdetC.at(3));

        fh_target_xx->Fill(XC, posdetC.at(0));
        fh_target_yy->Fill(YC, posdetC.at(3));
        fh_target_xx->Fill(XHe, posdetHe.at(0));
        fh_target_yy->Fill(YHe, posdetHe.at(3));

        Double_t sq_chi = sqrt(chiHex * chiHex + chiCx * chiCx);
        fh_chi2->Fill(sq_chi);
        fh_chiy_vs_chix_nc->Fill(chiHex, chiCx);

        LOG(debug) << "Entering Pair analysis***" << endl;
        // if (chiHex < fcut_chiX && chiHey < fcut_chiY && chiCx < fcut_chiX && chiCy < fcut_chiY)
        //		&&((alpha.Px() > 0 && carbon.Px() < 0) || (alpha.Px() < 0 && carbon.Px() > 0)))
        if (alpha.Pz() > 0 && carbon.Pz() > 0)
        {
            //   if(abs(XHe-XC)*100. > 0.03) return;
            //   if(abs(YHe-YC)*100. > 0.03) return;

            //  if(abs(XHe-XC)*100. > 0.2) return;
            //  if(abs(YHe-YC)*100. > 0.2) return;

            pa = alpha.Vect();
            pc = carbon.Vect();

            // Calculate angle between alpha and C
            if (alpha.Pz() == 0 || carbon.Pz() == 0)
                return;
            /*
                        if (alpha.Pz() < 3900 || alpha.Pz() > 4800)
                           return;
                       if (carbon.Pz() < 12600 || carbon.Pz() > 13500)
                           return;
           */

            theta_26 = alpha.Angle(carbon.Vect()) * TMath::RadToDeg(); // angle alpha carbon (theta)
            Double_t phi_26 = (alpha.Phi() - carbon.Phi()) * TMath::RadToDeg();

            // if(theta_26 < 0.5) return;

            /* * Calculate invariant mass and relative energy * */
            m_inva = (alpha + carbon).M(); // invariant mass

            Erel = m_inva - mHe - mC; // relative Energy

            R3BLOG(debug, "Erel has been calculated");

            /*
                        Double_t costh26_test = (alpha.Px() * carbon.Px() + alpha.Py() * carbon.Py() + alpha.Pz() *
               carbon.Pz()) / (sqrt(alpha.Px() * alpha.Px() + alpha.Py() * alpha.Py() + alpha.Pz() * alpha.Pz()) *
                               sqrt(carbon.Px() * carbon.Px() + carbon.Py() * carbon.Py() + carbon.Pz() * carbon.Pz()));

                        Double_t test_theta = acos(costh26_test) * TMath::RadToDeg(); // opening angle
                        Double_t m_inv_test = sqrt(mHe * mHe + mC * mC + 2. * sqrt(pc.Mag() * pc.Mag() + mC * mC) *
                                              sqrt(pa.Mag() * pa.Mag() + mHe * mHe) -
                                              2. * pa.Mag() * pc.Mag() * cos(test_theta * TMath::DegToRad()));

                        cout<<"Checking minv: "<<m_inva<<", "<<m_inv_test<<endl;
            */

            oxygen = alpha + carbon; // excited oxygen
            poexc = oxygen.Vect();

            // cout<<"gamma "<<oxygen.Gamma()<<endl;
            theta_16O = oxygen.Theta() * TMath::RadToDeg();
            Double_t phi_16O = oxygen.Phi() * TMath::RadToDeg();
            if (phi_16O < 0)
                phi_16O += 360.;

            // photon.SetPxPyPzE(oxygen.Px, oxygen.Py, 0., sqrt(pow(oxygen.Px, 2) + pow(oxygen.Py, 2)));

            // transfer to cm system and make some rotations
            // boost them to centre of mass
            alpha_cm = alpha;
            carbon_cm = carbon;
            oxygen_cm = oxygen;
            alpha_cm.Boost(-oxygen.BoostVector());
            carbon_cm.Boost(-oxygen.BoostVector());
            oxygen_cm.Boost(-oxygen.BoostVector());

            pa_cm = alpha_cm.Vect();
            pc_cm = carbon_cm.Vect();

            phi_bc_cm = (oxygen.Phi() - alpha_cm.Phi()) * TMath::RadToDeg();
            if (phi_bc_cm < 0)
                phi_bc_cm += 360.;

            TVector3 norm;
            Double_t normx = (-1.) * pbeam.Z() * poexc.Y();
            Double_t normy = (-1.) * pbeam.Z() * poexc.X();
            norm.SetXYZ(normx, normy, 0); // vector perpendicular to (pbeam,p16O*)
            Double_t thetanorm_alpha = pa_cm.Theta();
            Double_t thetanorm_carbon = pc_cm.Theta();
            Double_t pinplaneHe = pa_cm.Perp(norm); // project to plane defined by pbeam and p16O*
            Double_t pinplaneC = pc_cm.Perp(norm);  // project to plane defined by pbeam and p16O*
            Double_t ptransHe = pinplaneHe * pa_cm.Phi();
            Double_t ptransC = pinplaneC * pc_cm.Phi();

            // cout << " px: " << alpha_cm.Px() + carbon_cm.Px() << "  py:  " << alpha_cm.Py() + carbon_cm.Py()
            //     << "  pz:  " << alpha_cm.Pz() + carbon_cm.Pz() << endl;
            /*
                        cout << "lab " << alpha.Px() << "  " << alpha.Py() << "  " << alpha.Pz() << endl;
                        cout << "lab " << carbon.Px() << "  " << carbon.Py() << "  " << carbon.Pz() << endl;
                        cout << "lab " << oxygen.Px() << "  " << oxygen.Py() << "  " << oxygen.Pz() << endl;

                        cout << "cm " << alpha_cm.Px() << "  " << alpha_cm.Py() << "  " << alpha_cm.Pz() << endl;
                        cout << "cm " << carbon_cm.Px() << "  " << carbon_cm.Py() << "  " << carbon_cm.Pz() << endl;
                        cout << "cm " << oxygen_cm.Px() << "  " << oxygen_cm.Py() << "  " << oxygen_cm.Pz() << endl;
            */

            Double_t tthetaxC = pCx / pCz;    // * 1000.;
            Double_t tthetayC = pCy / pCz;    // * 1000.;
            Double_t tthetaxHe = pHex / pHez; // * 1000.;
            Double_t tthetayHe = pHey / pHez; // * 1000.;

            Double_t thetaxC = atan(tthetaxC) * 1000.;
            Double_t thetayC = atan(tthetayC) * 1000.;
            Double_t thetaxHe = atan(tthetaxHe) * 1000.;
            Double_t thetayHe = atan(tthetayHe) * 1000.;

            Double_t Ea = mHe * sqrt(1. + (pa.Mag() / mHe) * (pa.Mag() / mHe)) / 1000.; // in GeV
            Double_t Ec = mC * sqrt(1. + (pc.Mag() / mC) * (pc.Mag() / mC)) / 1000.;    // in GeV
            /*
                   cout<<setprecision(10)<<"12C einzelne: "<<pCx<<" - "<<carbon.Px()<<"; "<<pCy<<" - "<<carbon.Py()<<";
               "<<pCz<<" - "<<carbon.Pz()<<endl; cout<<setprecision(10)<<"12C total: "<<pc.Mag()<<" -
               "<<sqrt(pCx*pCx+pCy*pCy+pCz*pCz)<<endl; cout<<setprecision(10)<<"4He einzelne: "<<pHex<<" -
               "<<alpha.Px()<<"; "<<pHey<<" - "<<alpha.Py()<<"; "<<pHez<<" - "<<alpha.Pz()<<endl;
                   cout<<setprecision(10)<<"4He total: "<<pa.Mag()<<" - "<<sqrt(pHex*pHex+pHey*pHey+pHez*pHez)<<endl;
                   cout<<setprecision(10)<<"Sum of two: "<<(pa+pc).Mag()<<",
               "<<sqrt((pHex+pCx)*(pHex+pCx)+(pHey+pCy)*(pHey+pCy)+(pHez+pCz)*(pHez+pCz))<<", "<<
                                       sqrt(pa.Mag()*pa.Mag()+pc.Mag()*pc.Mag()+2. * pa.Mag() * pc.Mag() * cos(theta_26
               * TMath::DegToRad()))<<endl;
             */

            fh_target_xy_nc->Fill(XHe, YHe);
            fh_target_xy_nc->Fill(XC, YC);

            fh_p_vs_x0->Fill(XC, pc.Mag());
            fh_p_vs_x0->Fill(XHe, pa.Mag());
            fh_psum_vs_x0->Fill(XC, (pa + pc).Mag());

            fh_p_vs_y0->Fill(YC, pc.Mag());
            fh_p_vs_y0->Fill(YHe, pa.Mag());
            fh_psum_vs_y0->Fill(YC, (pa + pc).Mag());

            fh_mass_nc->Fill(mtrackC, mtrackHe);
            fh_energy_nc->Fill(Ec, Ea);
            fh_px_He_nc->Fill(alpha.Px());
            fh_py_He_nc->Fill(alpha.Py());
            fh_pz_He_nc->Fill(alpha.Pz());
            fh_p_He_nc->Fill(pa.Mag());

            fh_px_C_nc->Fill(carbon.Px());
            fh_py_C_nc->Fill(carbon.Py());
            fh_pz_C_nc->Fill(carbon.Pz());
            fh_p_C_nc->Fill(pc.Mag());

            fh_px_px_nc->Fill(alpha.Px(), carbon.Px());
            fh_py_py_nc->Fill(alpha.Py(), carbon.Py());
            fh_pz_pz_nc->Fill(alpha.Pz(), carbon.Pz());

            //  if((alpha.Px()>100. && carbon.Px()>300.) || (alpha.Px()<-100. && carbon.Px()<-300.)) cout<<"CHECK
            //  EVENTS: "<<fNEvents<<", "<<chiCx<<","<<
            // carbon.Px()<<"; "<<alpha.Px()<<", "<<Erel<<", "<<posdetC.at(0)-posdetHe.at(0)<<",
            // "<<posdetC.at(3)-posdetHe.at(3)<<endl;

            fh_psum_nc->Fill((pa + pc).Mag());
            fh_theta26_nc->Fill(theta_26);
            fh_Erel_nc->Fill(Erel); // relative Energy plot

            fh_psum_vs_event_nc->Fill(fNeventstrack, (pa + pc).Mag());
            fh_psum_vs_chi->Fill(sq_chi, (pa + pc).Mag());

            fh_theta26_vs_chi->Fill(sq_chi, theta_26);
            fh_Erel_vs_chi->Fill(sq_chi, Erel);
            fh_psum_vs_theta26_nc->Fill(theta_26, (pa + pc).Mag());

            if ((alpha.Px() > 0. && carbon.Px() < 0.) || (alpha.Px() < 0. && carbon.Px() > 0.))
            {

                fh_ErelB_nc->Fill(Erel);
            }

            fh_theta_bc_cm_nc->Fill(alpha_cm.Theta() * TMath::RadToDeg());
            fh_phi_bc_cm_nc->Fill(phi_bc_cm);

            fh_chiy_vs_chix->Fill(chiHex, chiCx);

            is_tracked = false;

            if (mtrackHe > 2.5 && mtrackHe < 3.5 && mtrackC > 11.5 && mtrackC < 12.5)
            {
                Double_t m_invaHe3 = (helium3 + carbon).M(); // invariant mass
                Double_t ErelHe3 = m_invaHe3 - m3He - mC;
                Double_t theta_15O = (helium3 + carbon).Theta() * TMath::RadToDeg();
                fh_Erel_vs_theta16O_3He12C->Fill(theta_15O, ErelHe3);
            }

            R3BLOG(debug, "Cuts to be set");

            //  Double_t Eacheck = 17.1514-0.665768*Ec;
            Double_t Eacheck = -75.6469 + 10.2155 * Ec - 0.318913 * Ec * Ec;
            Bool_t becut = false;
            if (fGraphCuts)
            {
                if (cut_EHe_EC->IsInside(Ec, Ea))
                    becut = true;
            }
            else
            {
                if (abs(Ea - Eacheck) < 0.15)
                    becut = true;
            }

            pcMC = carbonMC.Vect();
            if (fSimu)
            {
                fh_psum_vs_dx0->Fill(XHe - XHe_mc, (pa + pc).Mag());
            }
            if (mtrackC > 11.5 && mtrackC < 12.5 && abs(posdetC.at(0)) > fxfibcut && abs(posdetC.at(3)) > fyfibcut)
            {
                if (fSimu)
                    fh_pC_vs_dx0->Fill(XC - XC_mc, carbon.Pz());
                if (fSimu)
                    fh_pC_vs_dx0_mc->Fill(XC - XC_mc, carbonMC.Pz());
                if (fSimu)
                    fh_pC_vs_x0->Fill(XC, carbon.Pz());
                if (fSimu)
                    fh_pC_vs_x0_mc->Fill(XC_mc, carbonMC.Pz());
                if (fSimu)
                    fh_dpC_vs_dx0->Fill(XC - XC_mc, carbon.Pz() - carbonMC.Pz());
                if (fSimu)
                    fh_x0_vs_dx0->Fill(XC_mc, XC);

                for (Int_t idet = 0; idet < N_DET_MAX; idet++)
                {
                    fh_p_vs_xdet[idet]->Fill(posdetC.at(2 * idet), pc.Mag());
                }
            }
            if (mtrackHe > 3.5 && mtrackHe < 4.5 && abs(posdetHe.at(0)) > fxfibcut && abs(posdetHe.at(3)) > fyfibcut)
            {
                if (fSimu)
                    fh_pHe_vs_dx0->Fill(XHe - XHe_mc, alpha.Pz());
                if (fSimu)
                    fh_pHe_vs_dx0_mc->Fill(XHe - XHe_mc, alphaMC.Pz());
                if (fSimu)
                    fh_pHe_vs_x0->Fill(XHe, alpha.Pz());
                if (fSimu)
                    fh_pHe_vs_x0_mc->Fill(XHe_mc, alphaMC.Pz());
                if (fSimu)
                    fh_dpHe_vs_dx0->Fill(XHe - XHe_mc, alpha.Pz() - alphaMC.Pz());

                for (Int_t idet = 0; idet < N_DET_MAX; idet++)
                {
                    fh_psum_vs_xdet[idet]->Fill(posdetHe.at(2 * idet), (pa + pc).Mag());
                    fh_p_vs_xdet[idet]->Fill(posdetHe.at(2 * idet), pa.Mag());
                }
            }

            // Dima's tracker only one chi2;&& becut
            if (sq_chi < (fcut_chiX + fcut_chiY) && mtrackHe > 3.5 && mtrackHe < 4.5 && mtrackC > 11.5 &&
                mtrackC < 12.5) // && becut)
                                //  &&((alpha.Px() > 0. && carbon.Px() < 0.) ||
            // (alpha.Px() < 0. && carbon.Px() > 0.))) // abs(carbon.Pz()+alpha.Pz()-17390.) < 400. )
            {
                R3BLOG(debug, "Good track analysis start");

                fh_phi26_vs_chi->Fill(phi_16O, phi_bc_cm);

                // (0,1) = (xfi23a,yfi23a), .. (12,13) = (xtofd,ytofd)
                fh_xy_fib23->Fill(posdetHe.at(0), posdetHe.at(3));
                fh_xy_fib23->Fill(posdetC.at(0), posdetC.at(3));
                fh_xx_fib23->Fill(posdetHe.at(0), posdetC.at(0));
                fh_yy_fib23->Fill(posdetHe.at(3), posdetC.at(3));
                fh_xfi23_vs_xtofd->Fill(posdetHe.at(12), posdetHe.at(0));
                fh_xfi23_vs_xtofd->Fill(posdetC.at(12), posdetC.at(0));
                fh_yfi23_vs_ytofd->Fill(posdetHe.at(13), posdetHe.at(3));
                fh_yfi23_vs_ytofd->Fill(posdetC.at(13), posdetC.at(3));

                for (Int_t idet = 0; idet < N_DET_MAX; idet++)
                {
                    // fh_psum_vs_xdet[idet]->Fill(posdetC.at(2 * idet), (pa + pc).Mag());
                    // fh_psum_vs_xdet[idet]->Fill(posdetHe.at(2 * idet), (pa + pc).Mag());
                    fh_psum_vs_ydet[idet]->Fill(posdetC.at(2 * idet + 1), (pa + pc).Mag());
                    fh_psum_vs_ydet[idet]->Fill(posdetHe.at(2 * idet + 1), (pa + pc).Mag());
                    // fh_p_vs_xdet[idet]->Fill(posdetC.at(2 * idet), pc.Mag());
                    //  fh_p_vs_xdet[idet]->Fill(posdetHe.at(2 * idet), pa.Mag());
                }

                Double_t costheta_26_max =
                    (mHe * mHe + mC * mC +
                     2. * sqrt(13043.625 * 13043.625 + mC * mC) * sqrt(4347.875 * 4347.875 + mHe * mHe) -
                     (Erel + mHe + mC) * (Erel + mHe + mC)) /
                    (2. * 13043.625 * 4347.875);
                // With tofd positions
                Double_t dr_tofd = sqrt((posdetC.at(12) - posdetHe.at(12)) * (posdetC.at(12) - posdetHe.at(12)) +
                                        (posdetC.at(13) - posdetHe.at(13)) * (posdetC.at(13) - posdetHe.at(13)));
                fh_Erel_vs_dr_tofd->Fill(dr_tofd, Erel);
                Double_t theta_26_max_tofd = acos(costheta_26_max) * TMath::RadToDeg() + 0.7;
                Double_t dr_limit_tofd = 2. * 686.1278 * sin(0.5 * theta_26_max_tofd * TMath::DegToRad());
                fh_Erel_vs_dr_limit_tofd->Fill(dr_limit_tofd, Erel);
                if (dr_tofd < dr_limit_tofd)
                    fh_Erel_limit_tofd->Fill(Erel);

                // with fi23 positions
                Double_t dr_fi23 = sqrt((posdetC.at(0) - posdetHe.at(0)) * (posdetC.at(0) - posdetHe.at(0)) +
                                        (posdetC.at(3) - posdetHe.at(3)) * (posdetC.at(3) - posdetHe.at(3)));
                fh_Erel_vs_dr_fi23->Fill(dr_fi23, Erel);
                Double_t theta_26_max_fi23 = acos(costheta_26_max) * TMath::RadToDeg() + 0.3;
                Double_t dr_limit_fi23 = 2. * 91.02052 * sin(0.5 * theta_26_max_fi23 * TMath::DegToRad());
                fh_Erel_vs_dr_limit_fi23->Fill(dr_limit_fi23, Erel);
                if (dr_fi23 < dr_limit_fi23)
                    fh_Erel_limit_fi23->Fill(Erel);

                fh_dr_vs_dr->Fill(dr_fi23, dr_tofd);
                if (dr_tofd < dr_limit_tofd)
                    fh_dr_vs_dr_limit->Fill(dr_fi23, dr_tofd);

                R3BLOG(debug, "Starting Califa analysis");

                Bool_t califaHit_true = false;
                if (fHitItemsCalifa)
                {
                    // CALIFA
                    auto detCalifa = fHitItemsCalifa;
                    Int_t nHitsCalifaH = detCalifa->GetEntriesFast();

                    for (Int_t ihitC = 0; ihitC < nHitsCalifaH; ihitC++)
                    {
                        R3BCalifaClusterData* hitCalifa = (R3BCalifaClusterData*)detCalifa->At(ihitC);
                        if (!hitCalifa)
                            continue;

                        ULong64_t timeCalifa = hitCalifa->GetTime();
                        Double_t timerelCalifa = (double)(timeCalifa - timeTS);
                        Double_t tof = timerelCalifa - tCtofd;
                        Int_t CrystalNb = hitCalifa->GetCrystalList().size();
                        Double_t Energy = hitCalifa->GetEnergy();
                        Double_t theta = hitCalifa->GetTheta();
                        Double_t Energy_dc = Energy * GAMMA * (1 - BETA * TMath::Cos(theta));
                        Int_t motherId = hitCalifa->GetMotherCrystal();
                        std::vector<Int_t> clist = hitCalifa->GetCrystalList();
                        std::vector<Double_t> elist = hitCalifa->GetEnergyList();
                        std::vector<ULong64_t> tlist = hitCalifa->GetTimeList();

                        /*
                            Double_t GetEnergy() const { return fEnergy; }
                            Double_t GetNf() const { return fNf; }
                            Double_t GetNs() const { return fNs; }
                            Double_t GetTheta() const { return fTheta; }
                            Double_t GetPhi() const { return fPhi; }
                            ULong64_t GetTime() const { return fTime; }
                            uint32_t GetClusterType() const { return fClusterType; }
                            Int_t GetNbOfCrystalHits()  const { return fCrystalList.size(); }
                            std::vector<Int_t> GetCrystalList() const {return fCrystalList; }
                            Int_t GetMotherCrystal() const {return fCrystalList.at(0); }
                            std::vector<Double_t> GetEnergyList() const {return fEnergyList; }
                            std::vector<ULong64_t> GetTimeList() const {return fTimeList; }
                       */

                        if (hitCalifa->GetClusterType() == 1 &&
                            !(IS_NAN(Energy))) // gammas (IS_NAN(Energy) == overflow)
                        {
                            Double_t esumcluster = 0;
                            Int_t CrystalNb_new = 0;
                            Double_t time0 = 0;
                            for (Int_t icr = 0; icr < CrystalNb; icr++)
                            {
                                Double_t tcr = (double)(tlist.at(icr) - timeTS);
                                Double_t tofcr = tcr - tCtofd;
                                Int_t crystalId = clist.at(icr);

                                if (crystalId > 927 && crystalId < 1953) // barrel
                                {
                                    for (Int_t irange = 0; irange < 17; irange++)
                                    {
                                        Int_t low = 897 + irange * 64;
                                        Int_t high = 960 + irange * 64;
                                        if (crystalId >= low && crystalId <= high)
                                        {
                                            if (irange % 2 == 0)
                                            {
                                                tofcr = tofcr - tofCalifa_shift;
                                            }
                                            else
                                            {
                                                tofcr = tofcr + tofCalifa_shift;
                                            }

                                            break;
                                        }
                                    }
                                }
                                else // iphos
                                {
                                    if (crystalId >= 1953 && crystalId <= 1984)
                                        tofcr = tofcr - 118.8;
                                    if (crystalId >= 1985 && crystalId <= 2048)
                                        tofcr = tofcr + 118.8;
                                    if (crystalId >= 2049 && crystalId <= 2112)
                                        tofcr = tofcr - 118.8;
                                    if (crystalId >= 2113 && crystalId <= 2176)
                                        tofcr = tofcr + 118.8;
                                    if (crystalId >= 2177 && crystalId <= 2240)
                                        tofcr = tofcr - 118.8;
                                    if (crystalId >= 2241 && crystalId <= 2304)
                                        tofcr = tofcr + 118.8;
                                    if (crystalId >= 2305 && crystalId <= 2343)
                                        tofcr = tofcr - 118.8;
                                    if (crystalId >= 2344 && crystalId <= 2355)
                                        tofcr = tofcr + 118.8;
                                    if (crystalId >= 2356 && crystalId <= 2367)
                                        tofcr = tofcr - 118.8;
                                    if (crystalId >= 2368 && crystalId <= 2379)
                                        tofcr = tofcr + 118.8;
                                    if (crystalId >= 2380 && crystalId <= 2391)
                                        tofcr = tofcr - 118.8;
                                    if (crystalId >= 2392 && crystalId <= 2403)
                                        tofcr = tofcr + 118.8;
                                    if (crystalId >= 2404 && crystalId <= 2415)
                                        tofcr = tofcr - 118.8;
                                    if (crystalId >= 2416 && crystalId <= 2427)
                                        tofcr = tofcr + 118.8;
                                    if (crystalId >= 2428 && crystalId <= 2433)
                                        tofcr = tofcr - 118.8;
                                }

                                Double_t ecr = elist.at(icr);
                                if (icr == 0)
                                    time0 = tofcr;
                                fh_califa_dtof_ecr->Fill(ecr, tcr);
                                fh_califa_time_cId->Fill(crystalId, tofcr + tCtofd);

                                // if(icr > 0) fh_califa_dtof_ecr->Fill(ecr,tofcr);
                                if (icr > 0 && elist.at(0) < 6000.)
                                    fh_califa_dtime_ecr->Fill(ecr, timerelCalifa - tcr);
                                // if (ecr>=400.) {
                                if (cut_CalifaTof->IsInside(ecr, timerelCalifa - tcr))
                                {
                                    // if (std::abs(timerelCalifa-tcr) < 400.) {
                                    esumcluster = esumcluster + ecr;
                                    CrystalNb_new += 1;

                                    if (Energy_dc > 3200 && Energy_dc < 5000)
                                        fh_califa_tofd_ecr->Fill(ecr, tofcr);
                                }
                            }
                            fh_califa_ecl->Fill(elist.at(0), Energy_dc);
                            // cout<<"Energies: "<<Energy<<", "<<	Energy_dc<<", "<<elist.at(0)<<endl;

                            Double_t esumcluster_dc = esumcluster * GAMMA * (1 - BETA * TMath::Cos(theta));

                            if (motherId > 927 && motherId < 1953) // barrel
                            {
                                for (Int_t irange = 0; irange < 17; irange++)
                                {
                                    Int_t low = 897 + irange * 64;
                                    Int_t high = 960 + irange * 64;
                                    if (motherId >= low && motherId <= high)
                                    {
                                        if (irange % 2 == 0)
                                        {
                                            tof = tof - tofCalifa_shift;
                                        }
                                        else
                                        {
                                            tof = tof + tofCalifa_shift;
                                        }

                                        break;
                                    }
                                }
                            }
                            else // iphos
                            {
                                if (motherId >= 1953 && motherId <= 1984)
                                    tof = tof - 118.8;
                                if (motherId >= 1985 && motherId <= 2048)
                                    tof = tof + 118.8;
                                if (motherId >= 2049 && motherId <= 2112)
                                    tof = tof - 118.8;
                                if (motherId >= 2113 && motherId <= 2176)
                                    tof = tof + 118.8;
                                if (motherId >= 2177 && motherId <= 2240)
                                    tof = tof - 118.8;
                                if (motherId >= 2241 && motherId <= 2304)
                                    tof = tof + 118.8;
                                if (motherId >= 2305 && motherId <= 2343)
                                    tof = tof - 118.8;
                                if (motherId >= 2344 && motherId <= 2355)
                                    tof = tof + 118.8;
                                if (motherId >= 2356 && motherId <= 2367)
                                    tof = tof - 118.8;
                                if (motherId >= 2368 && motherId <= 2379)
                                    tof = tof + 118.8;
                                if (motherId >= 2380 && motherId <= 2391)
                                    tof = tof - 118.8;
                                if (motherId >= 2392 && motherId <= 2403)
                                    tof = tof + 118.8;
                                if (motherId >= 2404 && motherId <= 2415)
                                    tof = tof - 118.8;
                                if (motherId >= 2416 && motherId <= 2427)
                                    tof = tof + 118.8;
                                if (motherId >= 2428 && motherId <= 2433)
                                    tof = tof - 118.8;
                            }

                            if (Energy_dc > 3200 && Energy_dc < 5000)
                            {
                                fh_Erel_withCalifa_tof->Fill(tof, Erel);
                                fh_califa_tofd->Fill(tof, motherId);
                            }

                            if (motherId > 927 && motherId < 2433) // barrel+iphos
                            {
                                fh_califa_tofd_cut->Fill(Energy_dc, tof);
                                if (tof > 2130 && tof < 2470)
                                {
                                    // fh_califa_ecl->Fill(esumcluster_dc,Energy_dc);

                                    fh_esumcluster_dc_2d->Fill(CrystalNb_new, esumcluster_dc);
                                    if ((motherId < 1953 && CrystalNb_new > 1 && CrystalNb_new < 5) ||
                                        (motherId > 1952 && CrystalNb_new > 1 && CrystalNb_new < 5))
                                    {
                                        fh_califa_hitenergycorr_boost->Fill(esumcluster_dc);
                                    }
                                    fh_crystalNb->Fill(CrystalNb, CrystalNb_new);
                                    fh_Erel_withCalifa_2d->Fill(CrystalNb, Erel);
                                    fh_califaenergy_2d->Fill(CrystalNb, Energy_dc);
                                    if (motherId < 1953)
                                        fh_califaenergy_barrel_2d->Fill(CrystalNb, Energy_dc);
                                    if (motherId > 1952)
                                        fh_califaenergy_iphos_2d->Fill(CrystalNb, Energy_dc);
                                    if (motherId < 1953 && CrystalNb > 1 && CrystalNb < 5)
                                        fh_phi_vs_ecalifa_barrel->Fill(Energy_dc,
                                                                       hitCalifa->GetPhi() * TMath::RadToDeg());
                                    if (motherId > 1952 && CrystalNb > 1 && CrystalNb < 5)
                                        fh_phi_vs_ecalifa_iphos->Fill(Energy_dc,
                                                                      hitCalifa->GetPhi() * TMath::RadToDeg());

                                    if ((motherId < 1953 && CrystalNb > 1 && CrystalNb < 5) ||
                                        (motherId > 1952 && CrystalNb > 1 && CrystalNb < 5))
                                    {
                                        fh_Erel_withCalifa_motherId->Fill(motherId, Energy_dc);
                                        fh_califa_hitenergy_boost->Fill(Energy_dc);
                                        if (motherId < 1953)
                                            fh_califa_barrel_hitenergy_boost->Fill(Energy_dc);
                                        if (motherId > 1952)
                                            fh_califa_iphos_hitenergy_boost->Fill(Energy_dc);
                                        fh_Erel_withCalifa->Fill(Erel);
                                        if (Energy_dc > 3200 && Energy_dc < 5000)
                                        {
                                            fh_Erel_vs_theta16O_withcalifa->Fill(theta_16O, Erel);
                                            fh_califa_hitenergy_select->Fill(Energy_dc * 1.e-3 + Erel);
                                        }
                                        if (theta_16O < fThetaGrazing)
                                            fh_califa_hitenergy_bg->Fill(Erel, Energy_dc);
                                        else
                                            fh_califa_hitenergy_ag->Fill(Erel, Energy_dc);
                                    }
                                }
                                // else if((tof < 1967 && tof > 1665 ) || (tof > 2656 && tof < 2958))
                                else if (tof < 1967)
                                {
                                    fh_califa_hitenergy_boost_rand->Fill(Energy_dc);

                                    if (Energy_dc > 3200 && Energy_dc < 5000)
                                    {
                                        fh_Erel_vs_theta16O_withcalifa_rand->Fill(theta_16O, Erel);
                                    }
                                    if (theta_16O < fThetaGrazing)
                                        fh_califa_hitenergy_rand_bg->Fill(Erel, Energy_dc);
                                    else
                                        fh_califa_hitenergy_rand_ag->Fill(Erel, Energy_dc);
                                }
                            }

                            /*   if(motherId < 2433){
                               califaHit_true=true;
                               new ((*fCalifaHitItems)[fNofCalifaHitItems++])
                              R3BCalifaClusterData(hitCalifa->GetCrystalList(), hitCalifa->GetEnergyList(),
                                                                                           hitCalifa->GetTimeList(),
                                                                                           Energy_dc,
                                                                                           hitCalifa->GetNf(),
                                                                                           hitCalifa->GetNs(),
                                                                                           hitCalifa->GetTheta(),
                                                                                           hitCalifa->GetPhi(),
                                                                                           tof+tCtofd,
                                                                                           hitCalifa->GetClusterType());
                              }*/
                        }
                    }
                }

                {
                    if (fCalItemsCalifa)
                    {
                        // CALIFA
                        auto detCalifa = fCalItemsCalifa;
                        Int_t nCalCalifa = detCalifa->GetEntriesFast();

                        for (Int_t ical = 0; ical < nCalCalifa; ical++)
                        {
                            R3BCalifaCrystalCalData* calCalifa = (R3BCalifaCrystalCalData*)detCalifa->At(ical);
                            if (!calCalifa)
                                continue;

                            new ((*fCalifaCalItems)[fNofCalifaCalItems++])
                                R3BCalifaCrystalCalData(calCalifa->GetCrystalId(),
                                                        calCalifa->GetEnergy(),
                                                        calCalifa->GetNf(),
                                                        calCalifa->GetNs(),
                                                        calCalifa->GetTime(),
                                                        0.);
                        }
                    }
                    if (fMappedItemsCalifa)
                    {
                        // CALIFA
                        auto detCalifa = fMappedItemsCalifa;
                        Int_t nMappedCalifa = detCalifa->GetEntriesFast();

                        for (Int_t imap = 0; imap < nMappedCalifa; imap++)
                        {
                            R3BCalifaMappedData* mapCalifa = (R3BCalifaMappedData*)detCalifa->At(imap);
                            if (!mapCalifa)
                                continue;

                            new ((*fCalifaMappedItems)[fNofCalifaMappedItems++])
                                R3BCalifaMappedData(mapCalifa->GetCrystalId(),
                                                    mapCalifa->GetEnergy(),
                                                    mapCalifa->GetNf(),
                                                    mapCalifa->GetNs(),
                                                    mapCalifa->GetFebexTime(),
                                                    mapCalifa->GetWRTS(),
                                                    mapCalifa->GetOverflow(),
                                                    mapCalifa->GetPileup(),
                                                    mapCalifa->GetDiscard(),
                                                    mapCalifa->GetTot());
                        }
                    }
                }

                R3BLOG(debug, "Finished with Califa data");

                is_tracked = true;
                nTracksFound += 1;
                fh_energy->Fill(Ec, Ea);

                Double_t mche = mHe * mC / mO;
                Double_t deltaTheta26 =
                    0.028 / 91.2 * cos(theta_26 / TMath::RadToDeg()) * cos(theta_26 / TMath::RadToDeg());
                Double_t deltaErel = 2. * mche * sqrt(Ec / mC * Ea / mHe) * sin(theta_26 / TMath::RadToDeg()) *
                                     deltaTheta26 * 1000.; // Ea/Ec are in GeV
                //    cout<<"deltaErel: "<<deltaErel<<", deltaTheta26: "<<deltaTheta26<<endl;
                // fh_dErel->Fill(deltaErel);

                Double_t Ebeam = mO * sqrt(1. + (17410.6 / mO) * (17410.6 / mO)) - mO; // in MeV
                Double_t target_exc = Ebeam - (Ec * 1000. - mC) - (Ea * 1000. - mHe) - 7.16;
                // cout<<"exc: "<<Ebeam<<", "<<(Ec*1000. - mC)<<"; "<<(Ea*1000. - mHe) <<", "<<target_exc<<endl;
                fh_minv->Fill(target_exc);

                //  fh_chiy_vs_chix->Fill(chiHex,chiCx);
                fh_target_xy->Fill(XHe, YHe);
                fh_target_xy->Fill(XC, YC);

                fh_mass->Fill(mtrackC, mtrackHe);

                fh_Erel_vs_event->Fill(fNEvents, Erel);

                fh_px_He->Fill(alpha.Px());
                fh_py_He->Fill(alpha.Py());
                fh_pz_He->Fill(alpha.Pz());
                fh_p_He->Fill(pa.Mag());

                fh_px_C->Fill(carbon.Px());
                fh_py_C->Fill(carbon.Py());
                fh_pz_C->Fill(carbon.Pz());
                fh_p_C->Fill(pc.Mag());

                if (Erel < 12. && Erel > 2.)
                {
                    fh_px_px->Fill(alpha.Px(), carbon.Px());
                    fh_py_py->Fill(alpha.Py(), carbon.Py());
                    fh_pz_pz->Fill(alpha.Pz(), carbon.Pz());
                    fh_p_p->Fill(pa.Mag(), pc.Mag());
                    fh_px_xfi23->Fill(posdetC.at(0), pCx);
                    fh_px_xfi23->Fill(posdetHe.at(0), pHex);
                    fh_py_yfi23->Fill(posdetC.at(3), pCx);
                    fh_py_yfi23->Fill(posdetHe.at(3), pHex);
                }

                fh_psum->Fill((pa + pc).Mag());
                fh_pzsum->Fill(alpha.Pz() + carbon.Pz());

                fh_theta26->Fill(theta_26);
                fh_Erel->Fill(Erel); // relative Energy plot

                fh_phi_bc_cm->Fill(phi_bc_cm);

                if (alpha.Px() < 0. && carbon.Px() < 0.)
                    fh_ErelR->Fill(Erel);
                if (alpha.Px() > 0. && carbon.Px() > 0.)
                    fh_ErelL->Fill(Erel);
                if ((alpha.Px() > 0. && carbon.Px() < 0.) || (alpha.Px() < 0. && carbon.Px() > 0.))
                {
                    // if (theta_26 > 0.7)
                    fh_ErelB->Fill(Erel);
                    fh_ErelB_vs_theta16O->Fill(theta_16O, Erel);
                }

                fh_Erel_vs_theta16O->Fill(theta_16O, Erel);

                fh_pHe_vs_theta26->Fill(pa.Mag(), theta_26 * TMath::DegToRad() * 1000.);
                fh_psum_vs_theta26->Fill(theta_26, (pa + pc).Mag());

                fh_Erel_vs_psum->Fill((pa + pc).Mag(), Erel);
                fh_theta_16O->Fill(theta_16O); // theta oxygen
                fh_phi_16O->Fill(phi_16O);

                fh_Erel_vs_x0->Fill(XC, Erel);
                fh_Erel_vs_y0->Fill(YC, Erel);

                fh_Erel_vs_px->Fill(alpha.Px(), Erel);
                fh_Erel_vs_px->Fill(carbon.Px(), Erel);
                fh_Erel_vs_py->Fill(alpha.Py(), Erel);
                fh_Erel_vs_py->Fill(carbon.Py(), Erel);

                fh_Erel_vs_xfi23->Fill(posdetC.at(0) - posdetHe.at(0), Erel);
                fh_Erel_vs_yfi23->Fill(posdetC.at(3) - posdetHe.at(3), Erel);

                fh_py_pz->Fill((pa + pc).Mag(), (pa + pc).Py());

                Double_t alpha_cm_Phi = alpha_cm.Phi() * TMath::RadToDeg();
                if (alpha_cm_Phi < 0.)
                    alpha_cm_Phi += 360;

                Double_t carbon_cm_Phi = carbon_cm.Phi() * TMath::RadToDeg();
                if (carbon_cm_Phi < 0.)
                    carbon_cm_Phi += 360;

                fh_theta_bc_cm->Fill(alpha_cm.Theta() * TMath::RadToDeg(), Erel);

                fh_phi_bc_cm_polar->Fill(phi_bc_cm, 1);

                if ((alpha.Px() > 0. && carbon.Px() < 0.) || (alpha.Px() < 0. && carbon.Px() > 0.))
                {
                    if (theta_16O < fThetaGrazing)
                    {
                        fh_ErelB_vs_phibc_bg->Fill(Erel, phi_bc_cm);
                        fh_ErelB_vs_theta26_bg->Fill(theta_26, Erel);
                        fh_erelB_vs_ptransHe_bg->Fill(ptransHe, Erel);
                        fh_erelB_vs_ptransC_bg->Fill(ptransC, Erel);
                        fh_thetaB_4He_cm_bg->Fill(alpha_cm.CosTheta(), Erel);
                        fh_phiB_4He_cm_bg->Fill(alpha_cm_Phi, Erel);
                        fh_thetaB_12C_cm_bg->Fill(carbon_cm.CosTheta(), Erel);
                        fh_phiB_12C_cm_bg->Fill(carbon_cm_Phi, Erel);
                    }
                    else
                    {
                        fh_ErelB_vs_phibc_ag->Fill(Erel, phi_bc_cm);
                        fh_ErelB_vs_theta26_ag->Fill(theta_26, Erel);
                        fh_erelB_vs_ptransHe_ag->Fill(ptransHe, Erel);
                        fh_erelB_vs_ptransC_ag->Fill(ptransC, Erel);
                        fh_thetaB_4He_cm_ag->Fill(alpha_cm.CosTheta(), Erel);
                        fh_phiB_4He_cm_ag->Fill(alpha_cm_Phi, Erel);
                        fh_thetaB_12C_cm_ag->Fill(carbon_cm.CosTheta(), Erel);
                        fh_phiB_12C_cm_ag->Fill(carbon_cm_Phi, Erel);
                    }
                }

                if (theta_16O < fThetaGrazing)
                {
                    fh_Erel_vs_phibc_bg->Fill(Erel, phi_bc_cm);
                    fh_Erel_vs_theta26_bg->Fill(theta_26, Erel);
                    fh_erel_vs_ptransHe_bg->Fill(ptransHe, Erel);
                    fh_erel_vs_ptransC_bg->Fill(ptransC, Erel);
                    fh_theta_4He_cm_bg->Fill(alpha_cm.CosTheta(), Erel);
                    fh_phi_4He_cm_bg->Fill(alpha_cm_Phi, Erel);
                    fh_theta_12C_cm_bg->Fill(carbon_cm.CosTheta(), Erel);
                    fh_phi_12C_cm_bg->Fill(carbon_cm_Phi, Erel);
                }
                else
                {
                    fh_Erel_vs_phibc_ag->Fill(Erel, phi_bc_cm);
                    fh_Erel_vs_theta26_ag->Fill(theta_26, Erel);
                    fh_erel_vs_ptransHe_ag->Fill(ptransHe, Erel);
                    fh_erel_vs_ptransC_ag->Fill(ptransC, Erel);
                    fh_theta_4He_cm_ag->Fill(alpha_cm.CosTheta(), Erel);
                    fh_phi_4He_cm_ag->Fill(alpha_cm_Phi, Erel);
                    fh_theta_12C_cm_ag->Fill(carbon_cm.CosTheta(), Erel);
                    fh_phi_12C_cm_ag->Fill(carbon_cm_Phi, Erel);
                }

                fh_thetax_px->Fill(pHex, thetaxHe);
                fh_thetax_px->Fill(pCx, thetaxC);

                fh_thetay_py->Fill(pHey, thetayHe);
                fh_thetay_py->Fill(pCy, thetayC);

                fh_thetax_py->Fill(pHey, thetaxHe);
                fh_thetax_py->Fill(pCy, thetaxC);

                fh_thetay_px->Fill(pHex, thetayHe);
                fh_thetay_px->Fill(pCx, thetayC);
            }

            /*
                        // rotation in oxygen direction
                        //TLorentzRotation transform;
                        //transform.RotateX(oxygen.Px()).RotateY(oxygen.Py()).RotateZ(0.);
                        TVector3 transform = oxygen.Vect(); /// TODO: check if it needs to be normed
                        alpha *= transform; // apply transformation vector
                        carbon *= transform;

                        if(abs(alpha.Phi()*TMath::RadToDeg()) > 88 && abs(alpha.Phi()*TMath::RadToDeg()) < 92)
                        {
                        fh_target_xy->Fill(XHe * 100., YHe * 100.);

                        fh_px_He->Fill(pHex);
                        fh_py_He->Fill(pHey);
                        fh_pz_He->Fill(pHez);

                        fh_px_C->Fill(pCx);
                        fh_py_C->Fill(pCy);
                        fh_pz_C->Fill(pCz);

                            fh_theta26->Fill(theta26); ///angle alpha carbon (theta)
                            fh_theta_16O->Fill(theta16O); // theta oxygen


                            fh_Erel->Fill(Erel); // relative Energy plot
                        fh_cs_cd->Fill(Erel,Erel);


            */

        } // end if chi2
    }     // end if trackHits>1

    LOG(debug) << "finished with track" << endl;

    if (fTrack && is_oxygen && !fPairs)
    {

        fh_chi2->Fill(chix);
        fh_psum_vs_chi->Fill(chix, Pf_tot);
        fh_target_xy_nc->Fill(Xf, Yf);
        fh_px_O_nc->Fill(Pxf);
        fh_py_O_nc->Fill(Pyf);
        fh_pz_O_nc->Fill(Pzf);
        fh_p_O_nc->Fill(Pf_tot);
        fh_psum_vs_event_nc->Fill(fNEvents, Pf_tot);
        Double_t Radius = sqrt(Xf * Xf + Yf * Yf);

        if (chix < fcut_chiX && chiy < fcut_chiY)
        {
            is_tracked = true;

            fh_target_xy->Fill(Xf, Yf);

            fh_px_O->Fill(Pxf);
            fh_py_O->Fill(Pyf);
            fh_pz_O->Fill(Pzf);
            fh_p_O->Fill(Pf_tot);
            // fh_psum_vs_event->Fill(fNEvents, Pf_tot);
            fh_py_pz->Fill(Pf_tot, Pyf);

        } // end if chi2
    }     // end if trackHits>1

    if (fSimu && nHitsMCTrack > 0)
    {
        if (fPairs) //&& is_tracked)
        {

            Double_t theta_26MC = alphaMC.Angle(carbonMC.Vect()) * TMath::RadToDeg(); // angle alpha carbon (theta)
            fh_theta26_simu->Fill(theta_26MC);
            // Double_t randtheta = (std::rand() / (float)RAND_MAX) - 0.5;
            //** Calculate invariant mass and relative energy **/
            m_invaMC = (alphaMC + carbonMC).M(); // invariant mass
            ErelMC = m_invaMC - mHe - mC;        // relative Energy
            fh_Erel_vs_thetaMC->Fill(theta_26MC, ErelMC);

            /* * transfer to cm system and rotate * */

            oxygenMC = alphaMC + carbonMC; // excited oxygen
            // cout<<"gamma "<<oxygen.Gamma()<<endl;
            Double_t thetaMC_16O = oxygenMC.Theta() * TMath::RadToDeg();
            Double_t phiMC_16O = oxygenMC.Phi() * TMath::RadToDeg();
            if (phiMC_16O < 0)
                phiMC_16O += 360.;

            // transfer to cm system and make some rotations
            // boost them to centre of mass
            alphaMC_cm = alphaMC;
            carbonMC_cm = carbonMC;
            oxygenMC_cm = oxygenMC;
            alphaMC_cm.Boost(-oxygenMC.BoostVector());
            carbonMC_cm.Boost(-oxygenMC.BoostVector());
            oxygenMC_cm.Boost(-oxygenMC.BoostVector());

            Double_t phiMC = (alphaMC_cm.Phi() - carbonMC_cm.Phi()) * TMath::RadToDeg();
            Double_t thetaMC = alphaMC_cm.Angle(carbonMC_cm.Vect()) * TMath::RadToDeg();

            fh_psum_MC->Fill((paMC + pcMC).Mag());
            fh_pzsum_MC->Fill(alphaMC.Pz() + carbonMC.Pz());
            fh_pxsum_MC->Fill(alphaMC.Px() + carbonMC.Px());
            fh_pysum_MC->Fill(alphaMC.Py() + carbonMC.Py());

            phiMC_bc_cm = (oxygenMC.Phi() - alphaMC_cm.Phi()) * TMath::RadToDeg();
            if (phiMC_bc_cm < 0)
                phiMC_bc_cm += 360.;
            fh_phiMC_bc_cm->Fill(phiMC_bc_cm);

            fh_phiMC_bc_cm_polar->Fill(phiMC_bc_cm, 1);

            fh_minv_simu->Fill(m_invaMC);
            fh_Erel_vs_phibcMC->Fill(ErelMC, phiMC_bc_cm);
            fh_Erel_simu->Fill(ErelMC);

            fh_px_px_mc->Fill(pHex_mc, pCx_mc);
            fh_py_py_mc->Fill(pHey_mc, pCy_mc);
            fh_pz_pz_mc->Fill(pHez_mc, pCz_mc);

            paMC_cm = alphaMC_cm.Vect();
            pcMC_cm = carbonMC_cm.Vect();

            poexc = oxygenMC.Vect();
            TVector3 norm;
            Double_t normx = (-1.) * pbeam.Z() * poexc.Y();
            Double_t normy = (-1.) * pbeam.Z() * poexc.X();
            norm.SetXYZ(normx, normy, 0); // vector perpendicular to (pbeam,p16O*)
            Double_t thetanorm_alphaMC = paMC_cm.Theta();
            Double_t thetanorm_carbonMC = pcMC_cm.Theta();
            Double_t pinplaneHeMC = paMC_cm.Perp(norm); // project to plane defined by pbeam and p16O*
            Double_t pinplaneCMC = pcMC_cm.Perp(norm);  // project to plane defined by pbeam and p16O*
            Double_t ptransHeMC = pinplaneHeMC * paMC_cm.Phi();
            Double_t ptransCMC = pinplaneCMC * pcMC_cm.Phi();
            Double_t alphaMC_cm_Phi = alphaMC_cm.Phi() * TMath::RadToDeg();
            if (alphaMC_cm_Phi < 0.)
                alphaMC_cm_Phi += 360;

            Double_t carbonMC_cm_Phi = carbonMC_cm.Phi() * TMath::RadToDeg();
            if (carbonMC_cm_Phi < 0.)
                carbonMC_cm_Phi += 360;

            if ((alphaMC.Px() > 0. && carbonMC.Px() < 0.) || (alphaMC.Px() < 0. && carbonMC.Px() > 0.))
            {
                fh_ErelB_vs_theta16O_MC->Fill(thetaMC_16O, ErelMC);
                fh_ErelB_vs_phibc_MC->Fill(ErelMC, phiMC_bc_cm);
                fh_erelB_vs_ptransHe_MC->Fill(ptransHeMC, ErelMC);
                fh_erelB_vs_ptransC_MC->Fill(ptransCMC, ErelMC);
                fh_thetaB_4He_cm_MC->Fill(alphaMC_cm.CosTheta(), ErelMC);
                fh_phiB_4He_cm_MC->Fill(alphaMC_cm_Phi, ErelMC);
                fh_thetaB_12C_cm_MC->Fill(carbonMC_cm.CosTheta(), ErelMC);
                fh_phiB_12C_cm_MC->Fill(carbonMC_cm_Phi, ErelMC);
            }

            if (is_tracked)
            {
                fh_dtheta->Fill(theta_26MC - theta_26);

                fh_dx->Fill(XHe_mc - XHe);
                fh_dx->Fill(XC_mc - XC);
                fh_dy->Fill(YHe_mc - YHe);
                fh_dy->Fill(YC_mc - YC);
                fh_target_xy_MC->Fill(XHe_mc, YHe_mc);
                fh_target_xy_MC->Fill(XC_mc, YC_mc);

                //   cout << "X-momenta: MC_He: " << pHex_mc << ", tracked He: " << pHex << endl;
                //   cout << "X-momenta: MC_C: " << pCx_mc << ", tracked C: " << pCx << endl;
                //   cout << "Y-momenta: MC_He: " << pHey_mc << ", tracked He: " << pHey << endl;
                //   cout << "Y-momenta: MC_C: " << pCy_mc << ", tracked C: " << pCy << endl;
                //   cout << "Momenta: MC_He: " << paMC.Mag() << ", tracked He: " << pa.Mag() << endl;
                //   cout << "Momenta: MC_C: " << pcMC.Mag() << ", tracked C: " << pc.Mag() << endl;

                fh_dpxC->Fill((pCx_mc - pCx) / pCx_mc * 100.);
                fh_dpxHe->Fill((pHex_mc - pHex) / pHex_mc * 100.);
                fh_dpyC->Fill((pCy_mc - pCy) / pCy_mc * 100.);
                fh_dpyHe->Fill((pHey_mc - pHey) / pHey_mc * 100.);
                fh_dpzC->Fill((pCz_mc - pCz) / pCz_mc * 100.);
                fh_dpzHe->Fill((pHez_mc - pHez) / pHez_mc * 100.);
                fh_dpC->Fill((pcMC.Mag() - pc.Mag()) / pcMC.Mag() * 100.);
                fh_dpHe->Fill((paMC.Mag() - pa.Mag()) / paMC.Mag() * 100.);
                fh_dp->Fill(((paMC + pcMC).Mag() - (pa + pc).Mag()) / ((paMC + pcMC).Mag()) * 100.);
                fh_dErel->Fill((ErelMC - Erel) / ErelMC * 100.);

                Double_t diff =
                    sqrt(pCx_mc * pCx_mc + pCy_mc * pCy_mc + pCz_mc * pCz_mc) - sqrt(pCx * pCx + pCy * pCy + pCz * pCz);
                //	cout<<"Difference: "<<pcMC.Mag() - pc.Mag()<<", "<<diff<<endl;

                fh_x_dpx->Fill((pCx_mc - pCx) / pCx_mc * 100., pCx);
                fh_x_dpx->Fill((pHex_mc - pHex) / pHex_mc * 100., pHex);
                fh_y_dpy->Fill((pCy_mc - pCy) / pCy_mc * 100., pCy);
                fh_y_dpy->Fill((pHey_mc - pHey) / pHey_mc * 100., pHey);
                fh_dpy_dpx->Fill((pCx_mc - pCx) / pCx_mc * 100., (pCy_mc - pCy) / pCy_mc * 100.);
                fh_dpy_dpx->Fill((pHex_mc - pHex) / pHex_mc * 100., (pHey_mc - pHey) / pHey_mc * 100.);

                // Comparison with simulation
                Double_t tthetaxC_simu = pCx_mc / pCz_mc;    // * 1000.;
                Double_t tthetayC_simu = pCy_mc / pCz_mc;    // * 1000.;
                Double_t tthetaxHe_simu = pHex_mc / pHez_mc; // * 1000.;
                Double_t tthetayHe_simu = pHey_mc / pHez_mc; // *  1000.;

                Double_t thetaxC_simu = atan(tthetaxC_simu) * 1000.;
                Double_t thetayC_simu = atan(tthetayC_simu) * 1000.;
                Double_t thetaxHe_simu = atan(tthetaxHe_simu) * 1000.;
                Double_t thetayHe_simu = atan(tthetayHe_simu) * 1000.;

                fh_thetax_dpx_He->Fill((pHex_mc - pHex), thetaxHe_simu);
                fh_thetax_dpx_C->Fill((pCx_mc - pCx), thetaxC_simu);

                fh_thetay_dpy_He->Fill((pHey_mc - pHey), thetayHe_simu);
                fh_thetay_dpy_C->Fill((pCy_mc - pCy), thetayC_simu);

                fh_thetax_dpy->Fill((pHey_mc - pHey), thetaxHe_simu);
                fh_thetax_dpy->Fill((pCy_mc - pCy), thetaxC_simu);

                fh_thetay_dpx->Fill((pHex_mc - pHex), thetayHe_simu);
                fh_thetay_dpx->Fill((pCx_mc - pCx), thetayC_simu);

                fh_px_xfi23_MC->Fill(posdetC.at(0), pCx_mc);
                fh_px_xfi23->Fill(posdetC.at(0), pCx);
                fh_px_xfi23_MC->Fill(posdetHe.at(0), pHex_mc);
                fh_px_xfi23->Fill(posdetHe.at(0), pHex);

                fh_py_yfi23_MC->Fill(posdetC.at(3), pCx_mc);
                fh_py_yfi23->Fill(posdetC.at(3), pCx);
                fh_py_yfi23_MC->Fill(posdetHe.at(3), pHex_mc);
                fh_py_yfi23->Fill(posdetHe.at(3), pHex);
            }
        }
        else
        {

            if (is_tracked)
            {

                fh_dx->Fill((Xf_mc - Xf));
                fh_dy->Fill((Yf_mc - Yf));

                fh_dpxO->Fill((Pxf_mc - Pxf) / Pxf_mc * 100.);
                fh_dpyO->Fill((Pyf_mc - Pyf) / Pyf_mc * 100.);

                fh_dpzO->Fill((Pzf_mc - Pzf) / Pzf_mc * 100.);
                fh_dpO->Fill((Pf_tot_mc - Pf_tot) / Pf_tot_mc * 100.);

                // fh_x_dpx->Fill((Pxf_mc - Pxf) / Pxf_mc * 100., Xf);
                fh_x_dpx->Fill(posdetO.at(0), Pxf_mc / Pxf);
                // fh_y_dpy->Fill((Pyf_mc - Pyf) / Pyf_mc * 100., Pyf);
                fh_y_dpy->Fill(posdetO.at(3), Pyf_mc / Pyf);

                fh_py_yfi23_MC->Fill(posdetO.at(3), Pyf_mc);
                fh_py_yfi23->Fill(posdetO.at(3), Pyf);

                fh_dpx_2d->Fill(posdetO.at(0), posdetO.at(3), Pxf_mc / Pxf);
                fh_dpy_2d->Fill(posdetO.at(0), posdetO.at(3), Pyf_mc / Pyf);
                fh_dpz_2d->Fill(posdetO.at(0), posdetO.at(3), Pzf_mc / Pzf);
                fh_dp_2d->Fill(posdetO.at(0), posdetO.at(3), Pf_tot_mc / Pf_tot);

                //         cout<<"ptot, x, y: "<<Pf_tot<<", "<<posdetO.at(0)<<"; "<<posdetO.at(3)<<endl;

                fh_px_xfi23_MC->Fill(posdetO.at(0), Pxf_mc);
                fh_px_xfi23->Fill(posdetO.at(0), Pxf);

                // fh_dpy_dpx->Fill((Pxf_mc - Pxf) , (Pyf_mc - Pyf) );
                fh_dpy_dpx->Fill(Pzf, Pyf);

                fh_psum_MC->Fill(Pf_tot_mc);
                fh_pzsum_MC->Fill(Pzf_mc);
                fh_pxsum_MC->Fill(Pxf_mc);
                fh_pysum_MC->Fill(Pyf_mc);
                fh_py_pz_MC->Fill(Pf_tot_mc, Pyf_mc);
                fh_px_px->Fill(Pxf_mc, Pxf);
                fh_py_py->Fill(Pyf_mc, Pyf);
            }
        }
    }

    fNEvents += 1;

    return;
}

void R3BGlobalAnalysisS494::FinishEvent()
{
    fTrack->Clear();
    if (fMCTrack)
        fMCTrack->Clear();
    if (fMappedItemsCalifa)
        fMappedItemsCalifa->Clear();
    if (fCalItemsCalifa)
        fCalItemsCalifa->Clear();
    if (fHitItemsCalifa)
        fHitItemsCalifa->Clear();
    fNofCalifaMappedItems = 0;
    fCalifaMappedItems->Clear();
    fNofCalifaCalItems = 0;
    fCalifaCalItems->Clear();
    fNofCalifaHitItems = 0;
    fCalifaHitItems->Clear();
}

void R3BGlobalAnalysisS494::FinishTask()
{
    if (fvis)
    {
        if (fPairs)
        {
            TCanvas* check1 = new TCanvas("SelectedChi2", "SelectedChi2", 10, 10, 1200, 900);
            check1->Divide(4, 3);
            check1->cd(1);
            gPad->SetLogz();
            fh_px_px->Draw("colz");
            check1->cd(2);
            gPad->SetLogz();
            fh_py_py->Draw("colz");
            check1->cd(3);
            gPad->SetLogz();
            fh_pz_pz->Draw("colz");
            check1->cd(4);
            fh_p_He->Draw();
            check1->cd(5);
            fh_p_C->Draw();
            check1->cd(6);
            fh_psum->Draw();
            check1->cd(7);
            fh_Erel->Draw();
            check1->cd(8);
            fh_theta26->Draw();
            check1->cd(9);
            fh_Erel_vs_theta26_bg->Draw("colz");
            check1->cd(10);
            fh_erel_vs_ptransHe_bg->Draw("colz");
            check1->cd(11);
            fh_Erel_vs_theta16O->Draw("colz");
            check1->cd(12);
            fh_Erel_vs_phibc_bg->Draw("colz");
        }
        else
        {
            TCanvas* check = new TCanvas("SelectedChi2", "SelectedChi2", 10, 10, 900, 900);
            check->Divide(3, 3);
            check->cd(1);
            fh_chi2->Draw();
            check->cd(2);
            fh_p_O_nc->Draw();
            check->cd(3);
            fh_px_O_nc->Draw();
            check->cd(4);
            fh_py_O_nc->Draw();
            check->cd(5);
            fh_pz_O_nc->Draw();
            check->cd(6);
            fh_p_O->Draw();
            check->cd(7);
            // fh_px_O->Draw();
            fh_py_pz->Draw("colz");
            check->cd(8);
            fh_py_O->Draw();
            check->cd(9);
            fh_pz_O->Draw();
        }

        if (fSimu)
        {
            TCanvas* checkMC = new TCanvas("CheckingGlobalMC", "CheckingGlobalMC", 10, 10, 1200, 900);
            TCanvas* checkReso = new TCanvas("SimuReso", "SimuReso", 10, 10, 1500, 1000);
            if (fPairs)
            {
                checkMC->Divide(4, 3);
                checkMC->cd(1);
                fh_dx->Draw();
                cout << "dx mean&sigma: " << fh_dx->GetMean() << "  " << fh_dx->GetStdDev() << endl;
                checkMC->cd(2);
                fh_dy->Draw();
                cout << "dy mean&sigma: " << fh_dy->GetMean() << "  " << fh_dy->GetStdDev() << endl;
                checkMC->cd(3);
                fh_dp->Draw();
                cout << "dp mean&sigma: " << fh_dp->GetMean() << "  " << fh_dp->GetStdDev() << endl;
                checkMC->cd(4);
                fh_dpzC->Draw();
                checkMC->cd(5);
                fh_dpzHe->Draw();
                checkMC->cd(6);
                fh_dErel->Draw();
                cout << "dE mean&sigma: " << fh_dErel->GetMean() << "  " << fh_dErel->GetStdDev() << endl;
                checkMC->cd(7);
                fh_x_dpx->Draw("colz");
                checkMC->cd(8);
                fh_y_dpy->Draw("colz");
                checkMC->cd(9);
                fh_target_xy->Draw("colz");
                checkMC->cd(10);
                gPad->SetLogz();
                fh_target_xy_MC->Draw("colz");
                checkMC->cd(11);
                gPad->SetLogz();
                fh_thetay_dpy_C->Draw("colz");
                checkMC->cd(12);
                gPad->SetLogz();
                fh_thetay_dpy_He->Draw("colz");

                checkReso->Divide(3, 2);
                checkReso->cd(1);
                fh_psum_MC->Draw();
                checkReso->cd(2);
                fh_psum->Draw();
                checkReso->cd(3);
                fh_dp->Draw();
                checkReso->cd(4);
                fh_Erel_simu->Draw();
                checkReso->cd(5);
                fh_Erel->Draw();
                checkReso->cd(6);
                fh_dErel->Draw();
            }
            else
            {
                checkMC->Divide(4, 3);
                checkMC->cd(1);
                fh_psum_MC->Draw();
                checkMC->cd(2);
                fh_pxsum_MC->Draw();
                checkMC->cd(3);
                fh_pysum_MC->Draw();
                checkMC->cd(4);
                fh_pzsum_MC->Draw();
                checkMC->cd(5);
                fh_dpO->Draw();
                checkMC->cd(6);
                fh_dpxO->Draw();
                checkMC->cd(7);
                fh_dpyO->Draw();
                checkMC->cd(8);
                fh_dpzO->Draw();
                checkMC->cd(9);
                fh_dp_2d->Draw("colz");
                checkMC->cd(10);
                fh_px_px->Draw("colz");
                checkMC->cd(11);
                fh_py_py->Draw("colz");
                checkMC->cd(12);
                fh_dpz_2d->Draw("colz");
            }
        }
    } // for fvis

    fh_target_xy->Write();
    fh_target_xx->Write();
    fh_target_yy->Write();
    fh_target_xy_nc->Write();

    fh_chi2->Write();
    fh_psum_vs_chi->Write();
    fh_Erel_vs_chi->Write();
    fh_phi26_vs_chi->Write();
    fh_psum_vs_theta26->Write();
    fh_chiy_vs_chix->Write();
    fh_chiy_vs_chix_nc->Write();
    fh_theta26_vs_chi->Write();
    fh_px_He_nc->Write();
    fh_py_He_nc->Write();
    fh_pz_He_nc->Write();
    fh_p_He_nc->Write();
    fh_px_C_nc->Write();
    fh_py_C_nc->Write();
    fh_pz_C_nc->Write();
    fh_p_C_nc->Write();
    fh_px_O_nc->Write();
    fh_py_O_nc->Write();
    fh_pz_O_nc->Write();
    fh_p_O_nc->Write();
    fh_px_px_nc->Write();
    fh_py_py_nc->Write();
    fh_pz_pz_nc->Write();
    fh_psum_vs_event_nc->Write();
    fh_psum_nc->Write();
    fh_p_vs_x0->Write();
    fh_psum_vs_x0->Write();
    fh_p_vs_y0->Write();
    fh_psum_vs_y0->Write();

    fh_theta26_nc->Write();
    fh_ErelB_nc->Write();
    fh_Erel_nc->Write();
    fh_theta_bc_cm_nc->Write();
    fh_phi_bc_cm_nc->Write();

    fh_px_He->Write();
    fh_py_He->Write();
    fh_pz_He->Write();
    fh_p_He->Write();
    fh_px_C->Write();
    fh_py_C->Write();
    fh_pz_C->Write();
    fh_p_C->Write();
    fh_px_O->Write();
    fh_py_O->Write();
    fh_pz_O->Write();
    fh_p_O->Write();
    fh_px_px->Write();
    fh_py_py->Write();
    fh_pz_pz->Write();
    fh_p_p->Write();
    fh_py_pz->Write();
    fh_psum->Write();
    fh_pzsum->Write();
    fh_Erel_vs_event->Write();
    fh_pHe_vs_theta26->Write();
    fh_psum_vs_theta26_nc->Write();
    fh_py_yfi23->Write();
    fh_px_xfi23->Write();

    fh_theta26->Write();
    fh_minv->Write();
    fh_Erel->Write();
    fh_ErelB->Write();
    fh_ErelR->Write();
    fh_ErelL->Write();
    fh_Erel_limit_tofd->Write();
    fh_Erel_limit_fi23->Write();
    fh_Erel_vs_psum->Write();
    fh_Erel_vs_px->Write();
    fh_Erel_vs_py->Write();
    fh_Erel_vs_x0->Write();
    fh_Erel_vs_y0->Write();
    fh_Erel_vs_xfi23->Write();
    fh_Erel_vs_yfi23->Write();
    fh_Erel_vs_dr_tofd->Write();
    fh_Erel_vs_dr_limit_tofd->Write();
    fh_Erel_vs_dr_fi23->Write();
    fh_Erel_vs_dr_limit_fi23->Write();
    fh_dr_vs_dr->Write();
    fh_dr_vs_dr_limit->Write();
    fh_Erel_vs_theta16O->Write();
    fh_ErelB_vs_theta16O->Write();
    fh_Erel_vs_theta16O_3He12C->Write();
    fh_phi_bc_cm->Write();

    fh_thetax_px->Write();
    fh_thetay_py->Write();
    fh_thetax_py->Write();
    fh_thetay_px->Write();

    fh_theta_16O->Write();
    fh_phi_16O->Write();
    fh_theta_bc_cm->Write();
    fh_phi_bc_cm_polar->Write();

    fh_ErelB_vs_phibc_bg->Write();
    fh_ErelB_vs_theta26_bg->Write();
    fh_erelB_vs_ptransHe_bg->Write();
    fh_erelB_vs_ptransC_bg->Write();
    fh_thetaB_4He_cm_bg->Write();
    fh_phiB_4He_cm_bg->Write();
    fh_thetaB_12C_cm_bg->Write();
    fh_phiB_12C_cm_bg->Write();
    fh_ErelB_vs_phibc_ag->Write();
    fh_ErelB_vs_theta26_ag->Write();
    fh_erelB_vs_ptransHe_ag->Write();
    fh_erelB_vs_ptransC_ag->Write();
    fh_thetaB_4He_cm_ag->Write();
    fh_phiB_4He_cm_ag->Write();
    fh_thetaB_12C_cm_ag->Write();
    fh_phiB_12C_cm_ag->Write();
    fh_Erel_vs_phibc_bg->Write();
    // fh_Erel_vs_theta26->Write();
    // fh_Erel_vs_theta26_max->Write();
    fh_Erel_vs_theta26_bg->Write();
    fh_erel_vs_ptransHe_bg->Write();
    fh_erel_vs_ptransC_bg->Write();
    fh_theta_4He_cm_bg->Write();
    fh_phi_4He_cm_bg->Write();
    fh_theta_12C_cm_bg->Write();
    fh_phi_12C_cm_bg->Write();
    fh_Erel_vs_phibc_ag->Write();
    fh_Erel_vs_theta26_ag->Write();
    fh_erel_vs_ptransHe_ag->Write();
    fh_erel_vs_ptransC_ag->Write();
    fh_theta_4He_cm_ag->Write();
    fh_phi_4He_cm_ag->Write();
    fh_theta_12C_cm_ag->Write();
    fh_phi_12C_cm_ag->Write();

    fh_mass_nc->Write();
    fh_mass->Write();
    fh_energy_nc->Write();
    fh_energy->Write();

    if (fSimu)
    {
        fh_target_xy_MC->Write();
        fh_dx->Write();
        fh_dy->Write();
        fh_dz->Write();
        fh_dpxO->Write();
        fh_dpyO->Write();
        fh_dpzO->Write();
        fh_dpO->Write();
        fh_dpxC->Write();
        fh_dpyC->Write();
        fh_dpzC->Write();
        fh_dpC->Write();
        fh_dpxHe->Write();
        fh_dpyHe->Write();
        fh_dpzHe->Write();
        fh_dpHe->Write();
        fh_dp->Write();
        fh_py_pz_MC->Write();
        fh_thetax_dpx_C->Write();
        fh_thetay_dpy_C->Write();
        fh_thetax_dpx_He->Write();
        fh_thetay_dpy_He->Write();
        fh_thetax_dpy->Write();
        fh_thetay_dpx->Write();
        fh_x_dpx->Write();
        fh_y_dpy->Write();
        fh_dpy_dpx->Write();
        fh_Erel_vs_thetaMC->Write();
        fh_dtheta->Write();
        fh_psum_MC->Write();
        fh_pzsum_MC->Write();
        fh_pxsum_MC->Write();
        fh_pysum_MC->Write();
        fh_theta26_simu->Write();
        fh_Erel_simu->Write();
        fh_phiMC_bc_cm->Write();
        fh_phiMC_bc_cm_polar->Write();
        fh_Erel_vs_phibcMC->Write();
        fh_px_px_mc->Write();
        fh_py_py_mc->Write();
        fh_pz_pz_mc->Write();
        fh_minv_simu->Write();
        fh_ErelB_vs_phibc_MC->Write();
        fh_erelB_vs_ptransHe_MC->Write();
        fh_erelB_vs_ptransC_MC->Write();
        fh_thetaB_4He_cm_MC->Write();
        fh_phiB_4He_cm_MC->Write();
        fh_thetaB_12C_cm_MC->Write();
        fh_phiB_12C_cm_MC->Write();
        fh_ErelB_vs_theta16O_MC->Write();
        fh_py_yfi23_MC->Write();
        fh_px_xfi23_MC->Write();
        fh_dp_2d->Write();
        fh_dpx_2d->Write();
        fh_dpy_2d->Write();
        fh_dpz_2d->Write();
        fh_psum_vs_dx0->Write();
        fh_pC_vs_dx0->Write();
        fh_pC_vs_dx0_mc->Write();
        fh_x0_vs_dx0->Write();
        fh_pC_vs_x0->Write();
        fh_pC_vs_x0_mc->Write();
        fh_dpC_vs_dx0->Write();
        fh_pHe_vs_dx0->Write();
        fh_pHe_vs_dx0_mc->Write();
        fh_pHe_vs_x0->Write();
        fh_pHe_vs_x0_mc->Write();
        fh_dpHe_vs_dx0->Write();
    }
    fh_dErel->Write();

    if (fHitItemsCalifa)
    {
        fh_Erel_withCalifa->Write();
        fh_califa_hitenergy_boost->Write();
        fh_califa_barrel_hitenergy_boost->Write();
        fh_califa_iphos_hitenergy_boost->Write();
        fh_califa_hitenergy_ag->Write();
        fh_califa_hitenergy_bg->Write();
        fh_califa_hitenergy_rand_ag->Write();
        fh_califa_hitenergy_rand_bg->Write();
        fh_Erel_withCalifa_2d->Write();
        fh_califaenergy_2d->Write();
        fh_califaenergy_barrel_2d->Write();
        fh_califaenergy_iphos_2d->Write();
        fh_Erel_vs_theta16O_withcalifa->Write();

        fh_califa_hitenergy_select->Write();
        fh_Erel_withCalifa_tof->Write();
        fh_Erel_withCalifa_motherId->Write();
        fh_phi_vs_ecalifa_barrel->Write();
        fh_phi_vs_ecalifa_iphos->Write();
        fh_califa_tofd->Write();
        fh_califa_tofd_cut->Write();
        fh_califa_tofd_ecr->Write();
        fh_crystalNb->Write();
        fh_califa_ecl->Write();
        fh_califa_dtime_ecr->Write();
        fh_califa_dtof_ecr->Write();
        fh_califa_time_cId->Write();
        fh_califa_hitenergycorr_boost->Write();
        fh_esumcluster_dc_2d->Write();
        fh_califa_hitenergy_boost_rand->Write();
        fh_Erel_vs_theta16O_withcalifa_rand->Write();
    }
    for (int det = 0; det < N_DET_MAX; det++)
    {
        fh_psum_vs_xdet[det]->Write();
        fh_psum_vs_ydet[det]->Write();
        fh_p_vs_xdet[det]->Write();
    }
    fh_yfi23_vs_ytofd_bc->Write();
    fh_yfi23_vs_ytofd->Write();
    fh_xfi23_vs_xtofd_bc->Write();
    fh_xfi23_vs_xtofd->Write();
    fh_xy_fib23->Write();
    fh_xx_fib23->Write();
    fh_yy_fib23->Write();
    fh_xy_fib23_bc->Write();
    fh_dt_fib23_bc->Write();

    /*
        ofstream myfile("/u/kelic/R3BRoot/macros/r3b/unpack/s494/simulations/generators/checkpx.dat", std::ios::out);

        Int_t nchy = fh_px_xfi23->GetNbinsY();
        Int_t nchx = fh_px_xfi23->GetNbinsX();
        Int_t nch = nchx;
        Double_t px[nch],pxmc[nch],x[nch],xmc[nch];
        TAxis* yAxis = fh_px_xfi23->GetYaxis();
        TAxis* yAxis_MC = fh_px_xfi23_MC->GetYaxis();
        TAxis* xAxis = fh_px_xfi23->GetXaxis();
        TAxis* xAxis_MC = fh_px_xfi23_MC->GetXaxis();

        for(Int_t i=1;i<nch+1;i++)
        {
            TH1D* proj1 = fh_px_xfi23->ProjectionY("", i , i , 0);
            px[i] = proj1->GetMean();
            x[i] = xAxis->GetBinCenter(i);

            TH1D* proj2 = fh_px_xfi23_MC->ProjectionY("", i , i , 0);
            pxmc[i] = proj2->GetMean();
            xmc[i] = xAxis_MC->GetBinCenter(i);

            myfile<<px[i]<<" "<<x[i]<<" " <<pxmc[i]<<" "<<xmc[i]<<" "<<px[i]-pxmc[i]<<endl;
        }

        ofstream myfile1("/u/kelic/R3BRoot/macros/r3b/unpack/s494/simulations/generators/checkpy.dat", std::ios::out);

        nchy = fh_py_yfi23->GetNbinsY();
        nchx = fh_py_yfi23->GetNbinsX();
        nch = nchx;
        Double_t py[nch],pymc[nch],y[nch],ymc[nch];
        yAxis = fh_py_yfi23->GetYaxis();
        yAxis_MC = fh_py_yfi23_MC->GetYaxis();
        xAxis = fh_py_yfi23->GetXaxis();
        xAxis_MC = fh_py_yfi23_MC->GetXaxis();

        for(Int_t i=1;i<nch+1;i++)
        {
            TH1D* proj1 = fh_py_yfi23->ProjectionY("", i , i , 0);
            py[i] = proj1->GetMean();
            y[i] = xAxis->GetBinCenter(i);

            TH1D* proj2 = fh_py_yfi23_MC->ProjectionY("", i , i , 0);
            pymc[i] = proj2->GetMean();
            ymc[i] = xAxis_MC->GetBinCenter(i);

            myfile1<<py[i]<<" "<<y[i]<<" " <<pymc[i]<<" "<<ymc[i]<<" "<<py[i]-pymc[i]<<endl;
        }

    */
}
ClassImp(R3BGlobalAnalysisS494)
