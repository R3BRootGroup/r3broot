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
// -----                  R3BGlobalAnalysisS454                -----
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

#include "R3BBeamMonitorMappedData.h"

#include "R3BGlobalAnalysisS454.h"

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

R3BGlobalAnalysisS454::R3BGlobalAnalysisS454()
    : R3BGlobalAnalysisS454("GlobalAnalysis", 1)
{
}

R3BGlobalAnalysisS454::R3BGlobalAnalysisS454(const char* name, Int_t iVerbose)
    : FairTask(name, iVerbose)
    , fTrigger(-1)
    , fTpat(-1)
    , fCuts(0)
    , fGhost(0)
    , fPairs(0)
    , fB(-1672)
    , fNEvents(0)
{
}

R3BGlobalAnalysisS454::~R3BGlobalAnalysisS454() {}

InitStatus R3BGlobalAnalysisS454::Init()
{
    // Initialize random number:
    std::srand(std::time(0)); // use current time as seed for random generator

    LOG(info) << "R3BGlobalAnalysisS454::Init ";

    // try to get a handle on the EventHeader. EventHeader may not be
    // present though and hence may be null. Take care when using.

    FairRootManager* mgr = FairRootManager::Instance();
    if (NULL == mgr)
        LOG(fatal) << "FairRootManager not found";

    header = dynamic_cast<R3BEventHeader*>(mgr->GetObject("R3BEventHeader"));
    FairRunOnline* run = FairRunOnline::Instance();

    // Get objects for detectors on all levels
    fMCTrack = dynamic_cast<TClonesArray*>(mgr->GetObject("MCTrack"));
    fTrack = dynamic_cast<TClonesArray*>(mgr->GetObject("Track"));
    maxevent = mgr->CheckMaxEventNo();

    //------------------------------------------------------------------------
    // create histograms of all detectors
    //------------------------------------------------------------------------

    char strNameC[255];
    sprintf(strNameC, "Cave_C_position");
    fh_Cave_position = new TH2F(strNameC, "", 2100, -100., 2000., 1000, -500., 500.);

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

    fh_target_xy = new TH2F("target_xy", "target xy ", 1000, -0.5, 0.5, 1000, -0.5, 0.5);
    fh_target_xy->GetXaxis()->SetTitle("x / cm");
    fh_target_xy->GetYaxis()->SetTitle("y / cm");

    fh_chi2 = new TH1F("chi2", "chi2 ", 1000, 0., 100);
    fh_chi2->GetXaxis()->SetTitle("Chi2");
    fh_chi2->GetYaxis()->SetTitle("counts");

    fh_dx = new TH1F("tracker_dx", "tracker dx ", 2000, -5., 5);
    fh_dx->GetXaxis()->SetTitle("dx / cm");
    fh_dx->GetYaxis()->SetTitle("counts");

    fh_dy = new TH1F("tracker_dy", "tracker dy ", 2000, -5., 5.);
    fh_dy->GetXaxis()->SetTitle("dy / cm");
    fh_dy->GetYaxis()->SetTitle("counts");

    fh_dz = new TH1F("tracker_dz", "tracker dz ", 2000, -5., 5.);
    fh_dz->GetXaxis()->SetTitle("dz / cm");
    fh_dz->GetYaxis()->SetTitle("counts");

    fh_dpx = new TH1F("tracker_dpx", "tracker dpx ", 2000, -10, 10);
    fh_dpx->GetXaxis()->SetTitle("dpx / percent");
    fh_dpx->GetYaxis()->SetTitle("counts");

    fh_dpy = new TH1F("tracker_dpy", "tracker dpy ", 2000, -10, 10);
    fh_dpy->GetXaxis()->SetTitle("dpy / percent");
    fh_dpy->GetYaxis()->SetTitle("counts");

    fh_dpz = new TH1F("tracker_dpz", "tracker dpz ", 2000, -10, 10);
    fh_dpz->GetXaxis()->SetTitle("dpz / percent");
    fh_dpz->GetYaxis()->SetTitle("counts");

    fh_dp = new TH1F("tracker_dp", "tracker dp ", 2000, -10, 10);
    fh_dp->GetXaxis()->SetTitle("dp / percent");
    fh_dp->GetYaxis()->SetTitle("counts");

    fh_thetax_dpx = new TH2F("tracker_thethax_dpx", "tracker theta_x vs dpx ", 1000, -50, 50, 200, -100., 100.);
    fh_thetax_dpx->GetXaxis()->SetTitle("dpx / percent");
    fh_thetax_dpx->GetYaxis()->SetTitle("theta_x / mrad");

    fh_dpy_dpx = new TH2F("tracker_dpy_dpx", "tracker dpy vs dpx ", 200, -100, 100, 200, -100., 100.);
    fh_dpy_dpx->GetXaxis()->SetTitle("dpx / percent");
    fh_dpy_dpx->GetYaxis()->SetTitle("dpy / percent");

    fh_thetay_dpy = new TH2F("tracker_thethay_dpy", "tracker theta_y vs dpy ", 1000, -50, 50, 200, -100., 100.);
    fh_thetay_dpy->GetXaxis()->SetTitle("dpy / percent");
    fh_thetay_dpy->GetYaxis()->SetTitle("theta_y / mrad");

    fh_thetax_dpy = new TH2F("tracker_thethax_dpy", "tracker theta_x vs dpy ", 1000, -50, 50, 200, -100., 100.);
    fh_thetax_dpy->GetXaxis()->SetTitle("dpy / percent");
    fh_thetax_dpy->GetYaxis()->SetTitle("theta_x / mrad");

    fh_thetay_dpx = new TH2F("tracker_thethay_dpx", "tracker theta_y vs dpx ", 1000, -50, 50, 200, -100., 100.);
    fh_thetay_dpx->GetXaxis()->SetTitle("dpx / percent");
    fh_thetay_dpx->GetYaxis()->SetTitle("theta_y / mrad");

    fh_thetax_dpx_abs =
        new TH2F("tracker_thethax_dpx_abs", "tracker theta_x vs dpx abs ", 1000, -50, 50, 200, -100., 100.);
    fh_thetax_dpx_abs->GetXaxis()->SetTitle("dpx");
    fh_thetax_dpx_abs->GetYaxis()->SetTitle("theta_x / mrad");

    fh_x_dpx = new TH2F("tracker_x_dpx", "tracker x vs dpx ", 1000, -50, 50, 200, -10., 10.);
    fh_x_dpx->GetXaxis()->SetTitle("dpx / percent");
    fh_x_dpx->GetYaxis()->SetTitle("x / cm");

    fh_y_dpy = new TH2F("tracker_y_dpy", "tracker y vs dpy ", 1000, -50, 50, 200, -10., 10.);
    fh_y_dpy->GetXaxis()->SetTitle("dpy / percent");
    fh_y_dpy->GetYaxis()->SetTitle("y / cm");

    fh_px_px = new TH2F("px_px", "px He vs. px C", 1000, -1000, 1000, 1000, -1000., 1000.);
    fh_px_px->GetXaxis()->SetTitle("px / MeV/c");
    fh_px_px->GetYaxis()->SetTitle("px / MeV/c");

    fh_py_py = new TH2F("py_py", "py He vs. py C", 1000, -1000, 1000, 1000, -1000., 1000.);
    fh_py_py->GetXaxis()->SetTitle("py / MeV/c");
    fh_py_py->GetYaxis()->SetTitle("py / MeV/c");

    fh_pz_pz = new TH2F("pz_pz", "pz He vs. pz C", 1000, 0, 10000, 2000, 0., 20000.);
    fh_pz_pz->GetXaxis()->SetTitle("pz / MeV/c");
    fh_pz_pz->GetYaxis()->SetTitle("pz / MeV/c");

    fh_px_px_cm = new TH2F("px_px_cm", "px He vs. px C", 1000, -1000, 1000, 1000, -1000., 1000.);
    fh_px_px_cm->GetXaxis()->SetTitle("px / MeV/c");
    fh_px_px_cm->GetYaxis()->SetTitle("px / MeV/c");

    fh_py_py_cm = new TH2F("py_py_cm", "py He vs. py C", 1000, -1000, 1000, 1000, -1000., 1000.);
    fh_py_py_cm->GetXaxis()->SetTitle("py / MeV/c");
    fh_py_py_cm->GetYaxis()->SetTitle("py / MeV/c");

    fh_pz_pz_cm = new TH2F("pz_pz_cm", "pz He vs. pz C", 1000, -1000, 1000, 1000, -1000., 1000.);
    fh_pz_pz_cm->GetXaxis()->SetTitle("pz / MeV/c");
    fh_pz_pz_cm->GetYaxis()->SetTitle("pz / MeV/c");

    // Analysis
    fh_chiy_vs_chix_He = new TH2F("chiy_vs_chix_He", "chi y vs. chi x for He", 500, 0., 500., 500, 0., 500.);
    fh_chiy_vs_chix_He->GetXaxis()->SetTitle("chi x");
    fh_chiy_vs_chix_He->GetYaxis()->SetTitle("chi y");

    fh_chiy_vs_chix_C = new TH2F("chiy_vs_chix_C", "chi y vs. chi x for C", 500, 0., 500., 500, 0., 500.);
    fh_chiy_vs_chix_C->GetXaxis()->SetTitle("chi x");
    fh_chiy_vs_chix_C->GetYaxis()->SetTitle("chi y");

    fh_theta26 = new TH1F("theta26", "theta between alpha and carbon ", 500, 0., 5);
    fh_theta26->GetXaxis()->SetTitle("angle / degree");
    fh_theta26->GetYaxis()->SetTitle("counts");

    fh_theta_16O = new TH1F("theta_16O", "theta of excited 16O ", 500, 0., 5);
    fh_theta_16O->GetXaxis()->SetTitle("angle / degree");
    fh_theta_16O->GetYaxis()->SetTitle("counts");

    fh_phi_16O = new TH1F("phi_16O", "phi of excited 16O ", 360, 0., 360);
    fh_phi_16O->GetXaxis()->SetTitle("angle / degree");
    fh_phi_16O->GetYaxis()->SetTitle("counts");

    fh_theta26_cm = new TH1F("theta26_cm", "theta between alpha and carbon in cm-system ", 360, 0, 360.);
    fh_theta26_cm->GetXaxis()->SetTitle("angle / degree");
    fh_theta26_cm->GetYaxis()->SetTitle("counts");

    fh_phi26_cm = new TH1F("phi26_cm", "phi between alpha and carbon in cm-system", 360, 0., 360.);
    fh_phi26_cm->GetXaxis()->SetTitle("angle / degree");
    fh_phi26_cm->GetYaxis()->SetTitle("counts");

    fh_theta_4He_cm = new TH1F("theta_4He_cm", "Theta 4He in cm-system", 180, 0., 180.);
    fh_theta_4He_cm->GetXaxis()->SetTitle("angle / degree");
    fh_theta_4He_cm->GetYaxis()->SetTitle("counts");

    fh_phi_4He_cm = new TH1F("phi_4He_cm", "phi 4He in cm-system", 360, 0., 360.);
    fh_phi_4He_cm->GetXaxis()->SetTitle("angle / degree");
    fh_phi_4He_cm->GetYaxis()->SetTitle("counts");

    fh_theta_bc_cm = new TH1F("theta_bc_cm", "Theta bc in cm-system", 180, 0., 180.);
    fh_theta_bc_cm->GetXaxis()->SetTitle("angle / degree");
    fh_theta_bc_cm->GetYaxis()->SetTitle("counts");

    fh_phi_bc_cm = new TH1F("phi_bc_cm", "phi bc in cm-system", 360, 0., 360.);
    fh_phi_bc_cm->GetXaxis()->SetTitle("angle / degree");
    fh_phi_bc_cm->GetYaxis()->SetTitle("counts");

    fh_phi_bc_cm_polar = new TH2F("phi_bc_cm_polar", "phi_bc_cm_polar", 360, 0., 360., 100, 0., 1000.);
    fh_phi_bc_cm_polar->GetXaxis()->SetTitle("angle / degree");
    fh_phi_bc_cm_polar->GetYaxis()->SetTitle("counts");

    fh_theta_12C_cm = new TH1F("theta_12C_cm", "Theta 12C in cm-system", 180, 0., 180.);
    fh_theta_12C_cm->GetXaxis()->SetTitle("angle / degree");
    fh_theta_12C_cm->GetYaxis()->SetTitle("counts");

    fh_phi_12C_cm = new TH1F("phi_12C_cm", "phi 12C in cm-system", 360, 0., 360.);
    fh_phi_12C_cm->GetXaxis()->SetTitle("angle / degree");
    fh_phi_12C_cm->GetYaxis()->SetTitle("counts");

    fh_theta26_simu = new TH1F("theta26_simu", "theta 26 simulation", 500, 0., 5);
    fh_theta26_simu->GetXaxis()->SetTitle("angle / degree");
    fh_theta26_simu->GetYaxis()->SetTitle("counts");

    fh_Erel = new TH1F("Erel", "Erel ", 5000, 0., 100);
    fh_Erel->GetXaxis()->SetTitle("Erel / MeV");
    fh_Erel->GetYaxis()->SetTitle("counts");

    fh_ErelL = new TH1F("ErelL", "Erel links", 5000, 0., 100);
    fh_ErelL->GetXaxis()->SetTitle("Erel / MeV");
    fh_ErelL->GetYaxis()->SetTitle("counts");

    fh_ErelR = new TH1F("ErelR", "Erel rechts ", 5000, 0., 100);
    fh_ErelR->GetXaxis()->SetTitle("Erel / MeV");
    fh_ErelR->GetYaxis()->SetTitle("counts");

    fh_ErelB = new TH1F("ErelB", "Erel both sides ", 5000, 0., 100);
    fh_ErelB->GetXaxis()->SetTitle("Erel / MeV");
    fh_ErelB->GetYaxis()->SetTitle("counts");

    fh_Erel_simu = new TH1F("Erel_simu", "Erel simulation", 600, -10., 50);
    fh_Erel_simu->GetXaxis()->SetTitle("Erel / MeV");
    fh_Erel_simu->GetYaxis()->SetTitle("counts");

    fh_dErel_vs_x = new TH2F("dErel_vs_x", "delta Erel vs. x", 200, -100., 100., 100, -5., 5.);
    fh_dErel_vs_x->GetXaxis()->SetTitle("TofD x / cm");
    fh_dErel_vs_x->GetYaxis()->SetTitle("Delta Erel / MeV");

    fh_dErel_vs_y = new TH2F("dErel_vs_y", "delta Erel vs. y", 200, -100., 100., 100, -5., 5.);
    fh_dErel_vs_y->GetXaxis()->SetTitle("TofD y / cm");
    fh_dErel_vs_y->GetYaxis()->SetTitle("Delta Erel / MeV");

    return kSUCCESS;
}

void R3BGlobalAnalysisS454::Exec(Option_t* option)
{
    if (fNEvents / 100000. == (int)fNEvents / 100000)
        std::cout << "\rEvents: " << fNEvents << " / " << maxevent << " (" << (int)(fNEvents * 100. / maxevent)
                  << " %) " << std::flush;

    //    cout << "New event " << fNEvents << endl;

    FairRootManager* mgr = FairRootManager::Instance();
    if (NULL == mgr)
        LOG(fatal) << "FairRootManager not found";

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
    Double_t px, py, pz;
    Double_t theta_16, theta_26;
    Double_t costh26;
    Double_t chiHex = 0., chiHey = 0., chiCx = 0., chiCy = 0.;

    Bool_t is_alpha = false;
    Bool_t is_carbon = false;

    // Choose original MC data or exp. data
    Bool_t simData = true;

    if (fTrack && !simData)
    {

        // read in Track parameter
        Int_t nHitsTrack = fTrack->GetEntriesFast();
        // if(nHitsTrack>0) cout << "Track hits: " << nHitsTrack << endl;
        for (Int_t l = 0; l < nHitsTrack; l++)
        {
            R3BTrack* aTrack = dynamic_cast<R3BTrack*>(fTrack->At(l));

            LOG(debug) << "Charge " << aTrack->GetQ() << endl;

            if (aTrack->GetQ() == 2)
            {
                // 4He
                XHe = aTrack->GetX();
                YHe = aTrack->GetY();
                ZHe = aTrack->GetZ();

                pHex = aTrack->GetPx();
                pHey = aTrack->GetPy();
                pHez = aTrack->GetPz();

                chiHex = aTrack->GetChix();
                chiHey = aTrack->GetChiy();

                is_alpha = true;
                alpha.SetPxPyPzE(pHex, pHey, pHez, sqrt(pow(pHex, 2) + pow(pHey, 2) + pow(pHez, 2) + pow(mHe, 2)));

                // cout << "******************************************" << endl;
                // cout << "Track In 4He"
                //           << "x " << XHe << " y " << YHe << " z " << ZHe << endl;
                // cout << "px " << pHex << " py " << pHey << " z " << pHez << endl;
            }
            if (aTrack->GetQ() == 6)
            {
                // 12C
                XC = aTrack->GetX();
                YC = aTrack->GetY();
                ZC = aTrack->GetZ();

                pCx = aTrack->GetPx();
                pCy = aTrack->GetPy();
                pCz = aTrack->GetPz();

                chiCx = aTrack->GetChix();
                chiCy = aTrack->GetChiy();

                is_carbon = true;
                carbon.SetPxPyPzE(pCx, pCy, pCz, sqrt(pow(pCx, 2) + pow(pCy, 2) + pow(pCz, 2) + pow(mC, 2)));

                // cout << "******************************************" << endl;
                // cout << "Track In 12C"
                //           << "x " << XC << " y " << YC << " z " << ZC << endl;
                // cout << "px " << pCx << " py " << pCy << " z " << pCz << endl;
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

                // cout << "******************************************" << endl;
                // cout << "Track In 16O"
                //           << "x " << Xf << " y " << Yf << " z " << Zf << endl;
                // cout << "px " << Pxf << " py " << Pyf << " z " << Pzf << endl;
            }
        }
    }
    if (fMCTrack && simData)
    {
        // read in Monte Carlo Track parameter

        Int_t nHitsMCTrack = fMCTrack->GetEntriesFast();
        // cout << "MCTrack hits: " << nHitsMCTrack << endl;

        for (Int_t l = 0; l < nHitsMCTrack; l++)
        {
            // cout << "Original MC Data is analyzed" << endl;
            R3BMCTrack* aTrack = dynamic_cast<R3BMCTrack*>(fMCTrack->At(l));

            Int_t PID = aTrack->GetPdgCode();
            Int_t mother = aTrack->GetMotherId();
            LOG(debug) << "PID " << PID << endl;
            if (mother < 0)
            {
                if (PID == 1000020040)
                {
                    // 4He
                    XHe = aTrack->GetStartX();
                    YHe = aTrack->GetStartY();
                    ZHe = aTrack->GetStartZ();
                    THe = aTrack->GetStartT();

                    pHex = aTrack->GetPx() * 1000.;
                    pHey = aTrack->GetPy() * 1000.;
                    pHez = aTrack->GetPz() * 1000.;

                    is_alpha = true;
                    alpha.SetPxPyPzE(pHex, pHey, pHez, sqrt(pow(pHex, 2) + pow(pHey, 2) + pow(pHez, 2) + pow(mHe, 2)));

                    LOG(debug) << "******************************************" << endl;
                    LOG(debug) << "Track In 4He"
                               << "x " << XHe << " y " << YHe << " z " << ZHe << endl;
                    LOG(debug) << "px " << pHex << " py " << pHey << " pz " << pHez << endl;
                }
                if (PID == 1000060120)
                {
                    // 12C
                    XC = aTrack->GetStartX();
                    YC = aTrack->GetStartY();
                    ZC = aTrack->GetStartZ();
                    TC = aTrack->GetStartT();

                    pCx = aTrack->GetPx() * 1000.;
                    pCy = aTrack->GetPy() * 1000.;
                    pCz = aTrack->GetPz() * 1000.;

                    is_carbon = true;
                    carbon.SetPxPyPzE(pCx, pCy, pCz, sqrt(pow(pCx, 2) + pow(pCy, 2) + pow(pCz, 2) + pow(mC, 2)));

                    LOG(debug) << "******************************************" << endl;
                    LOG(debug) << "Track In 12C"
                               << "x " << XC << " y " << YC << " z " << ZC << endl;
                    LOG(debug) << "px " << pCx << " py " << pCy << " pz " << pCz << endl;
                }
                if (PID == 1000080160)
                {
                    // 16O
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

    // new version of analysis
    if (is_carbon && is_alpha)
    {
        fh_chiy_vs_chix_He->Fill(chiHex, chiHey);
        fh_chiy_vs_chix_C->Fill(chiCx, chiCy);
        Double_t cut_chi = 100;
        if (chiHex < cut_chi && chiHey < cut_chi && chiCx < cut_chi && chiCy < cut_chi)
        //			&&((alpha.Px() > 0 && carbon.Px() < 0) || (alpha.Px() < 0 && carbon.Px() > 0)))
        {
            fh_target_xy->Fill(XHe * 100., YHe * 100.);

            fh_px_He->Fill(alpha.Px());
            fh_py_He->Fill(alpha.Py());
            fh_pz_He->Fill(alpha.Pz());

            fh_px_C->Fill(carbon.Px());
            fh_py_C->Fill(carbon.Py());
            fh_pz_C->Fill(carbon.Pz());

            fh_px_px->Fill(alpha.Px(), carbon.Px());
            fh_py_py->Fill(alpha.Py(), carbon.Py());
            fh_pz_pz->Fill(alpha.Pz(), carbon.Pz());

            TVector3 pa = alpha.Vect();
            fh_p_He->Fill(pa.Mag());
            TVector3 pc = carbon.Vect();
            fh_p_C->Fill(pc.Mag());

            // Calculate angle between alpha and C
            if (alpha.Pz() == 0 || carbon.Pz() == 0)
                return;

            theta_26 = alpha.Angle(carbon.Vect()) * TMath::RadToDeg(); /// angle alpha carbon (theta)
            fh_theta26->Fill(theta_26);

            /** Calculate invariant mass and relative energy **/
            Double_t m_inva = (alpha + carbon).M(); // invariant mass
            Double_t Erel = m_inva - mHe - mC;      // relative Energy
            fh_Erel->Fill(Erel);                    // relative Energy plot

            if (alpha.Px() < 0. && carbon.Px() < 0.)
                fh_ErelR->Fill(Erel);
            if (alpha.Px() > 0. && carbon.Px() > 0.)
                fh_ErelL->Fill(Erel);
            if ((alpha.Px() > 0. && carbon.Px() < 0.) || (alpha.Px() < 0. && carbon.Px() > 0.))
                fh_ErelB->Fill(Erel);

            /** transfer to cm system and rotate **/

            oxygen = alpha + carbon; // excited oxygen
            // cout<<"gamma "<<oxygen.Gamma()<<endl;
            Double_t theta_16O = oxygen.Theta() * TMath::RadToDeg();
            fh_theta_16O->Fill(theta_16O); // theta oxygen
            Double_t phi_16O = oxygen.Phi() * TMath::RadToDeg();
            if (phi_16O < 0)
                phi_16O += 360.;
            fh_phi_16O->Fill(phi_16O);

            // photon.SetPxPyPzE(oxygen.Px, oxygen.Py, 0., sqrt(pow(oxygen.Px, 2) + pow(oxygen.Py, 2)));

            // transfer to cm system and make some rotations
            // boost them to centre of mass
            alpha_cm = alpha;
            carbon_cm = carbon;
            oxygen_cm = oxygen;
            alpha_cm.Boost(-oxygen.BoostVector());
            carbon_cm.Boost(-oxygen.BoostVector());
            oxygen_cm.Boost(-oxygen.BoostVector());
            cout << " px: " << alpha_cm.Px() + carbon_cm.Px() << "  py:  " << alpha_cm.Py() + carbon_cm.Py()
                 << "  pz:  " << alpha_cm.Pz() + carbon_cm.Pz() << endl;
            /*
                        cout << "lab " << alpha.Px() << "  " << alpha.Py() << "  " << alpha.Pz() << endl;
                        cout << "lab " << carbon.Px() << "  " << carbon.Py() << "  " << carbon.Pz() << endl;
                        cout << "lab " << oxygen.Px() << "  " << oxygen.Py() << "  " << oxygen.Pz() << endl;

                        cout << "cm " << alpha_cm.Px() << "  " << alpha_cm.Py() << "  " << alpha_cm.Pz() << endl;
                        cout << "cm " << carbon_cm.Px() << "  " << carbon_cm.Py() << "  " << carbon_cm.Pz() << endl;
                        cout << "cm " << oxygen_cm.Px() << "  " << oxygen_cm.Py() << "  " << oxygen_cm.Pz() << endl;
            */

            Double_t phi = (alpha_cm.Phi() - carbon_cm.Phi()) * TMath::RadToDeg();
            fh_phi26_cm->Fill(phi); // phi between alpha + carbon in cm
            Double_t theta = alpha_cm.Angle(carbon_cm.Vect()) * TMath::RadToDeg();
            fh_theta26_cm->Fill(theta); // angle between alpha + carbon in cm

            fh_theta_4He_cm->Fill(alpha_cm.Theta() * TMath::RadToDeg());
            fh_phi_4He_cm->Fill(alpha_cm.Phi() * TMath::RadToDeg());
            fh_theta_12C_cm->Fill(carbon_cm.Theta() * TMath::RadToDeg());
            fh_phi_12C_cm->Fill(carbon_cm.Phi() * TMath::RadToDeg());

            fh_theta_bc_cm->Fill(alpha_cm.Theta() * TMath::RadToDeg());
            Double_t phi_bc_cm = (oxygen.Phi() - alpha_cm.Phi()) * TMath::RadToDeg();
            if (phi_bc_cm < 0)
                phi_bc_cm += 360.;
            fh_phi_bc_cm->Fill(phi_bc_cm);
            fh_phi_bc_cm_polar->Fill(phi_bc_cm, 1);

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

    // second, old version of analysis
    if (is_carbon && is_alpha && 1 == 0)
    {
        fh_chiy_vs_chix_He->Fill(chiHex, chiHey);
        fh_chiy_vs_chix_C->Fill(chiCx, chiCy);

        if (chiHex < 1e10 && chiHey < 1e10 && chiCx < 1e10 && chiCy < 1e10)
        {
            fh_target_xy->Fill(XHe * 100., YHe * 100.);

            fh_px_He->Fill(pHex);
            fh_py_He->Fill(pHey);
            fh_pz_He->Fill(pHez);

            fh_px_C->Fill(pCx);
            fh_py_C->Fill(pCy);
            fh_pz_C->Fill(pCz);

            fh_px_px->Fill(pHex, pCx);
            fh_py_py->Fill(pHey, pCy);
            fh_pz_pz->Fill(pHez, pCz);

            // Calculate angle between alphs and C
            if (pCz == 0 || pHez == 0)
                return;

            costh26 = (pHex * pCx + pHey * pCy + pHez * pCz) /
                      (sqrt(pHex * pHex + pHey * pHey + pHez * pHez) * sqrt(pCx * pCx + pCy * pCy + pCz * pCz));

            theta_26 = acos(costh26) * 180. / 3.14159; // opening angle

            fh_theta26->Fill(theta_26);

            Double_t pHe = sqrt(pow(pHex, 2) + pow(pHey, 2) + pow(pHez, 2));
            Double_t pC = sqrt(pow(pCx, 2) + pow(pCy, 2) + pow(pCz, 2));
            Double_t eHe = sqrt(pow(pHe, 2) + pow(mHe, 2)) - mHe;
            Double_t eC = sqrt(pow(pC, 2) + pow(mC, 2)) - mC;

            fh_p_He->Fill(pHe);
            fh_p_C->Fill(pC);

            Double_t Erela = sqrt(pow((mHe + mC + eHe + eC), 2) - pow(pHe, 2) - pow(pC, 2) -
                                  2 * pHe * pC * cos(theta_26 * 3.1415 / 180.)) -
                             (mHe + mC); // Erel

            Double_t ptota = sqrt(pow(pHe, 2) + pow(pC, 2) + 2 * pHe * pC * cos(theta_26 * 3.1415 / 180.));

            Double_t m_inv = sqrt(mHe * mHe + mC * mC + 2. * sqrt(pC * pC + mC * mC) * sqrt(pHe * pHe + mHe * mHe) -
                                  2. * pHe * pC * cos(theta_26 * 3.1415 / 180.));
            Double_t Erelb = m_inv - mHe - mC;

            fh_Erel->Fill(Erela);
            if (pHex < 0. && pCx < 0.)
                fh_ErelR->Fill(Erela);
            if (pHex > 0. && pCx > 0.)
                fh_ErelL->Fill(Erela);
            if ((pHex > 0. && pCx < 0.) || (pHex < 0. && pCx > 0.))
                fh_ErelB->Fill(Erela);

            //	fh_dErel_vs_x->Fill(xTest * 100., Erela - 1.);
            //	fh_dErel_vs_y->Fill(yTest * 100., Erela - 1.);
            // cout << "Theta 26: " << theta_26 << " Erel: " << Erela << " " << Erelb << endl;

            // transfer to cm system and make some rotations

            px = pHex + pCx;
            py = pHey + pCy;
            pz = pHez + pCz;

            if (pz != 0.)
            {
                theta_16 = atan2(sqrt(pow(px, 2) + pow(py, 2)), pz) * 180. / 3.14159265;
            }

            fh_theta_16O->Fill(theta_16);

            Double_t n2[10], eaa, vaa, g2, v2[10], lt2[10][10];
            Double_t tmp1[10], tmp2[10], tmp3[10], tmp4[10];
            Double_t pxO_1, pyO_1, pzO_1, esum_3;
            Double_t psum_1, thcm, phcm, tanphicm;

            esum_3 = sqrt(pow((sqrt(pow(pHe, 2) + pow(pC, 2) + 2 * pHe * pC * cos(theta_26 * 3.14159 / 180.))), 2) +
                          pow(mO, 2)) -
                     mO;

            psum_1 = sqrt(pow(px, 2) + pow(py, 2) + pow(pz, 2));
            thcm = acos(pz / psum_1);
            tanphicm = py / px;

            if (px > 0)
            {
                phcm = atan(tanphicm) * 180. / 3.14159265;
            }

            if (px < 0. && py > 0.)
            {
                phcm = (atan(tanphicm) * 180. / 3.14159265) + 180.;
            }

            if (px < 0. && py < 0.)
            {
                phcm = (atan(tanphicm) * 180. / 3.14159265) - 180.;
            }

            if (phcm < 0)
            {
                phcm = phcm + 360.;
            }
            phcm = phcm * 3.14159265 / 180.;

            //            fh_theta_16O_cm->Fill(thcm);
            //            fh_phi_16O_cm->Fill(phcm);

            n2[1] = px / ptota;
            n2[2] = py / ptota;
            n2[3] = pz / ptota;

            //       eaa   = esum_3;
            //     velocity of 6Li
            //       vaa=sqrt(eaa*eaa+2*mO*eaa)/(eaa+mO);

            //       cout<< " va1 "<< vaa<<endl;

            eaa = sqrt(mO * mO + psum_1 * psum_1);
            vaa = psum_1 / eaa;

            //       cout<< " va2 "<< vaa<<endl;

            g2 = 1.0 / sqrt(1.0 - vaa * vaa);

            v2[1] = vaa * n2[1];
            v2[2] = vaa * n2[2];
            v2[3] = vaa * n2[3];

            // Lorentz transformation
            lt2[1][1] = g2;
            lt2[1][2] = -1. * g2 * v2[1];
            lt2[1][3] = -1. * g2 * v2[2];
            lt2[1][4] = -1. * g2 * v2[3];
            lt2[2][1] = lt2[1][2];
            lt2[2][2] = 1.0 + (g2 - 1.0) * n2[1] * n2[1];
            lt2[2][3] = (g2 - 1.0) * n2[1] * n2[2];
            lt2[2][4] = (g2 - 1.0) * n2[1] * n2[3];
            lt2[3][1] = lt2[1][3];
            lt2[3][2] = lt2[2][3];
            lt2[3][3] = 1.0 + (g2 - 1.0) * n2[2] * n2[2];
            lt2[3][4] = (g2 - 1.0) * n2[2] * n2[3];
            lt2[4][1] = lt2[1][4];
            lt2[4][2] = lt2[2][4];
            lt2[4][3] = lt2[3][4];
            lt2[4][4] = 1.0 + (g2 - 1.0) * n2[3] * n2[3];

            //      four momentum of nucleus b
            tmp1[1] = sqrt(mHe * mHe + pHe * pHe);
            tmp1[2] = pHex;
            tmp1[3] = pHey;
            tmp1[4] = pHez;

            //     four momentum of nucleus c
            tmp2[1] = sqrt(mC * mC + pC * pC);
            tmp2[2] = pCx;
            tmp2[3] = pCy;
            tmp2[4] = pCz;

            //     four momentum of nucleus a for a check of the transformation
            tmp4[1] = sqrt(mO * mO + psum_1 * psum_1);
            tmp4[2] = px;
            tmp4[3] = py;
            tmp4[4] = pz;
            //       cout <<"compare " << tmp4[1]<<"  "<< esum_3+mO<<endl;

            /*
                            cout <<"in " << tmp4[1]<<endl;
                            cout <<"in " << tmp4[2]<<endl;
                            cout <<"in " << tmp4[3]<<endl;
                            cout <<"in " << tmp4[4]<<endl;
            */
            for (Int_t i1 = 1; i1 < 5; i1++)
            {
                tmp3[i1] = 0.0;
                for (Int_t i2 = 1; i2 < 5; i2++)
                {
                    tmp3[i1] = tmp3[i1] + lt2[i1][i2] * tmp1[i2];
                }
            }

            Double_t pHex_1, pHey_1, pHez_1, pHe_1, pCx_1, pCy_1, pCz_1, pC_1;

            pHex_1 = tmp3[2];
            pHey_1 = tmp3[3];
            pHez_1 = tmp3[4];

            for (Int_t i1 = 1; i1 < 5; i1++)
            {
                tmp3[i1] = 0.0;
                for (Int_t i2 = 1; i2 < 5; i2++)
                {
                    tmp3[i1] = tmp3[i1] + lt2[i1][i2] * tmp2[i2];
                }
            }

            pCx_1 = tmp3[2];
            pCy_1 = tmp3[3];
            pCz_1 = tmp3[4];

            for (Int_t i1 = 1; i1 < 5; i1++)
            {
                tmp3[i1] = 0.0;
                for (Int_t i2 = 1; i2 < 5; i2++)
                {
                    tmp3[i1] = tmp3[i1] + lt2[i1][i2] * tmp4[i2];
                }
            }

            pxO_1 = tmp3[2];
            pyO_1 = tmp3[3];
            pzO_1 = tmp3[4];

            /*
                            cout <<"out " << tmp3[1]<<endl;
                            cout <<"out " << tmp3[2]<<endl;
                            cout <<"out " << tmp3[3]<<endl;
                            cout <<"out " << tmp3[4]<<endl;
            */

            Double_t pHex_cm, pHey_cm, pHez_cm, pHe_cm;
            Double_t pCx_cm, pCy_cm, pCz_cm, pC_cm;
            Double_t thet_cm_He, thet_cm_C, thet;
            Double_t theta_cm_He, theta_cm_C, theta;
            Double_t pxyd_cm, costhetacm_He, tanphi_d, phicm_He, pxyhe_cm;
            Double_t costhetacm_C, tanphi_he, phicm_C, theta_cm_26;
            Double_t phi_cm_26, pxO_cm, pyO_cm, pzO_cm;

            //     Vector for rotation
            n2[1] = py / sqrt(px * px + py * py);
            n2[2] = -px / sqrt(px * px + py * py);
            n2[3] = 0;

            //       cout <<"winkel " << thcm<<"  "<<theta6/180.*3.14159265<<endl;

            lt2[1][1] = cos(thcm) + n2[1] * n2[1] * (1. - cos(thcm));
            lt2[1][2] = n2[1] * n2[2] * (1. - cos(thcm)) - n2[3] * sin(thcm);
            lt2[1][3] = n2[1] * n2[3] * (1. - cos(thcm)) + n2[2] * sin(thcm);
            lt2[2][1] = n2[2] * n2[1] * (1. - cos(thcm)) + n2[3] * sin(thcm);
            lt2[2][2] = cos(thcm) + n2[2] * n2[2] * (1. - cos(thcm));
            lt2[2][3] = n2[2] * n2[3] * (1. - cos(thcm)) - n2[1] * sin(thcm);
            lt2[3][1] = n2[3] * n2[1] * (1. - cos(thcm)) - n2[2] * sin(thcm);
            lt2[3][2] = n2[3] * n2[2] * (1. - cos(thcm)) + n2[1] * sin(thcm);
            lt2[3][3] = cos(thcm) + n2[3] * n2[3] * (1. - cos(thcm));

            //      four momentum of nucleus b
            tmp1[1] = pHex_1;
            tmp1[2] = pHey_1;
            tmp1[3] = pHez_1;

            //     four momentum of nucleus c
            tmp2[1] = pCx_1;
            tmp2[2] = pCy_1;
            tmp2[3] = pCz_1;

            //     four momentum of nucleus a for a check of the transformation
            tmp4[1] = px;
            tmp4[2] = py;
            tmp4[3] = pz;
            //       cout <<"compare " << tmp4[1]<<"  "<< esum_3+mO<<endl;

            /*
                            cout <<"in1 " << tmp4[1]<<endl;
                            cout <<"in1 " << tmp4[2]<<endl;
                            cout <<"in1 " << tmp4[3]<<endl;
            */
            Double_t pHex_2, pHey_2, pHez_2, pHe_2, pCx_2, pCy_2, pCz_2, pC_2;
            Double_t pxO_2, pyO_2, pzO_2;

            for (Int_t i1 = 1; i1 < 4; i1++)
            {
                tmp3[i1] = 0.0;
                for (Int_t i2 = 1; i2 < 4; i2++)
                {
                    tmp3[i1] = tmp3[i1] + lt2[i1][i2] * tmp1[i2];
                }
            }

            pHex_cm = tmp3[1];
            pHey_cm = tmp3[2];
            pHez_cm = tmp3[3];

            for (Int_t i1 = 1; i1 < 4; i1++)
            {
                tmp3[i1] = 0.0;
                for (Int_t i2 = 1; i2 < 4; i2++)
                {
                    tmp3[i1] = tmp3[i1] + lt2[i1][i2] * tmp2[i2];
                }
            }

            pCx_cm = tmp3[1];
            pCy_cm = tmp3[2];
            pCz_cm = tmp3[3];

            for (Int_t i1 = 1; i1 < 4; i1++)
            {
                tmp3[i1] = 0.0;
                for (Int_t i2 = 1; i2 < 4; i2++)
                {
                    tmp3[i1] = tmp3[i1] + lt2[i1][i2] * tmp4[i2];
                }
            }

            pxO_cm = tmp3[1];
            pyO_cm = tmp3[2];
            pzO_cm = tmp3[3];

            /*
                            cout <<"out1 " << tmp3[1]<<endl;
                            cout <<"out1 " << tmp3[2]<<endl;
                            cout <<"out1 " << tmp3[3]<<endl;
            */

            // 4He fragment

            pHe_cm = sqrt(pHex_cm * pHex_cm + pHey_cm * pHey_cm + pHez_cm * pHez_cm);

            costhetacm_He = pHez_cm / pHe_cm;
            theta_cm_He = acos(costhetacm_He) * 180. / 3.14159265;
            tanphi_d = pHey_cm / pHex_cm;

            if (pHex_cm > 0)
            {
                phicm_He = atan(tanphi_d) * 180. / 3.14159265;
            }

            if (pHex_cm < 0. && pHey_cm > 0.)
            {
                phicm_He = (atan(tanphi_d) * 180. / 3.14159265) + 180.;
            }

            if (pHex_cm < 0. && pHey_cm < 0.)
            {
                phicm_He = (atan(tanphi_d) * 180. / 3.14159265) - 180.;
            }

            if (phicm_He < 0)
            {
                phicm_He = phicm_He + 360.;
            }
            fh_theta_4He_cm->Fill(theta_cm_He);
            fh_phi_4He_cm->Fill(phicm_He);

            // 12C fragment

            pC_cm = sqrt(pCx_cm * pCx_cm + pCy_cm * pCy_cm + pCz_cm * pCz_cm);

            costhetacm_C = pCz_cm / pC_cm;
            theta_cm_C = acos(costhetacm_C) * 180. / 3.14159265;
            tanphi_he = pCy_cm / pCx_cm;

            if (pCx_cm > 0.)
            {
                phicm_C = atan(tanphi_he) * 180. / 3.14159265;
            }

            if (pCx_cm < 0. && pCy_cm > 0.)
            {
                phicm_C = (atan(tanphi_he) * 180. / 3.14159265) + 180.;
            }

            if (pCx_cm < 0. && pCy_cm < 0.)
            {
                phicm_C = (atan(tanphi_he) * 180. / 3.14159265) - 180.;
            }

            if (phicm_C < 0.)
            {
                phicm_C = phicm_C + 360.;
            }
            fh_theta_12C_cm->Fill(theta_cm_C);
            fh_phi_12C_cm->Fill(phicm_C);

            theta_cm_26 = theta_cm_He + theta_cm_C;
            phi_cm_26 = phicm_He - phicm_C;

            fh_theta26_cm->Fill(theta_cm_26);
            fh_phi26_cm->Fill(phi_cm_26);

            fh_px_px_cm->Fill(pHex_cm, pCx_cm);
            fh_py_py_cm->Fill(pHey_cm, pCy_cm);
            fh_pz_pz_cm->Fill(pHez_cm, pCz_cm);

        } // end if chi2
    }     // end if trackHits>1

    fNEvents += 1;
    return;
}

void R3BGlobalAnalysisS454::FinishEvent()
{
    fTrack->Clear();
    if (fMCTrack)
        fMCTrack->Clear();
}

void R3BGlobalAnalysisS454::FinishTask()
{

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

    fh_theta26->Write();
    fh_theta26_cm->Write();
    fh_theta_16O->Write();
    fh_phi_16O->Write();
    fh_phi26_cm->Write();
    fh_theta_4He_cm->Write();
    fh_phi_4He_cm->Write();
    fh_theta_bc_cm->Write();
    fh_phi_bc_cm->Write();
    fh_phi_bc_cm_polar->Write();
    fh_theta_12C_cm->Write();
    fh_phi_12C_cm->Write();
    fh_Erel->Write();
    fh_ErelB->Write();
    fh_ErelR->Write();
    fh_ErelL->Write();
    fh_theta26_simu->Write();
    fh_Erel_simu->Write();
    fh_chiy_vs_chix_He->Write();
    fh_chiy_vs_chix_C->Write();
    fh_dErel_vs_x->Write();
    fh_dErel_vs_y->Write();
    fh_px_px->Write();
    fh_py_py->Write();
    fh_pz_pz->Write();
    fh_px_px_cm->Write();
    fh_py_py_cm->Write();
    fh_pz_pz_cm->Write();
}
ClassImp(R3BGlobalAnalysisS454)
