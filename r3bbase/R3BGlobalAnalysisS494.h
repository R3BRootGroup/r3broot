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
// -----            Created 13-04-2016 by M.Heil          -----
// -----               Fill online histograms             -----
// ------------------------------------------------------------

#ifndef R3BGLOBALANALYSISS494
#define R3BGLOBALANALYSISS494

#include "FairTask.h"
#include <array>
#include <fstream>
#include <iostream>
#include <sstream>

#include "TClonesArray.h"
#include "TMath.h"
#include "TCutG.h"
#include <TLorentzVector.h>
#include <cstdlib>

// For 16O@500AMeV
#define E_0 (15.99065084 * 931.494028)
#define E_Beam  (E_0 + 16. * 500.)
#define BETA  (sqrt(1. - (E_0 * E_0) / (E_Beam * E_Beam)))
#define GAMMA  (1. / (sqrt(1. - BETA * BETA)))
// 7 detectors in s494 fi23a,b,fi30-33,tofd
#define N_DET_MAX 7

class TClonesArray;
class TH1F;
class TH2F;
class R3BEventHeader;

/**
 * This taks reads all detector data items and plots histograms
 * for online checks.
 */
class R3BGlobalAnalysisS494 : public FairTask
{

  public:
    /**
     * Default constructor.
     * Creates an instance of the task with default parameters.
     */
    R3BGlobalAnalysisS494();

    /**
     * Standard constructor.
     * Creates an instance of the task.
     * @param name a name of the task.
     * @param iVerbose a verbosity level.
     */
    R3BGlobalAnalysisS494(const char* name, Int_t iVerbose = 1);

    /**
     * Destructor.
     * Frees the memory used by the object.
     */
    virtual ~R3BGlobalAnalysisS494();

    /**
     * Method for task initialization.
     * This function is called by the framework before
     * the event loop.
     * @return Initialization status. kSUCCESS, kERROR or kFATAL.
     */
    virtual InitStatus Init();

    /**
     * Method for event loop implementation.
     * Is called by the framework every time a new event is read.
     * @param option an execution option.
     */
    virtual void Exec(Option_t* option);

    /**
     * A method for finish of processing of an event.
     * Is called by the framework for each event after executing
     * the tasks.
     */
    virtual void FinishEvent();

    /**
     * Method for finish of the task execution.
     * Is called by the framework after processing the event loop.
     */
    virtual void FinishTask();


    /**
     * Method for setting the trigger value.
     * @param trigger 1 - onspill, 2 - offspill, -1 - all events.
     */
    inline void SetTrigger(Int_t trigger) { fTrigger = trigger; }
    inline void SetTpat(Int_t tpat) { fTpat = tpat; }

    /**
     * Methods for setting reset and readout times for Bmon
     */
    inline void SetBmon(Int_t sens_SEE, Int_t sens_IC)
    {
        fsens_SEE = sens_SEE;
        fsens_IC = sens_IC;
    }

    /**
     * Methods for setting cuts
     */
    inline void SetCuts(Bool_t cuts)
    {
        fCuts = cuts;
    }

    inline void SetGhost(Bool_t ghost)
    {
        fGhost = ghost;
    }

    inline void SetPairs(Bool_t p)
    {
        fPairs = p;
    }

    inline void SetSimu(Bool_t simu)
    {
        fSimu = simu;
    }

    inline void SetBfield(Int_t B)
    {
        fB = B;
    }
    inline void SetGraphicalCuts(Bool_t graphCuts)
    {
        fGraphCuts = graphCuts;
    }
    
    inline void SetChi2(Double_t cut_chiX, Double_t cut_chiY)
    {
        fcut_chiX = cut_chiX;
        fcut_chiY = cut_chiY;
    }
    inline void SetThetaGrazing(Double_t theta)
    {
		fThetaGrazing = theta;
	}
	inline void SetFib23Cuts(Double_t xfibcut, Double_t yfibcut)
	{
		fxfibcut = xfibcut;
		fyfibcut = yfibcut;
	}
	inline void SetVisual(Bool_t vis)
	{
		fvis = vis;
	}
	
	inline void SetEvsECut(TString file){fEvsE = file;}
	inline void SetCalifaCut(TString file){fCalifaCut = file;}
  
  private:
    TClonesArray* fMCTrack;
    TClonesArray* fTrack;
    TClonesArray* fMappedItemsCalifa;
    TClonesArray* fCalItemsCalifa;
    TClonesArray* fHitItemsCalifa;
    
    TClonesArray* fCalifaMappedItems;
    TClonesArray* fCalifaCalItems;
    TClonesArray* fCalifaHitItems; 
    Int_t fNofCalifaMappedItems = 0; 
    Int_t fNofCalifaCalItems = 0;     
    Int_t fNofCalifaHitItems = 0;
    
    const char* fDetectorNames[N_DET_MAX + 1] = {"Fi23a", "Fi23b",  "Fi30", "Fi31","Fi32", "Fi33", "Tofd", NULL };
    // check for trigger should be done globablly (somewhere else)
    R3BEventHeader* header; /**< Event header. */
    Int_t fTrigger;         /**< Trigger value. */
    Int_t fTpat;
	Bool_t fCuts;
	Bool_t fGhost;
	Bool_t fPairs;
	Bool_t fGraphCuts;
	Bool_t fSimu;
	Bool_t fvis;
	Int_t fB;
	Double_t fThetaGrazing;
	Bool_t tracker = true;
	Double_t fxfibcut, fyfibcut;
	TString fEvsE;
	TString fCalifaCut;
	
    unsigned long long time_start = 0, time = 0;
    unsigned long ic_start = 0, see_start = 0, tofdor_start = 0;
    unsigned long fNEvents = 0, fNEvents_start = 0, fNeventstrack = 0; /**< Event counter. */

    Int_t maxevent;

    Int_t fsens_SEE, fsens_IC; // SEETRAM and IC sensitivity, between -4 and -10
    Double_t calib_SEE = 1.;   // SEETRAM calibration factor
    Double_t see_offset = 7.1; // SEETRAM offset in kHz
    Double_t counts_SEE = 0;
    Double_t counts_IC = 0;
    Double_t counts_TofD = 0;

	Double_t XHe, YHe, ZHe, XC, YC, ZC, THe, TC;
	Double_t mtrackHe, mtrackC, mtrackO;	
	Double_t pHex, pHey, pHez, pCx, pCy, pCz;
	Double_t Pxf, Pyf, Pzf, Xf, Yf, Zf, Pf_tot;
	Double_t XHe_mc, YHe_mc, ZHe_mc, XC_mc, YC_mc, ZC_mc, THe_mc, TC_mc;	
	Double_t pHex_mc, pHey_mc, pHez_mc, pCx_mc, pCy_mc, pCz_mc;
	Double_t Pxf_mc, Pyf_mc, Pzf_mc, Xf_mc, Yf_mc, Zf_mc, Pf_tot_mc;

     Double_t fcut_chiX ;
     Double_t fcut_chiY ;	
	
	Double_t tofCalifa_shift=119.52; //ns

	Double_t Erel, m_inva, phi_bc_cm;
	Double_t ErelMC, m_invaMC, phiMC_bc_cm;
	TVector3 pa, pc, paMC, pcMC;
	TVector3 pa_cm, pc_cm, poexc, paMC_cm, pcMC_cm;
	TVector3 pbeam;
	ULong64_t timeTS ;
    ULong64_t eventTS ;
    
    TLorentzVector alpha, carbon, oxygen, helium3;
    TLorentzVector alphaMC, carbonMC, oxygenMC;
    TLorentzVector alpha_cm, carbon_cm, oxygen_cm;
    TLorentzVector alphaMC_cm, carbonMC_cm, oxygenMC_cm;
    TLorentzVector photon, photonMC;
    TLorentzVector oxygenBeam;

	Int_t Q = 0;
	Double_t tPrev[10];
	Int_t detPrev[10];
	
	Int_t counter1 = 0;
	Int_t counter2 = 0;
	Int_t counter3 = 0;
	Int_t counter4 = 0;
	Int_t countdet;
	TCutG *cut_EHe_EC;
	TCutG *cut_CalifaTof;
	Int_t clocal = 0;
	
    UInt_t num_spills = 0;

//    TH1F* fh_Tpat;
//    TH1F* fh_Trigger;
//    TH1F* fh_IC;
//    TH1F* fh_SEE;
 //   TH1F* fh_TOFDOR;
//    TH2F* fh_Cave_position;
    TH2F* fh_target_xy;
    TH2F* fh_target_xx;
    TH2F* fh_target_yy;
    TH1F* fh_chi2;
    TH1F* fh_px_He;
    TH1F* fh_py_He;
    TH1F* fh_pz_He;
    TH1F* fh_p_He;
    TH1F* fh_px_C;
    TH1F* fh_py_C;
    TH1F* fh_pz_C;
    TH1F* fh_p_C;
    TH1F* fh_px_O;
    TH1F* fh_py_O;
    TH1F* fh_pz_O;
    TH1F* fh_p_O;

    TH2F* fh_chiy_vs_chix_nc;
    TH2F* fh_chiy_vs_chix;
    TH2F* fh_phi26_vs_chi;
    TH2F* fh_psum_vs_theta26_nc;
    TH2F* fh_Erel_vs_psum;
    TH2F* fh_phiMC_bc_cm_polar;

    TH1F* fh_dx;
    TH1F* fh_dy;
    TH1F* fh_dz;
    TH1F* fh_dpxO;
    TH1F* fh_dpyO;
    TH1F* fh_dpzO;
    TH1F* fh_dpO;
    TH1F* fh_dpxC;
    TH1F* fh_dpyC;
    TH1F* fh_dpzC;
    TH1F* fh_dpC;
    TH1F* fh_dpxHe;
    TH1F* fh_dpyHe;
    TH1F* fh_dpzHe;
    TH1F* fh_dpHe;
    TH1F* fh_dp;
    TH2F* fh_thetax_dpx_C;
    TH2F* fh_thetay_dpy_C;
    TH2F* fh_thetax_dpx_He;
    TH2F* fh_thetay_dpy_He;
    TH2F* fh_x_dpx;
    TH2F* fh_y_dpy;
    TH2F* fh_thetax_dpy;
    TH2F* fh_thetay_dpx;
    TH2F* fh_dpy_dpx;
    TH2F* fh_px_px_mc;
    TH2F* fh_py_py_mc;
    TH2F* fh_pz_pz_mc;
    TH2F* fh_px_px;
    TH2F* fh_py_py;
    TH2F* fh_pz_pz;
    TH2F* fh_p_p;
    TH2F* fh_Erel_vs_phibcMC;
    TH2F* fh_Erel_vs_event;  
    
    TH2F* fh_thetax_px;
    TH2F* fh_thetay_py;
    TH2F* fh_thetax_py;
    TH2F* fh_thetay_px;
    
    TH2F* fh_pHe_vs_theta26;
    TH2F* fh_psum_vs_theta26;
    
    TH2F* fh_theta26_vs_chi;  
    TH2F* fh_Erel_vs_chi; 
    TH2F* fh_psum_vs_chi; 
    TH2F* fh_target_xy_MC;
    
    TH2F* fh_target_xy_nc;
	TH1F* fh_px_He_nc;
	TH1F* fh_py_He_nc;
    TH1F* fh_pz_He_nc;
    TH1F* fh_p_He_nc;           
	TH1F* fh_px_C_nc;
	TH1F* fh_py_C_nc;
	TH1F* fh_pz_C_nc;
	TH1F* fh_p_C_nc;
	TH1F* fh_px_O_nc;
	TH1F* fh_py_O_nc;
    TH1F* fh_pz_O_nc;
    TH1F* fh_p_O_nc;      	
	TH2F* fh_px_px_nc;
	TH2F* fh_py_py_nc;
	TH2F* fh_pz_pz_nc;
	TH1F* fh_psum_nc;
	TH2F* fh_psum_vs_event_nc;				
	TH1F* fh_theta26_nc;
	TH1F* fh_Erel_nc;
	TH1F* fh_ErelB_nc;
	TH1F* fh_theta_bc_cm_nc;
	TH1F* fh_phi_bc_cm_nc;
	TH2F* fh_p_vs_x0;
    TH2F* fh_psum_vs_x0;
	TH2F* fh_p_vs_y0;
    TH2F* fh_psum_vs_y0;

	TH1F* fh_theta26_simu;
	TH1F* fh_Erel_simu;
	TH1F* fh_theta26;
	TH1F* fh_theta_16O;	
	TH1F* fh_phi_16O;
	TH1F* fh_phi_bc_cm;
	TH1F* fh_phiMC_bc_cm;
	TH2F* fh_theta_bc_cm;
	TH1F* fh_Erel;
	TH1F* fh_ErelL;
	TH1F* fh_ErelR;
	TH1F* fh_ErelB;
	TH2F* fh_phi_bc_cm_polar;
	TH2F* fh_Erel_vs_theta16O;
	TH2F* fh_ErelB_vs_theta16O;
	TH2F* fh_Erel_vs_theta16O_3He12C;
	TH2F* fh_Erel_vs_thetaMC;
	TH2F* fh_Erel_vs_y0;
	TH2F* fh_Erel_vs_x0;
	TH2F* fh_Erel_vs_yfi23;
	TH2F* fh_Erel_vs_xfi23;
	TH1F* fh_psum;
	TH1F* fh_psum_MC;
	TH1F* fh_pzsum;
	TH1F* fh_pzsum_MC;
	TH1F* fh_pxsum_MC;
	TH1F* fh_pysum_MC;
	TH1F* fh_dErel;
	TH1F* fh_dtheta;
	TH2F* fh_py_pz;
	TH2F* fh_py_pz_MC;
	TH2F* fh_py_yfi23;
	TH2F* fh_py_yfi23_MC;
	TH2F* fh_px_xfi23;
	TH2F* fh_px_xfi23_MC;
	TH2F* fh_dp_2d;
	TH2F* fh_dpx_2d;
	TH2F* fh_dpy_2d;
	TH2F* fh_dpz_2d;
	TH2F* fh_psum_vs_dx0;
	TH2F* fh_pC_vs_dx0;
	TH2F* fh_pC_vs_dx0_mc;
	TH2F* fh_pC_vs_x0;
	TH2F* fh_pC_vs_x0_mc;
	TH2F* fh_pHe_vs_dx0;
	TH2F* fh_pHe_vs_dx0_mc;
	TH2F* fh_pHe_vs_x0;
	TH2F* fh_pHe_vs_x0_mc;
	TH2F* fh_dpC_vs_dx0;
	TH2F* fh_dpHe_vs_dx0;
	TH2F* fh_x0_vs_dx0;
	TH2F* fh_Erel_vs_px;
	TH2F* fh_Erel_vs_py;
	TH2F* fh_Erel_vs_dr_tofd;
	TH2F* fh_Erel_vs_dr_limit_tofd;
	TH1F* fh_Erel_limit_tofd;
	TH2F* fh_Erel_vs_dr_fi23;
	TH2F* fh_Erel_vs_dr_limit_fi23;
	TH1F* fh_Erel_limit_fi23;
	TH2F* fh_dr_vs_dr;
	TH2F* fh_dr_vs_dr_limit;
	
	
	TH2F* fh_ErelB_vs_phibc_bg;
	TH2F* fh_ErelB_vs_theta26_bg;
	TH2F* fh_erelB_vs_ptransHe_bg;
	TH2F* fh_erelB_vs_ptransC_bg;
	TH2F* fh_thetaB_4He_cm_bg;
	TH2F* fh_phiB_4He_cm_bg;
	TH2F* fh_thetaB_12C_cm_bg;
	TH2F* fh_phiB_12C_cm_bg  ;
	TH2F* fh_ErelB_vs_phibc_ag;
	TH2F* fh_ErelB_vs_theta26_ag;
	TH2F* fh_erelB_vs_ptransHe_ag;
	TH2F* fh_erelB_vs_ptransC_ag;
	TH2F* fh_thetaB_4He_cm_ag;
	TH2F* fh_phiB_4He_cm_ag;
	TH2F* fh_thetaB_12C_cm_ag;
	TH2F* fh_phiB_12C_cm_ag;
	TH2F* fh_Erel_vs_phibc_bg;
	//TH2F* fh_Erel_vs_theta26;
	//TH2F* fh_Erel_vs_theta26_max;
	TH2F* fh_Erel_vs_theta26_bg;
	TH2F* fh_erel_vs_ptransHe_bg;
	TH2F* fh_erel_vs_ptransC_bg;
	TH2F* fh_theta_4He_cm_bg;
	TH2F* fh_phi_4He_cm_bg;
	TH2F* fh_theta_12C_cm_bg;
	TH2F* fh_phi_12C_cm_bg;
	TH2F* fh_Erel_vs_phibc_ag;
	TH2F* fh_Erel_vs_theta26_ag;
	TH2F* fh_erel_vs_ptransHe_ag;
	TH2F* fh_erel_vs_ptransC_ag;
	TH2F* fh_theta_4He_cm_ag;
	TH2F* fh_phi_4He_cm_ag;
	TH2F* fh_theta_12C_cm_ag;
	TH2F* fh_phi_12C_cm_ag;
	
    TH2F* fh_ErelB_vs_phibc_MC;
    TH2F* fh_erelB_vs_ptransHe_MC;
    TH2F* fh_erelB_vs_ptransC_MC;
    TH2F* fh_thetaB_4He_cm_MC;
    TH2F* fh_phiB_4He_cm_MC;
    TH2F* fh_thetaB_12C_cm_MC;
    TH2F* fh_phiB_12C_cm_MC;
    TH2F* fh_ErelB_vs_theta16O_MC;
	
	TH2F* fh_mass_nc;
	TH2F* fh_mass;
	TH2F* fh_energy_nc;
	TH2F* fh_energy;
	    
    TH1F* fh_califa_hitenergy_boost;
    TH1F* fh_califa_barrel_hitenergy_boost;
    TH1F* fh_califa_iphos_hitenergy_boost;
    TH1F* fh_Erel_withCalifa;
    TH2F* fh_califaenergy_2d;
    TH2F* fh_califaenergy_barrel_2d;
    TH2F* fh_califaenergy_iphos_2d;
    TH2F* fh_Erel_withCalifa_2d;
    TH2F* fh_Erel_vs_theta16O_withcalifa;
    TH2F* fh_califa_hitenergy_bg;
    TH2F* fh_califa_hitenergy_ag; 
    TH2F* fh_califa_hitenergy_rand_bg;
    TH2F* fh_califa_hitenergy_rand_ag; 
    TH2F* fh_califa_time_cId;
	
	TH2F* fh_califa_tofd;
	TH2F* fh_califa_tofd_cut;
	TH1F* fh_califa_hitenergy_select;
	TH1F* fh_califa_hitenergy_nc;
	TH2F* fh_Ecalifa_vs_theta;
	TH2F* fh_Erel_withCalifa_tof;
	TH2F* fh_Erel_withCalifa_motherId;
	TH2F* fh_phi_vs_ecalifa_barrel;
	TH2F* fh_phi_vs_ecalifa_iphos;
	TH2F* fh_califa_tofd_ecr;
	TH2F* fh_califa_ecl;
	TH2F* fh_califa_dtime_ecr;
	TH2F* fh_califa_dtof_ecr;
	TH2F* fh_esumcluster_dc_2d;
	TH1F* fh_califa_hitenergycorr_boost;
	TH2F* fh_crystalNb;
	TH1F* fh_califa_hitenergy_boost_rand;
	TH2F* fh_Erel_vs_theta16O_withcalifa_rand;
	
	TH1F* fh_minv_simu;
	TH1F* fh_minv;
	
	TH2F* fh_psum_vs_xdet[N_DET_MAX];		                                		
	TH2F* fh_p_vs_xdet[N_DET_MAX];		                                		
	TH2F* fh_psum_vs_ydet[N_DET_MAX];  
	TH2F* fh_yfi23_vs_ytofd_bc;
	TH2F* fh_yfi23_vs_ytofd;
	TH2F* fh_xfi23_vs_xtofd_bc;
	TH2F* fh_xfi23_vs_xtofd;
	TH2F* fh_xy_fib23;
	TH2F* fh_xx_fib23;
	TH2F* fh_yy_fib23;
	TH2F* fh_xy_fib23_bc;
	TH1F* fh_dt_fib23_bc;
	
	
	
  public:
    ClassDef(R3BGlobalAnalysisS494, 1)
};

#endif
