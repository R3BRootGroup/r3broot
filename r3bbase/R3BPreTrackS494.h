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
// -----                  R3BPreTrackS494                 -----
// -----            Created 04.01.2022 by A.Kelic-Heil    -----
// -----               Prepare hit level data for tracker -----
// ------------------------------------------------------------

#ifndef PRER3BTRACKS494
#define PRER3BTRACKS494
#define N_PLANE_MAX_TOFD 4
#define N_PADDLE_MAX_TOFD 50
#define N_PADDLE_MAX_PTOF 100
#define N_PSPX_S494 4

#include "FairTask.h"
#include <array>
#include <fstream>
#include <iostream>
#include <sstream>
#include "TCutG.h"

#include "TClonesArray.h"
#include "TMath.h"
#include <cstdlib>

class TClonesArray;
class TH1F;
class TH2F;
class R3BEventHeader;

/**
 * This taks reads all detector data items and plots histograms
 * for online checks.
 */
class R3BPreTrackS494 : public FairTask
{

  public:
    /**
     * Default constructor.
     * Creates an instance of the task with default parameters.
     */
    R3BPreTrackS494();

    /**
     * Standard constructor.
     * Creates an instance of the task.
     * @param name a name of the task.
     * @param iVerbose a verbosity level.
     */
    R3BPreTrackS494(const char* name, Int_t iVerbose = 1);

    /**
     * Destructor.
     * Frees the memory used by the object.
     */
    virtual ~R3BPreTrackS494();

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
     inline void SetTpat(Int_t tpat1, Int_t tpat2) { fTpat1 = tpat1;  fTpat2 = tpat2;}

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

    inline void SetGraphicalCuts(Bool_t graphCuts)
    {
        fGraphCuts = graphCuts;
    }

    inline void SetPairs(Bool_t p)
    {
        fPairs = p;
    }

    inline void SetBfield(Int_t B)
    {
        fB = B;
    }
    inline void SetSimu(Int_t simu)
    {
        fSimu = simu;
    }
    inline void SetTracker(Bool_t trackerType)
    {
        ftrackerType = trackerType;   // 0 == Rene, 1 == Dima
    }
    inline void SetAverage(Bool_t averageValues)
    {
		fAverage = averageValues;
	}
	inline void SetIdeal(Bool_t idealData)
    {
		fidealData = idealData;
	}
	inline void SetXlimits(Double_t xmin1, Double_t xmax1, Double_t xmin2, Double_t xmax2){
		fX1min=xmin1;
		fX1max=xmax1;
		fX2min=xmin2;
		fX2max=xmax2;
	}
	inline void SetFi30Fi32Cut(TString file){fFi30Fi32Cut = file;}
	inline void SetFi31Fi33Cut(TString file){fFi31Fi33Cut = file;}
	inline void SetTofiCut(TString file){fTofiCut = file;}
    inline void SetHisto(Bool_t histo)
    {
        fHisto = histo;
    }
	    
  private:
    std::vector<TClonesArray*> fMappedItems;
    std::vector<TClonesArray*> fCalItems;
    std::vector<TClonesArray*> fHitItems;
    TClonesArray* fMCTrack;
    
    TClonesArray* fTofdHitItems;       
    Int_t fNofTofdHitItems = 0;
    TClonesArray* fFi23aHitItems;       
    Int_t fNofFi23aHitItems = 0;
    TClonesArray* fFi23bHitItems;       
    Int_t fNofFi23bHitItems = 0;
    TClonesArray* fFi30HitItems;       
    Int_t fNofFi30HitItems = 0;
    TClonesArray* fFi31HitItems;       
    Int_t fNofFi31HitItems = 0;
    TClonesArray* fFi32HitItems;       
    Int_t fNofFi32HitItems = 0;
    TClonesArray* fFi33HitItems;       
    Int_t fNofFi33HitItems = 0;
    TClonesArray* fCalifaHitItems;       
    Int_t fNofCalifaHitItems = 0;

    enum DetectorInstances
    {
        DET_CALIFA,
        DET_BMON,
        DET_ROLU,
        DET_FI_FIRST,
        DET_FI23A = DET_FI_FIRST,
        DET_FI23B,
        DET_FI30,
        DET_FI31,
        DET_FI32,
        DET_FI33,
        DET_FI_LAST = DET_FI33,
        DET_TOFD,
        DET_MAX
    };

#define NOF_FIB_DET (DET_FI_LAST - DET_FI_FIRST + 1)

    const char* fDetectorNames[DET_MAX + 1] = { "Califa", "BeamMonitor", "Rolu", "Fi23a", "Fi23b", "Fi30",
                                                "Fi31",   "Fi32",        "Fi33", "Tofd", NULL };

    // If FiberI is present or not:
    Int_t ifibdet;
    // Number of fibers per detector
    Double_t n_fiber[NOF_FIB_DET] = { 512., 512., 2048., 2048., 2048., 2048. };

    // check for trigger should be done globablly (somewhere else)
    R3BEventHeader* header; /**< Event header. */
    Int_t fTrigger;         /**< Trigger value. */
    Int_t fTpat1,fTpat2;
	Bool_t fCuts;
	Bool_t fPairs;
	Bool_t fGraphCuts;
	Bool_t fSimu;
	Bool_t fidealData;
	Int_t fB;
	Bool_t tracker = true;
	Double_t fX1min, fX1max, fX2min, fX2max;
	Double_t delta;
	Bool_t ftrackerType;
	Bool_t fAverage;
	Bool_t fHisto;
	TString fFi30Fi32Cut, fFi31Fi33Cut, fTofiCut;
    
	unsigned long IcountwriteOut1[23]={0}, IcountwriteOut2[23]={0}, IcountwriteOut1mem = 0, IcountwriteOut2mem=0,countdet_written[10]={0};
	

	TCutG *cut_qfi31_qfi33;
	TCutG *cut_qfi30_qfi32;
	TCutG *cut_qtofi_bar14;
	
	
    unsigned long long time_start = 0, time = 0;
    unsigned long long ic_start = 0, see_start = 0, tofdor_start = 0;
    unsigned long fNEvents = 0, fNEvents_start = 0, fNeventstrack=0, fNeventselect=0; /**< Event counter. */

    Int_t maxevent;
    Int_t ncorellx = 0, nanticorellx = 0;
    Int_t ncorelly = 0, nanticorelly = 0;

    Int_t fsens_SEE, fsens_IC; // SEETRAM and IC sensitivity, between -4 and -10
    Double_t calib_SEE = 1.;   // SEETRAM calibration factor
    Double_t see_offset = 7.1; // SEETRAM offset in kHz
    Double_t counts_SEE = 0;
    Double_t counts_IC = 0;
    Double_t counts_TofD = 0;
    
    Double_t ytofd_offsetZ2[2][44]={0},ytofd_offsetZ6[2][44]={0},ytofd_offsetZ8[2][44]={0};
    
    Double_t qdet_sum[10] ={ 0};
    Double_t xdet_sum[10] ={ 0};
    Double_t ydet_sum[10] ={ 0};
    Double_t tdet_sum[10] ={ 0};
    Double_t zdet_sum[10] ={ 0};
   
    Int_t nsum[10] = { 0 };

	Double_t XHes, YHes, ZHes, XCs, YCs, ZCs, THes, TCs;
	Double_t pHexs, pHeys, pHezs, pCxs, pCys, pCzs, pHes, pCs;
	
	Double_t amu = 0.931494028e3; //931.49410242;
	Double_t pHex, pHey, pHez, pCx, pCy, pCz;
	Double_t Pxf, Pyf, Pzf, Xf, Yf, Zf, Pf_tot;
//	Double_t mHe = 4.00260325413*amu;
//	Double_t mC = 12. * amu;
	Double_t mHe = 3.7273791;//3727.409;
	Double_t mC = 11.174862;//11174.950;
                           
	
	Int_t Q = 0;
	
	Int_t counter1 = 0, counterwo=0;
	Int_t counter1mem = 0;
	Int_t counter2 = 0;
	Int_t counter3 = 0;
	Int_t counter4 = 0;
	Int_t counterTofd = 0;
	Int_t counterTofdMulti = 0;
	Int_t counterCalifa = 0;
	Int_t counterWrongTpat = 0;
	Int_t counterWrongTrigger = 0;
	Int_t counterRolu = 0;
	Int_t counterTracker = 0;
	Int_t countdet;
	Int_t countdet_s;
	Int_t countdet_ss;
	Int_t countdet50 = 0;
	Double_t hits1 = 0;
	Double_t hits30 = 0;
	Double_t hits30bc = 0;
	Double_t hits31 = 0;
	Double_t hits31bc = 0;
	Double_t hits32 = 0;
	Double_t hits32bc = 0;
	Double_t hits33 = 0;
	Double_t hits33bc = 0;
	Bool_t increment1 = false;
    Bool_t increment2 = false;
    Int_t Nevent_singles=0;  
    Int_t Nevent_wrongcharge =0;
    Int_t fNeventpair=0;
    Int_t Nevent_zerofi23a_corel = 0;
    Int_t Nevent_notzerofi23a_corel = 0;
    Int_t Nevent_zerofi23b_corel = 0;
    Int_t Nevent_notzerofi23b_corel = 0;
    Int_t Nevent_zerofi3032_corel = 0;
    Int_t Nevent_notzerofi3032_corel = 0;
    Int_t Nevent_zerofi3133_corel = 0;
    Int_t Nevent_notzerofi3133_corel = 0;
    Int_t Nevent_goodQtofd = 0, Nevent_goodQ = 0, Nevent_zerofi3x_corel =0,Neventafterfibcor=0;
    Int_t Nhits_before_fibcor = 0, Nhits_after_fibcor = 0;

    UInt_t num_spills = 0;

	Int_t ndet = 10;
	
    TH1F* fh_Tpat;
    TH1F* fh_Trigger;
    TH1F* fh_IC;
    TH1F* fh_SEE;
    TH1F* fh_TOFDOR;

    TH2F* fh_califa_energy;
	TH2F* fh_califa_energy_select_barrel;
	TH2F* fh_califa_energy_select_iphos;
	TH2F* fh_califa_tofd_barrel;
	TH2F* fh_califa_tofd_iphos;
	TH2F* fh_califa_time;
	TH2F* fh_califa_energy_dc;
	TH2F* fh_califa_cluster;
	TH2F* fh_califa_overflow_rolu;
	TH2F* fh_califa_overflow_antirolu;
	TH2F* fh_califa_cluster_rolu;
	TH2F* fh_califa_overflow;
	TH2F* fh_califa_good;
	TH2F* fh_califa_good_rolu;
	TH2F* fh_califa_good_antirolu;
	TH2F* fh_califa_motherid_clustertype;
	
    TH2F* fh_xy_Fib[NOF_FIB_DET];
    TH2F* fh_xy_Fib_ac[NOF_FIB_DET];
    TH1F* fh_mult_Fib[NOF_FIB_DET];
    TH1F* fh_mult_Fib_ac[NOF_FIB_DET];
    TH2F* fh_Fib_ToF[NOF_FIB_DET];
    TH2F* fh_Fib_ToF_ac[NOF_FIB_DET];
    TH2F* fh_Fib_Time[NOF_FIB_DET];
    TH2F* fh_Fib_Time_ac[NOF_FIB_DET];
    TH2F* fh_ToT_Fib[NOF_FIB_DET];
    TH2F* fh_ToT_Fib_ac[NOF_FIB_DET];
    TH2F* fh_Fib_vs_Events[NOF_FIB_DET];
    TH2F* fh_Fib_vs_Events_ac[NOF_FIB_DET];
    TH2F* fh_Fibs_vs_Tofd[NOF_FIB_DET];
    TH2F* fh_Fibs_vs_Tofd_ac[NOF_FIB_DET];
    TH2F* fh_ToF_vs_Events[NOF_FIB_DET];
    TH2F* fh_ToF_vs_Events_ac[NOF_FIB_DET];
    TH2F* fh_qtofd_vs_qFib[NOF_FIB_DET];
    TH2F* fh_qtofd_vs_qFib_ac[NOF_FIB_DET];
    TH2F* fh_Fib_ToF_vs_Qtofd[NOF_FIB_DET];
    TH2F* fh_Fib32_vs_Fib30_tot;
    TH2F* fh_Fib33_vs_Fib31_tot;
    
    TH2F* fh_Fi23aToF_Q[21];
    TH2F* fh_Fi23bToF_Q[21];
    
    
    
	TH2F* fh_Fib30_vs_Fib23a;
	TH2F* fh_Fib30_vs_Fib23a_dx;
	TH2F* fh_Fib30_vs_Fib23a_dt;
	TH2F* fh_tofdright_vs_Fib23b;
	TH2F* fh_Fib31_vs_Fib23b_dx;
	TH2F* fh_Fib31_vs_Fib23b_dt;

	TH2F* fh_Fib33_vs_Fib31;
	TH2F* fh_Fib33_vs_Fib31_dx;
	TH2F* fh_Fib33_vs_Fib31_dt;
	TH2F* fh_Fib31_vs_Fib23a;
	TH2F* fh_Fib31_vs_Fib23a_dx;
	TH2F* fh_Fib31_vs_Fib23a_dt;
	TH2F* fh_Fib32_vs_Fib30;
	TH2F* fh_Fib32_vs_Fib30_dx;
	TH2F* fh_Fib32_vs_Fib30_dt;
	TH2F* fh_tofdleft_vs_Fib23b;
	TH2F* fh_Fib30_vs_Fib23b_dx;
	TH2F* fh_Fib30_vs_Fib23b_dt;
	TH2F* fh_Fib23a_vs_Fib3X_tot;
	TH2F* fh_Fib23b_vs_Fib3X_tot;	
	TH2F* fh_tofd_vs_Fib31;
    TH2F* fh_ytofd_vs_yFib31;
    TH2F* fh_tofd_vs_Fib31_ac;
    TH2F* fh_tofd_vs_Fib30;
    TH2F* fh_ytofd_vs_yFib30;
    TH2F* fh_tofd_vs_Fib30_ac;

	TH2F* fh_Fib33_vs_Fib31_ac;
	TH2F* fh_Fib33_vs_Fib31_dx_ac;
	TH2F* fh_Fib31_vs_Fib23a_ac;
	TH2F* fh_Fib31_vs_Fib23a_dx_ac;
	TH2F* fh_tofdright_vs_Fib23b_ac;
	TH2F* fh_Fib31_vs_Fib23b_dx_ac;
	TH2F* fh_Fib32_vs_Fib30_ac;
	TH2F* fh_Fib32_vs_Fib30_dx_ac;
	TH2F* fh_tofdleft_vs_Fib23b_ac;
	TH2F* fh_Fib30_vs_Fib23b_dx_ac;
	TH2F* fh_Fib30_vs_Fib23a_ac;
	TH2F* fh_Fib30_vs_Fib23a_dx_ac;
	TH2F* fh_Fi23aX_TofdX_He;
	TH2F* fh_Fi23aX_TofdX_C;
	TH2F* fh_Fi23bY_TofdY_He;
	TH2F* fh_Fi23bY_TofdY_C;
	TH2F* fh_Fi301X_TofdX_He;
	TH2F* fh_Fi301X_TofdX_C;
	TH2F* fh_Fi301Y_TofdY_He;
	TH2F* fh_Fi301Y_TofdY_C;
	TH2F* fh_Fi323X_TofdX_He;
	TH2F* fh_Fi323X_TofdX_C;
	TH2F* fh_Fi323Y_TofdY_He;
	TH2F* fh_Fi323Y_TofdY_C;
	TH2F* fh_Fi23aX_TofdX_bw1;
	TH2F* fh_Fi23aX_TofdX_bw2;
	TH2F* fh_yTofd_vs_yTofd;
	TH2F* fh_yTofdHe_vs_yTofdC   ;
	TH2F* fh_xTofdHe_vs_xTofdC;
	TH2F* fh_xy_target;
	TH2F* fh_TimePreviousEventFi23a;
    TH2F* fh_TimePreviousEventFi23b;

    TH2F* fh_Cave_position;

    TH2F* fh_xy_tofd;
    TH2F* fh_xy_tofd_ac;
    TH2F* fh_tofd_charge;
    TH2F* fh_tofd_charge_ac;
    TH1F* fh_tofd_time;
    TH1F* fh_tofd_time_ac;
    TH2F* fh_TimePreviousEvent;
    TH1F* fh_tofd_mult;
    TH1F* fh_tofd_mult_ac;
	TH2F* fh_tofd_q2_vs_q1;
	TH2F* fh_tofd_q2_vs_q1_ac;
	TH2F* fh_check_QvsX[7];
	TH2F* fh_check_TvsX[7];
	TH2F* fh_check_XvsY[7];
	TH2F* fh_tofd_x_vs_y_z[6];
	TH2F* fh_qsum_mult_fi23a;
	TH2F* fh_qsum_mult_fi23b;
	TH2F* fh_qtofi_vs_qtofd;


	
  public:
    ClassDef(R3BPreTrackS494, 1)
};

#endif
