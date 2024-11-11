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
// -----          Created 04.01.2022 by A.Kelic-Heil      -----
// ------------------------------------------------------------

/*
 * This task should fill histograms with detector variables which allow
 * to test the detectors online
 *
 */

#include "R3BCalifaClusterData.h"
#include "R3BCalifaCrystalCalData.h"
#include "R3BCalifaMappedData.h"
#include "R3BLosCalData.h"
#include "R3BLosHitData.h"
#include "R3BLosMappedData.h"

#include "R3BBeamMonitorMappedData.h"

#include "R3BPreTrackS494.h"

#include "R3BSci8CalData.h"
#include "R3BSci8MappedData.h"

#include "R3BTofdCalData.h"
#include "R3BTofdHitData.h"
#include "R3BTofdMappedData.h"

#include "R3BTofiCalData.h"
#include "R3BTofiHitData.h"
#include "R3BTofiMappedData.h"

#include "R3BFiberMAPMTCalData.h"
#include "R3BFiberMAPMTHitData.h"
#include "R3BFiberMAPMTMappedData.h"

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

#include "R3BWRData.h"

#include "R3BMCTrack.h"
#include "R3BTrack.h"

#include "FairLogger.h"
#include "FairRootManager.h"
#include "FairRunAna.h"
#include "FairRunOnline.h"
#include "FairRuntimeDb.h"
#include "R3BLogger.h"

#include "TCanvas.h"
#include "TH1F.h"
#include "TH2F.h"

#include "TCutG.h"
#include "tracker_routines.h"

#include "TClonesArray.h"
#include "TFile.h"
#include "TMath.h"
#include "TRandom.h"
#include <TRandom3.h>
#include <TRandomGen.h>
#include <array>
#include <cstdlib>
#include <ctime>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#define IS_NAN(x) TMath::IsNaN(x)
using namespace std;

R3BPreTrackS494::R3BPreTrackS494()
    : R3BPreTrackS494("PreTrack", 1)
{
}

R3BPreTrackS494::R3BPreTrackS494(const char* name, Int_t iVerbose)
    : FairTask(name, iVerbose)
    , fTrigger(-1)
    , fTpat1(-1)
    , fTpat2(-1)
    , fCuts(0)
    , fGraphCuts(0)
    , fPairs(0)
    , fB(-1250)
    , fSimu(0)
    , ftrackerType(0)
    , fAverage(0)
    , fidealData(1)
    , fX1min(-100)
    , fX1max(100)
    , fX2min(-100)
    , fX2max(100)
    , fNEvents(0)
    , fFi31Fi33Cut("")
    , fFi30Fi32Cut("")
    , fTofiCut("")
    , fHisto(1)
    , fTofdHitItems(new TClonesArray("R3BTofdHitData"))
    , fFi23aHitItems(new TClonesArray("R3BFiberMAPMTHitData"))
    , fFi23bHitItems(new TClonesArray("R3BFiberMAPMTHitData"))
    , fFi30HitItems(new TClonesArray("R3BFiberMAPMTHitData"))
    , fFi31HitItems(new TClonesArray("R3BFiberMAPMTHitData"))
    , fFi32HitItems(new TClonesArray("R3BFiberMAPMTHitData"))
    , fFi33HitItems(new TClonesArray("R3BFiberMAPMTHitData"))
    , fCalifaHitItems(new TClonesArray("R3BCalifaClusterData"))
    , fNofTofdHitItems(0)
    , fNofFi23aHitItems(0)
    , fNofFi23bHitItems(0)
    , fNofFi30HitItems(0)
    , fNofFi31HitItems(0)
    , fNofFi32HitItems(0)
    , fNofFi33HitItems(0)
    , fNofCalifaHitItems(0)
{
}

R3BPreTrackS494::~R3BPreTrackS494()
{
    for (int i = 0; i < NOF_FIB_DET; i++)
    {
        delete fh_xy_Fib[i];
        delete fh_mult_Fib[i];
        delete fh_ToT_Fib[i];
    }
}

InitStatus R3BPreTrackS494::Init()
{

    // Initialize random number:
    std::srand(std::time(0)); // use current time as seed for random generator

    LOG(info) << "R3BPreTrackS494::Init ";

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
    if (fMCTrack)
        mgr->Register("MCTrack", "Monte Carlo Tracks", fMCTrack, kTRUE);

    assert(DET_MAX + 1 == sizeof(fDetectorNames) / sizeof(fDetectorNames[0]));

    printf("Have %d detectors.\n", DET_MAX);
    printf("Have %d fiber detectors.\n", NOF_FIB_DET);

    for (int det = 0; det < DET_MAX; det++)
    {
        cout << "Reading detector " << det << ", " << fDetectorNames[det] << endl;

        if (det == 0) // CALIFA
            fMappedItems.push_back((TClonesArray*)mgr->GetObject(Form("%sMappedData", fDetectorNames[0])));
        else
            fMappedItems.push_back((TClonesArray*)mgr->GetObject(Form("%sMapped", fDetectorNames[det])));
        if (NULL == fMappedItems.at(det))
        {
            printf("Could not find mapped data for '%s'.\n", fDetectorNames[det]);
        }
        if (det == 9)
            maxevent = mgr->CheckMaxEventNo();
        if (det == 0) // CALIFA
            fCalItems.push_back((TClonesArray*)mgr->GetObject(Form("%sCrystalCalData", fDetectorNames[0])));
        else
            fCalItems.push_back((TClonesArray*)mgr->GetObject(Form("%sCal", fDetectorNames[det])));
        if (NULL == fCalItems.at(det))
        {
            printf("Could not find Cal data for '%s'.\n", fDetectorNames[det]);
        }
        if (det == 0) // CALIFA
            fHitItems.push_back((TClonesArray*)mgr->GetObject(Form("%sClusterData", fDetectorNames[0])));
        else
            fHitItems.push_back((TClonesArray*)mgr->GetObject(Form("%sHit", fDetectorNames[det])));

        if (NULL == fHitItems.at(det))
        {
            printf("Could not find hit data for '%s'.\n", fDetectorNames[det]);
        }
    }

    mgr->Register("TofdHit", "Land", fTofdHitItems, kTRUE);
    mgr->Register("Fi23aHit", "Land", fFi23aHitItems, kTRUE);
    mgr->Register("Fi23bHit", "Land", fFi23bHitItems, kTRUE);
    mgr->Register("Fi30Hit", "Land", fFi30HitItems, kTRUE);
    mgr->Register("Fi31Hit", "Land", fFi31HitItems, kTRUE);
    mgr->Register("Fi32Hit", "Land", fFi32HitItems, kTRUE);
    mgr->Register("Fi33Hit", "Land", fFi33HitItems, kTRUE);
    if (fHitItems.at(DET_CALIFA) && !fSimu)
        mgr->Register("CalifaClusterData", "Land", fCalifaHitItems, kTRUE);

    // *******
    //------------------------------------------------------------------------
    // graphical cuts
    //------------------------------------------------------------------------
    if (fGraphCuts)
    {
        cout << "READING GRAPHICAL CUTS " << endl;
        cut_qfi31_qfi33 = NULL;
        cut_qfi30_qfi32 = NULL;
        cut_qtofi_bar14 = NULL;

        if (fFi31Fi33Cut)
        {
            TFile* f31_33 = TFile::Open(fFi31Fi33Cut, "read");
            cut_qfi31_qfi33 = dynamic_cast<TCutG*>(f31_33->Get("CutFi31Fi33_12C"));
        }
        else
            R3BLOG(warn, "fFi31Fi33Cut file not defined");

        if (fFi30Fi32Cut)
        {
            TFile* f30_32 = TFile::Open(fFi30Fi32Cut, "read");
            cut_qfi30_qfi32 = dynamic_cast<TCutG*>(f30_32->Get("CutFi30Fi32_12C"));
        }
        else
            R3BLOG(warn, "fFi30Fi32Cut file not defined");

        if (fTofiCut)
        {
            TFile* qtofi = TFile::Open(fTofiCut, "read");
            cut_qtofi_bar14 = dynamic_cast<TCutG*>(qtofi->Get("qtofi_bar14"));
        }
        else
            R3BLOG(warn, "fTofiCut file not defined");

        cout << "GRAPHICAL CUTS ARE READ" << endl;
    }
    //------------------------------------------------------------------------
    // create histograms of all detectors
    //------------------------------------------------------------------------

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
    // Fiber Detectors 1-NOF_FIB_DET

    char canvName[255];
    UInt_t Nmax = 1e7;
    for (Int_t ifibcount = 0; ifibcount < NOF_FIB_DET; ifibcount++)
    {
        if (fHitItems.at(DET_FI_FIRST + ifibcount))
        {

            const char* detName;
            const char* detName2;
            detName = fDetectorNames[DET_FI_FIRST + ifibcount];

            LOG(debug) << "I am creating canvas " << detName << endl;

            // xy:
            if (ifibcount == 0)
            {
                fh_xy_Fib[ifibcount] =
                    new TH2F(Form("%s_xy", detName), Form("%s xy", detName), 429, -6.006, 6.006, 1000, -50., 50);
                fh_xy_Fib[ifibcount]->GetXaxis()->SetTitle("x / cm ");
                fh_xy_Fib[ifibcount]->GetYaxis()->SetTitle("y / cm");

                fh_xy_Fib_ac[ifibcount] = new TH2F(
                    Form("%s_xy_ac", detName), Form("%s xy after cuts", detName), 429, -6.006, 6.006, 1000, -50., 50);
                fh_xy_Fib_ac[ifibcount]->GetXaxis()->SetTitle("x / cm ");
                fh_xy_Fib_ac[ifibcount]->GetYaxis()->SetTitle("y / cm");
            }
            else if (ifibcount == 1)
            {
                fh_xy_Fib[ifibcount] =
                    new TH2F(Form("%s_xy", detName), Form("%s xy", detName), 1000, -50., 50, 429, -6.006, 6.006);
                fh_xy_Fib[ifibcount]->GetXaxis()->SetTitle("x / cm ");
                fh_xy_Fib[ifibcount]->GetYaxis()->SetTitle("y / cm");

                fh_xy_Fib_ac[ifibcount] = new TH2F(
                    Form("%s_xy_ac", detName), Form("%s xy after cuts", detName), 1000, -50., 50, 429, -6.006, 6.006);
                fh_xy_Fib_ac[ifibcount]->GetXaxis()->SetTitle("x / cm ");
                fh_xy_Fib_ac[ifibcount]->GetYaxis()->SetTitle("y / cm");
            }
            else
            {
                fh_xy_Fib[ifibcount] =
                    new TH2F(Form("%s_xy", detName), Form("%s xy", detName), 1000, -50., 50., 1000, -50., 50);
                fh_xy_Fib[ifibcount]->GetXaxis()->SetTitle("x / cm ");
                fh_xy_Fib[ifibcount]->GetYaxis()->SetTitle("y / cm");

                fh_xy_Fib_ac[ifibcount] = new TH2F(
                    Form("%s_xy_ac", detName), Form("%s xy after cuts", detName), 1000, -50., 50, 1000, -50., 50);
                fh_xy_Fib_ac[ifibcount]->GetXaxis()->SetTitle("x / cm ");
                fh_xy_Fib_ac[ifibcount]->GetYaxis()->SetTitle("y / cm");
            }

            // Multiplicity (number of hit fibers):
            fh_mult_Fib[ifibcount] = new TH1F(Form("%s_mult", detName), Form("%s # of fibers", detName), 200, 0., 200.);
            fh_mult_Fib[ifibcount]->GetXaxis()->SetTitle("Multiplicity");
            fh_mult_Fib[ifibcount]->GetYaxis()->SetTitle("Counts");

            fh_mult_Fib_ac[ifibcount] =
                new TH1F(Form("%s_mult_ac", detName), Form("%s # of fibers after cuts", detName), 20, 0., 20.);
            fh_mult_Fib_ac[ifibcount]->GetXaxis()->SetTitle("Multiplicity");
            fh_mult_Fib_ac[ifibcount]->GetYaxis()->SetTitle("Counts");

            // ToT MAPMT:
            fh_ToT_Fib[ifibcount] =
                new TH2F(Form("%s_tot_m", detName), Form("%s ToT ", detName), 600, -30., 30, 400, 0., 40.);
            fh_ToT_Fib[ifibcount]->GetXaxis()->SetTitle("Fiber x / cm");
            fh_ToT_Fib[ifibcount]->GetYaxis()->SetTitle("ToT / ns");

            fh_ToT_Fib_ac[ifibcount] =
                new TH2F(Form("%s_tot_ac", detName), Form("%s ToT after cuts", detName), 600, -30., 30, 400, 0., 40);
            fh_ToT_Fib_ac[ifibcount]->GetXaxis()->SetTitle("Fiber x / cm");
            fh_ToT_Fib_ac[ifibcount]->GetYaxis()->SetTitle("ToT / ns");

            // ToF Tofd -> Fiber:
            fh_Fib_ToF[ifibcount] = new TH2F(
                Form("%s_tof", detName), Form("%s ToF Tofd to Fiber", detName), 600, -30., 30, 2000, -1000., 1000.);
            fh_Fib_ToF[ifibcount]->GetYaxis()->SetTitle("ToF / ns");
            fh_Fib_ToF[ifibcount]->GetXaxis()->SetTitle("x / cm");

            fh_Fib_ToF_ac[ifibcount] = new TH2F(Form("%s_tof_ac", detName),
                                                Form("%s ToF Tofd to Fiber after cuts", detName),
                                                1000,
                                                -50.,
                                                50,
                                                2000,
                                                -1000.,
                                                1000.);
            fh_Fib_ToF_ac[ifibcount]->GetYaxis()->SetTitle("ToF / ns");
            fh_Fib_ToF_ac[ifibcount]->GetXaxis()->SetTitle("x / cm");

            // Time:
            fh_Fib_Time[ifibcount] =
                new TH2F(Form("%s_time", detName), Form("%s Time", detName), 1000, -50., 50, 1000, -1000., 1000.);
            fh_Fib_Time[ifibcount]->GetYaxis()->SetTitle("Time / ns");
            fh_Fib_Time[ifibcount]->GetXaxis()->SetTitle("x / cm");

            fh_Fib_Time_ac[ifibcount] = new TH2F(
                Form("%s_time_ac", detName), Form("%s Time after cuts", detName), 1000, -50., 50, 1000, -1000., 1000.);
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
                                                  1000,
                                                  -50.,
                                                  50.);
            fh_Fibs_vs_Tofd[ifibcount]->GetYaxis()->SetTitle("Fiber x / cm");
            fh_Fibs_vs_Tofd[ifibcount]->GetXaxis()->SetTitle("Tofd x / cm");

            fh_Fibs_vs_Tofd_ac[ifibcount] = new TH2F(Form("%s_fib_vs_TofdX_ac", detName),
                                                     Form("%s Fiber # vs. Tofd x-pos after cuts", detName),
                                                     200,
                                                     -100,
                                                     100,
                                                     1000,
                                                     -50.,
                                                     50.);
            fh_Fibs_vs_Tofd_ac[ifibcount]->GetYaxis()->SetTitle("Fiber x / cm");
            fh_Fibs_vs_Tofd_ac[ifibcount]->GetXaxis()->SetTitle("Tofd x / cm");

            fh_qtofd_vs_qFib[ifibcount] = new TH2F(
                Form("%s_fibQ_vs_TofdQ", detName), Form("%s Fiber Q vs. Tofd Q", detName), 200, 0, 20, 100, 0, 10);
            fh_qtofd_vs_qFib[ifibcount]->GetXaxis()->SetTitle("Fiber Q");
            fh_qtofd_vs_qFib[ifibcount]->GetYaxis()->SetTitle("Tofd Q");

            fh_qtofd_vs_qFib_ac[ifibcount] = new TH2F(Form("%s_fibQ_vs_TofdQ_ac", detName),
                                                      Form("%s Fiber Q vs. Tofd Q after cuts", detName),
                                                      200,
                                                      0,
                                                      20,
                                                      100,
                                                      0,
                                                      10);
            fh_qtofd_vs_qFib_ac[ifibcount]->GetXaxis()->SetTitle("Fiber Q");
            fh_qtofd_vs_qFib_ac[ifibcount]->GetYaxis()->SetTitle("Tofd Q");

            fh_Fib_ToF_vs_Qtofd[ifibcount] = new TH2F(
                Form("%s_tof_vs_TofdQ", detName), Form("%s ToF vs. Tofd Q", detName), 200, 0, 20, 2000, -1000., 1000.);
            fh_Fib_ToF_vs_Qtofd[ifibcount]->GetYaxis()->SetTitle("ToF / ns");
            fh_Fib_ToF_vs_Qtofd[ifibcount]->GetXaxis()->SetTitle("Tofd Q");

            // hit fiber vs. fiber position:

        } // end if(Mapped)

    } // end for(ifibcount)
    fh_Fib33_vs_Fib31 = new TH2F("fib33_vs_fib31", "Fiber 33 vs. Fiber 31", 1000, -50, 50, 1000, -50., 50.);
    fh_Fib33_vs_Fib31->GetYaxis()->SetTitle("Fiber33");
    fh_Fib33_vs_Fib31->GetXaxis()->SetTitle("Fiber31");

    fh_Fib31_vs_Fib23a = new TH2F("fib31_vs_fib23a", "Fiber 31 vs. Fiber 23a", 500, -10, 10, 1000, -50., 50.);
    fh_Fib31_vs_Fib23a->GetYaxis()->SetTitle("Fiber31");
    fh_Fib31_vs_Fib23a->GetXaxis()->SetTitle("Fiber23a");

    fh_Fib32_vs_Fib30 = new TH2F("fib32_vs_fib30", "Fiber 32 vs. Fiber 30", 1000, -50, 50, 1000, -50., 50.);
    fh_Fib32_vs_Fib30->GetYaxis()->SetTitle("Fiber32");
    fh_Fib32_vs_Fib30->GetXaxis()->SetTitle("Fiber30");

    fh_tofdleft_vs_Fib23b = new TH2F("tofdleft_vs_Fib23b", "Tofd left 30 vs. Fiber 23b", 500, -10, 10, 1000, -50., 50.);
    fh_tofdleft_vs_Fib23b->GetYaxis()->SetTitle("Tofd left");
    fh_tofdleft_vs_Fib23b->GetXaxis()->SetTitle("Fiber23b");

    fh_Fib30_vs_Fib23a = new TH2F("fib30_vs_fib23a", "Fiber 30 vs. Fiber 23a", 500, -10, 10, 1000, -50., 50.);
    fh_Fib30_vs_Fib23a->GetYaxis()->SetTitle("Fiber30");
    fh_Fib30_vs_Fib23a->GetXaxis()->SetTitle("Fiber23a");

    fh_tofdright_vs_Fib23b = new TH2F("tofdright_vs_fib23b", "Tofd right vs. Fiber 23b", 500, -10, 10, 1000, -50., 50.);
    fh_tofdright_vs_Fib23b->GetYaxis()->SetTitle("Tofd right");
    fh_tofdright_vs_Fib23b->GetXaxis()->SetTitle("Fiber23b");

    // dx between fibers vs x
    fh_Fib33_vs_Fib31_dx = new TH2F("fib33_fib31_dx", "dx of Fiber 33 and Fiber 31", 1000, -50, 50, 1000, -50., 50.);
    fh_Fib33_vs_Fib31_dx->GetYaxis()->SetTitle("xFi33 - xFi31 / cm");
    fh_Fib33_vs_Fib31_dx->GetXaxis()->SetTitle("x Fi31 / cm");

    fh_Fib33_vs_Fib31_dt = new TH2F("fib33_fib31_dt", "dt of Fiber 33 and Fiber 31", 1000, -50, 50, 1000, -50., 50.);
    fh_Fib33_vs_Fib31_dt->GetYaxis()->SetTitle("tFi33 - tFi31 / ns");
    fh_Fib33_vs_Fib31_dt->GetXaxis()->SetTitle("x Fi31 / cm");

    fh_Fib31_vs_Fib23a_dx = new TH2F("fib31_fib23a_dx", "dx of Fiber 31 and Fiber 23a", 1000, -50, 50, 1000, -50., 50.);
    fh_Fib31_vs_Fib23a_dx->GetYaxis()->SetTitle("xFi31 - xFi23a / cm");
    fh_Fib31_vs_Fib23a_dx->GetXaxis()->SetTitle("x Fi23a / cm");

    fh_Fib31_vs_Fib23a_dt =
        new TH2F("fib31_fib23a_dt", "dt of Fiber 31 and Fiber 23a", 1000, -50, 50, 5500, -500., 5000.);
    fh_Fib31_vs_Fib23a_dt->GetYaxis()->SetTitle("tFi31 - tFi23a / ns");
    fh_Fib31_vs_Fib23a_dt->GetXaxis()->SetTitle("x Fi23a / cm");

    fh_Fib32_vs_Fib30_dx = new TH2F("fib32_fib30_dx", "dx of Fiber 32 and Fiber 30", 1000, -50, 50, 1000, -50., 50.);
    fh_Fib32_vs_Fib30_dx->GetYaxis()->SetTitle("xFi32 - xFi30 / cm");
    fh_Fib32_vs_Fib30_dx->GetXaxis()->SetTitle("x Fi30 / cm");

    fh_Fib32_vs_Fib30_dt = new TH2F("fib32_fib30_dt", "dt of Fiber 32 and Fiber 30", 1000, -50, 50, 1000, -50., 50.);
    fh_Fib32_vs_Fib30_dt->GetYaxis()->SetTitle("tFi32 - tFi30 / ns");
    fh_Fib32_vs_Fib30_dt->GetXaxis()->SetTitle("x Fi30 / cm");

    fh_Fib30_vs_Fib23b_dx = new TH2F("fib30_fib23b_dx", "dx of Fiber 30 and Fiber 23b", 1000, -50, 50, 1000, -50., 50.);
    fh_Fib30_vs_Fib23b_dx->GetYaxis()->SetTitle("xFi30 - xFi23b / cm");
    fh_Fib30_vs_Fib23b_dx->GetXaxis()->SetTitle("x Fi23b / cm");

    fh_Fib30_vs_Fib23b_dt =
        new TH2F("fib30_fib23b_dt", "dt of Fiber 30 and Fiber 23b", 1000, -10, 10, 5500, -500., 5000.);
    fh_Fib30_vs_Fib23b_dt->GetYaxis()->SetTitle("tFi30 - tFi23b / ns");
    fh_Fib30_vs_Fib23b_dt->GetXaxis()->SetTitle("x Fi23b / cm");

    fh_Fib31_vs_Fib23b_dx = new TH2F("fib31_fib23b_dx", "dx of Fiber 31 and Fiber 23b", 1000, -50, 50, 1000, -50., 50.);
    fh_Fib31_vs_Fib23b_dx->GetYaxis()->SetTitle("xFi31 - xFi23b / cm");
    fh_Fib31_vs_Fib23b_dx->GetXaxis()->SetTitle("x Fi23b / cm");

    fh_Fib31_vs_Fib23b_dt =
        new TH2F("fib31_fib23b_dt", "dt of Fiber 31 and Fiber 23b", 1000, -10, 10, 5500, -500., 5000.);
    fh_Fib31_vs_Fib23b_dt->GetYaxis()->SetTitle("tFi31 - tFi23b / ns");
    fh_Fib31_vs_Fib23b_dt->GetXaxis()->SetTitle("x Fi23b / cm");

    fh_Fib30_vs_Fib23a_dx = new TH2F("fib30_fib23a_dx", "dx of Fiber 30 and Fiber 23a", 1000, -50, 50, 1000, -50., 50.);
    fh_Fib30_vs_Fib23a_dx->GetYaxis()->SetTitle("xFi30 - xFi23a / cm");
    fh_Fib30_vs_Fib23a_dx->GetXaxis()->SetTitle("x Fi23a / cm");

    fh_Fib30_vs_Fib23a_dt =
        new TH2F("fib30_fib23a_dt", "dt of Fiber 30 and Fiber 23a", 1000, -50, 50, 5500, -500., 5000.);
    fh_Fib30_vs_Fib23a_dt->GetYaxis()->SetTitle("tFi30 - tFi23a / ns");
    fh_Fib30_vs_Fib23a_dt->GetXaxis()->SetTitle("x Fi23a / cm");

    fh_tofd_vs_Fib31_ac = new TH2F("tofd_vs_fib31_ac", "tofd vs. Fiber 31_ac", 1000, -50, 50, 50, -50., 50.);
    fh_tofd_vs_Fib31_ac->GetYaxis()->SetTitle("tofd");
    fh_tofd_vs_Fib31_ac->GetXaxis()->SetTitle("Fiber31");

    fh_tofd_vs_Fib31 = new TH2F("tofd_vs_fib31", "tofd vs. Fiber 31", 150, -30, 30, 50, -50., 50.);
    fh_tofd_vs_Fib31->GetYaxis()->SetTitle("tofd");
    fh_tofd_vs_Fib31->GetXaxis()->SetTitle("Fiber31");

    fh_ytofd_vs_yFib31 = new TH2F("ytofd_vs_yfib31", "ytofd vs. yFiber 31", 150, -30, 30, 50, -50., 50.);
    fh_ytofd_vs_yFib31->GetYaxis()->SetTitle("tofd");
    fh_ytofd_vs_yFib31->GetXaxis()->SetTitle("Fiber31");

    fh_tofd_vs_Fib30_ac = new TH2F("tofd_vs_fib30_ac", "tofd vs. Fiber 30_ac", 150, -30, 30, 50, -50., 50.);
    fh_tofd_vs_Fib30_ac->GetYaxis()->SetTitle("tofd");
    fh_tofd_vs_Fib30_ac->GetXaxis()->SetTitle("Fiber30");

    fh_tofd_vs_Fib30 = new TH2F("tofd_vs_fib30", "tofd vs. Fiber 30", 150, -30, 30, 50, -50., 50.);
    fh_tofd_vs_Fib30->GetYaxis()->SetTitle("tofd");
    fh_tofd_vs_Fib30->GetXaxis()->SetTitle("Fiber30");

    fh_ytofd_vs_yFib30 = new TH2F("ytofd_vs_yfib30", "ytofd vs. yFiber 30", 150, -30, 30, 50, -50., 50.);
    fh_ytofd_vs_yFib30->GetYaxis()->SetTitle("tofd");
    fh_ytofd_vs_yFib30->GetXaxis()->SetTitle("Fiber30");

    fh_Fib33_vs_Fib31_ac = new TH2F("fib33_vs_fib31_ac", "Fiber 33 vs. Fiber 31_ac", 150, -30, 30, 150, -30, 30);
    fh_Fib33_vs_Fib31_ac->GetYaxis()->SetTitle("Fiber33");
    fh_Fib33_vs_Fib31_ac->GetXaxis()->SetTitle("Fiber31");

    fh_Fib31_vs_Fib23a_ac = new TH2F("fib31_vs_fib23a_ac", "Fiber 31 vs. Fiber 23a_ac", 500, -10, 10, 150, -30, 30);
    fh_Fib31_vs_Fib23a_ac->GetYaxis()->SetTitle("Fiber31");
    fh_Fib31_vs_Fib23a_ac->GetXaxis()->SetTitle("Fiber23a");

    fh_Fib32_vs_Fib30_ac = new TH2F("fib32_vs_fib30_ac", "Fiber 32 vs. Fiber 30_ac", 150, -30, 30, 150, -30, 30);
    fh_Fib32_vs_Fib30_ac->GetYaxis()->SetTitle("Fiber32");
    fh_Fib32_vs_Fib30_ac->GetXaxis()->SetTitle("Fiber30");

    fh_tofdleft_vs_Fib23b_ac =
        new TH2F("tofdleft_vs_Fib23b_ac", "Tofd left vs. Fiber 23b_ac", 500, -10, 10, 1000, -50., 50.);
    fh_tofdleft_vs_Fib23b_ac->GetYaxis()->SetTitle("Tofd left");
    fh_tofdleft_vs_Fib23b_ac->GetXaxis()->SetTitle("Fiber23b");

    fh_Fib30_vs_Fib23a_ac = new TH2F("fib30_vs_fib23a_ac", "Fiber 30 vs. Fiber 23a_ac", 500, -10, 10, 150, -30, 30);
    fh_Fib30_vs_Fib23a_ac->GetYaxis()->SetTitle("Fiber30");
    fh_Fib30_vs_Fib23a_ac->GetXaxis()->SetTitle("Fiber23a");

    fh_tofdright_vs_Fib23b_ac =
        new TH2F("tofdright_vs_fib23b_ac", "Tofd right vs. Fiber 23b_ac", 500, -10, 10, 1000, -50., 50.);
    fh_tofdright_vs_Fib23b_ac->GetYaxis()->SetTitle("Tofd right");
    fh_tofdright_vs_Fib23b_ac->GetXaxis()->SetTitle("Fiber23b");

    // dx between fibers vs x
    fh_Fib33_vs_Fib31_dx_ac =
        new TH2F("fib33_fib31_dx_ac", "dx of Fiber 33 and Fiber 31_ac", 1000, -50, 50, 1000, -50., 50.);
    fh_Fib33_vs_Fib31_dx_ac->GetYaxis()->SetTitle("xFi33 - xFi31 / cm");
    fh_Fib33_vs_Fib31_dx_ac->GetXaxis()->SetTitle("x Fi31 / cm");

    fh_Fib31_vs_Fib23a_dx_ac =
        new TH2F("fib31_fib23a_dx_ac", "dx of Fiber 31 and Fiber 23a_ac", 1000, -50, 50, 1000, -50., 50.);
    fh_Fib31_vs_Fib23a_dx_ac->GetYaxis()->SetTitle("xFi31 - xFi23a / cm");
    fh_Fib31_vs_Fib23a_dx_ac->GetXaxis()->SetTitle("x Fi23a / cm");

    fh_Fib32_vs_Fib30_dx_ac =
        new TH2F("fib32_fib30_dx_ac", "dx of Fiber 32 and Fiber 30_ac", 1000, -50, 50, 1000, -50., 50.);
    fh_Fib32_vs_Fib30_dx_ac->GetYaxis()->SetTitle("xFi32 - xFi30 / cm");
    fh_Fib32_vs_Fib30_dx_ac->GetXaxis()->SetTitle("x Fi30 / cm");

    fh_Fib30_vs_Fib23b_dx_ac =
        new TH2F("fib32_fib23b_dx_ac", "dx of Fiber 32 and Fiber 23b_ac", 1000, -50, 50, 1000, -50., 50.);
    fh_Fib30_vs_Fib23b_dx_ac->GetYaxis()->SetTitle("xFi32 - xFi23b / cm");
    fh_Fib30_vs_Fib23b_dx_ac->GetXaxis()->SetTitle("x Fi23b / cm");

    fh_Fib31_vs_Fib23b_dx_ac =
        new TH2F("fib31_fib23b_dx_ac", "dx of Fiber 31 and Fiber 23b_ac", 1000, -50, 50, 1000, -50., 50.);
    fh_Fib31_vs_Fib23b_dx_ac->GetYaxis()->SetTitle("xFi31 - xFi23b / cm");
    fh_Fib31_vs_Fib23b_dx_ac->GetXaxis()->SetTitle("x Fi23b / cm");

    fh_Fib30_vs_Fib23a_dx_ac =
        new TH2F("fib30_fib23a_dx_ac", "dx of Fiber 30 and Fiber 23a_ac", 1000, -50, 50, 1000, -50., 50.);
    fh_Fib30_vs_Fib23a_dx_ac->GetYaxis()->SetTitle("xFi30 - xFi23a / cm");
    fh_Fib30_vs_Fib23a_dx_ac->GetXaxis()->SetTitle("x Fi23a / cm");

    fh_Fib33_vs_Fib31_tot = new TH2F("fib33_fib31_tot", "Fiber 33 vs Fiber 31 tot", 600, 0, 30, 600, 0, 30);
    fh_Fib33_vs_Fib31_tot->GetXaxis()->SetTitle("Q fi33");
    fh_Fib33_vs_Fib31_tot->GetYaxis()->SetTitle("Q fi31");

    fh_Fib32_vs_Fib30_tot = new TH2F("fib32_fib30_tot", "Fiber 32 vs Fiber 30 tot", 600, 0, 30, 600, 0, 30);
    fh_Fib32_vs_Fib30_tot->GetXaxis()->SetTitle("Q fi32");
    fh_Fib32_vs_Fib30_tot->GetYaxis()->SetTitle("Q fi30");

    fh_Fib23a_vs_Fib3X_tot = new TH2F("fib23a_fib3X_tot", "Fiber 3X vs Fiber 23a tot", 600, 0, 30, 600, 0, 30);
    fh_Fib23a_vs_Fib3X_tot->GetXaxis()->SetTitle("Q fi23a");
    fh_Fib23a_vs_Fib3X_tot->GetYaxis()->SetTitle("Q fi3X");

    fh_Fib23b_vs_Fib3X_tot = new TH2F("fib23b_fib3X_tot", "Fiber 3X vs Fiber 23b tot", 600, 0, 30, 600, 0, 30);
    fh_Fib23b_vs_Fib3X_tot->GetXaxis()->SetTitle("Q fi23b");
    fh_Fib23b_vs_Fib3X_tot->GetYaxis()->SetTitle("Q fi3X");

    fh_yTofd_vs_yTofd = new TH2F("yTofd_vs_yTofd", "yTofd_vs_yTofd", 80, -60., 60., 80, -60., 60.);
    fh_yTofd_vs_yTofd->GetYaxis()->SetTitle("Tofd 2 y / cm");
    fh_yTofd_vs_yTofd->GetXaxis()->SetTitle("Tofd 1 y / cm");

    fh_qtofi_vs_qtofd = new TH2F("qtofi_vs_qtofd", "fh_qtofi_vs_qtofd", 200, 0, 10, 200, 0, 10);
    fh_qtofi_vs_qtofd->GetYaxis()->SetTitle("Q tofi");
    fh_qtofi_vs_qtofd->GetXaxis()->SetTitle("Q tofd");

    fh_xy_target = new TH2F("fh_xy_target", "fh_xy_target", 420, -6, 6, 420, -6, 6);
    fh_xy_target->GetYaxis()->SetTitle(" y / cm");
    fh_xy_target->GetXaxis()->SetTitle("x / cm");

    fh_Fi23bY_TofdY_He = new TH2F("Fi23bY_TofdY_He", "Fi23bY_TofdY_He", 80, -60., 60., 420, -6, 6);
    fh_Fi23bY_TofdY_He->GetYaxis()->SetTitle("Fiber y / cm");
    fh_Fi23bY_TofdY_He->GetXaxis()->SetTitle("Tofd y / cm");

    fh_Fi23bY_TofdY_C = new TH2F("Fi23bY_TofdY_C", "fh_Fi23bY_TofdY_C", 80, -60., 60., 420, -6, 6);
    fh_Fi23bY_TofdY_C->GetYaxis()->SetTitle("Fiber y / cm");
    fh_Fi23bY_TofdY_C->GetXaxis()->SetTitle("Tofd y / cm");

    fh_Fi23aX_TofdX_bw1 = new TH2F("Fi23aX_TofdX_bw1", "fh_Fi23aX_TofdX_Ebene1", 80, -60., 60., 420, -6, 6);
    fh_Fi23aX_TofdX_bw1->GetYaxis()->SetTitle("Fiber x / cm");
    fh_Fi23aX_TofdX_bw1->GetXaxis()->SetTitle("Tofd x / cm");

    fh_Fi23aX_TofdX_bw2 = new TH2F("Fi23aX_TofdX_bw2", "fh_Fi23aX_TofdX_Ebene2", 80, -60., 60., 420, -6, 6);
    fh_Fi23aX_TofdX_bw2->GetYaxis()->SetTitle("Fiber x / cm");
    fh_Fi23aX_TofdX_bw2->GetXaxis()->SetTitle("Tofd x / cm");

    fh_Fi23aX_TofdX_He = new TH2F("Fi23aX_TofdX_He", "fh_Fi23aX_TofdX_He", 80, -60., 60., 420, -6, 6);
    fh_Fi23aX_TofdX_He->GetYaxis()->SetTitle("Fiber x / cm");
    fh_Fi23aX_TofdX_He->GetXaxis()->SetTitle("Tofd x / cm");

    fh_Fi23aX_TofdX_C = new TH2F("Fi23aX_TofdX_C", "fh_Fi23aX_TofdX_C", 80, -60., 60., 420, -6, 6);
    fh_Fi23aX_TofdX_C->GetYaxis()->SetTitle("Fiber x / cm");
    fh_Fi23aX_TofdX_C->GetXaxis()->SetTitle("Tofd x / cm");

    fh_yTofdHe_vs_yTofdC = new TH2F("fh_yTofdHe_vs_yTofdC", "fh_yTofdHe_vs_yTofdC", 80, -60., 60., 80, -60, 60);
    fh_yTofdHe_vs_yTofdC->GetXaxis()->SetTitle("He y / cm");
    fh_yTofdHe_vs_yTofdC->GetYaxis()->SetTitle("C y / cm");

    fh_xTofdHe_vs_xTofdC = new TH2F("fh_xTofdHe_vs_xTofdC", "fh_xTofdHe_vs_xTofdC", 80, -60., 60., 80, -60, 60);
    fh_xTofdHe_vs_xTofdC->GetXaxis()->SetTitle("He x / cm");
    fh_xTofdHe_vs_xTofdC->GetYaxis()->SetTitle("C x / cm");

    fh_Fi301Y_TofdY_He = new TH2F("Fi301Y_TofdY_He", "Fi301Y_TofdY_He", 80, -60., 60., 1200, -60, 60);
    fh_Fi301Y_TofdY_He->GetYaxis()->SetTitle("Fiber y / cm");
    fh_Fi301Y_TofdY_He->GetXaxis()->SetTitle("Tofd y / cm");

    fh_Fi301Y_TofdY_C = new TH2F("Fi301Y_TofdY_C", "fh_Fi301Y_TofdY_C", 80, -60., 60., 1200, -60, 60);
    fh_Fi301Y_TofdY_C->GetYaxis()->SetTitle("Fiber y / cm");
    fh_Fi301Y_TofdY_C->GetXaxis()->SetTitle("Tofd y / cm");

    fh_Fi301X_TofdX_He = new TH2F("Fi301X_TofdX_He", "fh_Fi301X_TofdX_He", 80, -60., 60., 1200, -60, 60);
    fh_Fi301X_TofdX_He->GetYaxis()->SetTitle("Fiber x / cm");
    fh_Fi301X_TofdX_He->GetXaxis()->SetTitle("Tofd x / cm");

    fh_Fi301X_TofdX_C = new TH2F("Fi301X_TofdX_C", "fh_Fi301X_TofdX_C", 80, -60., 60., 1200, -60, 60);
    fh_Fi301X_TofdX_C->GetYaxis()->SetTitle("Fiber x / cm");
    fh_Fi301X_TofdX_C->GetXaxis()->SetTitle("Tofd x / cm");

    fh_Fi323Y_TofdY_He = new TH2F("Fi323Y_TofdY_He", "Fi323Y_TofdY_He", 80, -60., 60., 1200, -60, 60);
    fh_Fi323Y_TofdY_He->GetYaxis()->SetTitle("Fiber y / cm");
    fh_Fi323Y_TofdY_He->GetXaxis()->SetTitle("Tofd y / cm");

    fh_Fi323Y_TofdY_C = new TH2F("Fi323Y_TofdY_C", "fh_Fi323Y_TofdY_C", 80, -60., 60., 1200, -60, 60);
    fh_Fi323Y_TofdY_C->GetYaxis()->SetTitle("Fiber y / cm");
    fh_Fi323Y_TofdY_C->GetXaxis()->SetTitle("Tofd y / cm");

    fh_Fi323X_TofdX_He = new TH2F("Fi323X_TofdX_He", "fh_Fi323X_TofdX_He", 80, -60., 60., 1200, -60, 60);
    fh_Fi323X_TofdX_He->GetYaxis()->SetTitle("Fiber x / cm");
    fh_Fi323X_TofdX_He->GetXaxis()->SetTitle("Tofd x / cm");

    fh_Fi323X_TofdX_C = new TH2F("Fi323X_TofdX_C", "fh_Fi323X_TofdX_C", 80, -60., 60., 1200, -60, 60);
    fh_Fi323X_TofdX_C->GetYaxis()->SetTitle("Fiber x / cm");
    fh_Fi323X_TofdX_C->GetXaxis()->SetTitle("Tofd x / cm");

    for (Int_t i = 0; i < 21; i++)
    {
        fh_Fi23aToF_Q[i] = new TH2F(
            Form("Fi23a_ToF_vs_Q_mult%d", i), Form("Fi23a ToF vs_Q for mult %d", i), 20, 0., 20, 1000, -1000., 1000.);
        fh_Fi23bToF_Q[i] = new TH2F(
            Form("Fi23b_ToF_vs_Q_mult%d", i), Form("Fi23b ToF vs_Q for mult %d", i), 20, 0., 20, 1000, -1000., 1000.);
    }

    fh_TimePreviousEventFi23a =
        new TH2F("TimePreviousEvent23a ", "Time between 2 particles 23a", 3000, -1500, 1500, 3000, -1500, 1500);
    fh_TimePreviousEventFi23a->GetYaxis()->SetTitle("dtime fi23a/ ns");
    fh_TimePreviousEventFi23a->GetXaxis()->SetTitle("dtime tofd / ns");

    fh_TimePreviousEventFi23b =
        new TH2F("TimePreviousEvent23b", "Time between 2 particles 23b", 3000, -1500, 1500, 3000, -1500, 1500);
    fh_TimePreviousEventFi23b->GetYaxis()->SetTitle("dtime fi23b/ ns");
    fh_TimePreviousEventFi23b->GetXaxis()->SetTitle("dtime tofd / ns");

    //---------------------------------------------------------------------------------------------------
    // TofD detector

    if (fHitItems.at(DET_TOFD))
    {

        // xy:
        fh_xy_tofd = new TH2F("tofd_xy", "tofd xy", 240, -60., 60., 200, -100., 100.);
        fh_xy_tofd->GetXaxis()->SetTitle("x / cm ");
        fh_xy_tofd->GetYaxis()->SetTitle("y / cm");

        fh_xy_tofd_ac = new TH2F("tofd_xy_ac", "tofd xy after cuts", 240, -60., 60., 200, -100., 100.);
        fh_xy_tofd_ac->GetXaxis()->SetTitle("x / cm ");
        fh_xy_tofd_ac->GetYaxis()->SetTitle("y / cm");

        fh_tofd_charge = new TH2F("tofd_Q", "Charge of Tofd vs yToFD", 200, -100., 100., 100, 0., 10.);
        fh_tofd_charge->GetXaxis()->SetTitle("y / cm ");
        fh_tofd_charge->GetYaxis()->SetTitle("Q");

        fh_tofd_charge_ac = new TH2F("tofd_Q_ac", "Charge of Tofd vs yToFD after cuts", 200, -100., 100., 100, 0., 10.);
        fh_tofd_charge_ac->GetXaxis()->SetTitle("y / cm ");
        fh_tofd_charge_ac->GetYaxis()->SetTitle("Q");

        fh_tofd_mult = new TH1F("tofd_mult", "ToFD multiplicits ", 100, 0, 100);
        fh_tofd_mult->GetXaxis()->SetTitle("multiplicity");
        fh_tofd_mult->GetYaxis()->SetTitle("counts");

        fh_tofd_mult_ac = new TH1F("tofd_mult_ac", "ToFD multiplicits after cuts", 100, 0, 100);
        fh_tofd_mult_ac->GetXaxis()->SetTitle("multiplicity");
        fh_tofd_mult_ac->GetYaxis()->SetTitle("counts");

        fh_TimePreviousEvent =
            new TH2F("TimePreviousEvent", "Time between 2 particles vs q1+q2 ", 3000, -1500, 1500, 200, 0, 20);
        fh_TimePreviousEvent->GetXaxis()->SetTitle("t2-t1 / ns");
        fh_TimePreviousEvent->GetYaxis()->SetTitle("q1+q2");

        fh_tofd_time = new TH1F("tofd_time", "Tofd times ", 40000, -2000, 2000);
        fh_tofd_time->GetXaxis()->SetTitle("time / ns");
        fh_tofd_time->GetYaxis()->SetTitle("counts");

        fh_tofd_time_ac = new TH1F("tofd_time_ac", "Tofd times after cut", 40000, -2000, 2000);
        fh_tofd_time_ac->GetXaxis()->SetTitle("time / ns");
        fh_tofd_time_ac->GetYaxis()->SetTitle("counts");

        fh_tofd_q2_vs_q1 = new TH2F("tofd_q2_vs_q1", "tofd q2 vs. q1", 100, 0., 10., 100, 0., 10.);
        fh_tofd_q2_vs_q1->GetXaxis()->SetTitle("q1");
        fh_tofd_q2_vs_q1->GetYaxis()->SetTitle("q2");

        fh_tofd_q2_vs_q1_ac = new TH2F("tofd_q2_vs_q1_ac", "tofd q2 vs. q1 after cut", 100, 0., 10., 100, 0., 10.);
        fh_tofd_q2_vs_q1_ac->GetXaxis()->SetTitle("q1");
        fh_tofd_q2_vs_q1_ac->GetYaxis()->SetTitle("q2");

        fh_qsum_mult_fi23a = new TH2F("qsum_mult_fi23a", "qsum_mult_fi23a", 100, 0., 100., 40, 0, 40);
        fh_qsum_mult_fi23a->GetXaxis()->SetTitle("ncouts[23a]");
        fh_qsum_mult_fi23a->GetYaxis()->SetTitle("Sum of Q in Tofd");

        fh_qsum_mult_fi23b = new TH2F("qsum_mult_fi23b", "qsum_mult_fi23b", 100, 0., 100., 40, 0, 40);
        fh_qsum_mult_fi23b->GetXaxis()->SetTitle("ncouts[23b]");
        fh_qsum_mult_fi23b->GetYaxis()->SetTitle("Sum of Q in Tofd");

        for (int i = 0; i < 3; i++)
        {
            const char* ZZ;
            if (i == 0)
                ZZ = Form("%01d", 2);
            if (i == 1)
                ZZ = Form("%01d", 6);
            if (i == 2)
                ZZ = Form("%01d", 8);
            fh_tofd_x_vs_y_z[i] = new TH2F(
                Form("tofd_plane1_x_vs_y_Z%s", ZZ), Form("tofd plane1 BarId vs_y_Z=%s", ZZ), 45, 0, 45, 200, -100, 100);
            fh_tofd_x_vs_y_z[i]->GetXaxis()->SetTitle("BarId");
            fh_tofd_x_vs_y_z[i]->GetYaxis()->SetTitle("y / cm");
        }
        for (int i = 3; i < 6; i++)
        {
            const char* ZZ;
            if (i == 3)
                ZZ = Form("%01d", 2);
            if (i == 4)
                ZZ = Form("%01d", 6);
            if (i == 5)
                ZZ = Form("%01d", 8);
            fh_tofd_x_vs_y_z[i] = new TH2F(
                Form("tofd_plane2_x_vs_y_Z%s", ZZ), Form("tofd plane2 BarId vs_y_Z=%s", ZZ), 45, 0, 45, 200, -100, 100);
            fh_tofd_x_vs_y_z[i]->GetXaxis()->SetTitle("BarId");
            fh_tofd_x_vs_y_z[i]->GetYaxis()->SetTitle("y / cm");
        }
    }

    if (fHitItems.at(DET_CALIFA))
    {
        fh_califa_energy = new TH2F("fh_califa_energy", "Califa E vs Nb crystals", 50, 0, 50, 1000, 0., 10000.);
        fh_califa_energy->GetYaxis()->SetTitle("Energy / keV");
        fh_califa_energy->GetXaxis()->SetTitle("Num crystals in clustre");

        fh_califa_energy_dc =
            new TH2F("fh_califa_energy_dc", "Califa E  Doppler-corrected vs Nb crystals", 50, 0, 50, 1000, 0., 10000.);
        fh_califa_energy_dc->GetYaxis()->SetTitle("Energy / keV");
        fh_califa_energy_dc->GetXaxis()->SetTitle("Num crystals in clustre");

        fh_califa_time = new TH2F("fh_califa_time", "Califa time vs Nb crystals", 50, 0, 50, 1500, 1000., 4000.);
        fh_califa_time->GetYaxis()->SetTitle("time / ns");
        fh_califa_time->GetXaxis()->SetTitle("Num crystals in clustre");

        fh_califa_tofd_barrel =
            new TH2F("fh_califa_tofd_barrel", "Califa time vs tofd time", 1500, 1000., 4000., 200, -200, 200);
        fh_califa_tofd_barrel->GetYaxis()->SetTitle("tofd time / ns");
        fh_califa_tofd_barrel->GetXaxis()->SetTitle("califa tome / ns");

        fh_califa_tofd_iphos =
            new TH2F("fh_califa_tofd_iphos", "Califa time vs tofd time", 1500, 1000., 4000., 200, -200, 200);
        fh_califa_tofd_iphos->GetYaxis()->SetTitle("tofd time / ns");
        fh_califa_tofd_iphos->GetXaxis()->SetTitle("califa tome / ns");

        fh_califa_energy_select_barrel = new TH2F(
            "fh_califa_energy_select_barrel", "Califa E Doppler-corrected after cuts", 50, 0, 50, 1000, 0., 10000.);
        fh_califa_energy_select_barrel->GetYaxis()->SetTitle("Energy / MeV");
        fh_califa_energy_select_barrel->GetXaxis()->SetTitle("Num crystals in clustre");

        fh_califa_energy_select_iphos = new TH2F(
            "fh_califa_energy_select_iphos", "Califa E Doppler-corrected after cuts", 50, 0, 50, 1000, 0., 10000.);
        fh_califa_energy_select_iphos->GetYaxis()->SetTitle("Energy / MeV");
        fh_califa_energy_select_iphos->GetXaxis()->SetTitle("Num crystals in clustre");

        fh_califa_cluster =
            new TH2F("fh_califa_cluster", "Califa MotherId vs Nb crystals rolu=false", 50, 0, 50, 1600, 900., 2500.);
        fh_califa_cluster->GetYaxis()->SetTitle("MotherId");
        fh_califa_cluster->GetXaxis()->SetTitle("Num crystals in clustre");

        fh_califa_overflow_antirolu = new TH2F(
            "fh_califa_overflow_antirolu", "Califa overflow is_nan(en) rolu=false", 5, 0, 5, 1600, 900., 2500.);
        fh_califa_overflow_antirolu->GetXaxis()->SetTitle("ClusterType");
        fh_califa_overflow_antirolu->GetYaxis()->SetTitle("motherId");

        fh_califa_overflow_rolu =
            new TH2F("fh_califa_overflow_rolu", "Califa overflow is_nan(en) rolu=true", 5, 0, 5, 1600, 900., 2500.);
        fh_califa_overflow_rolu->GetXaxis()->SetTitle("ClusterType");
        fh_califa_overflow_rolu->GetYaxis()->SetTitle("motherId");

        fh_califa_cluster_rolu = new TH2F(
            "fh_califa_cluster_rolu", "Califa MotherId vs Nb crystals rolu=true", 50, 0, 50, 1600, 900., 2500.);
        fh_califa_cluster_rolu->GetYaxis()->SetTitle("MotherId");
        fh_califa_cluster_rolu->GetXaxis()->SetTitle("Num crystals in clustre");

        fh_califa_overflow =
            new TH2F("fh_califa_overflow", "Califa overflow is_nan(en) after cuts", 5, 0, 5, 1600, 900., 2500.);
        fh_califa_overflow->GetXaxis()->SetTitle("ClusterType");
        fh_califa_overflow->GetYaxis()->SetTitle("motherId");

        fh_califa_good_rolu =
            new TH2F("fh_califa_good_rolu", "Califa good !(is:nan(en)) rolu=true", 5, 0, 5, 1600, 900., 2500.);
        fh_califa_good_rolu->GetXaxis()->SetTitle("ClusterType");
        fh_califa_good_rolu->GetYaxis()->SetTitle("motherId");

        fh_califa_good_antirolu =
            new TH2F("fh_califa_good_antirolu", "Califa good !(is:nan(en)) rolu=false", 5, 0, 5, 1600, 900., 2500.);
        fh_califa_good_antirolu->GetXaxis()->SetTitle("ClusterType");
        fh_califa_good_antirolu->GetYaxis()->SetTitle("motherId");

        fh_califa_good = new TH2F("fh_califa_good", "Califa good !(is:nan(en)) after cuts", 5, 0, 5, 1600, 900., 2500.);
        fh_califa_good->GetXaxis()->SetTitle("ClusterType");
        fh_califa_good->GetYaxis()->SetTitle("motherId");

        fh_califa_motherid_clustertype =
            new TH2F("fh_califa_motherid_clustertype", "ClusterType vs MotherId", 1600, 900., 2500., 5, 0, 5);
        fh_califa_motherid_clustertype->GetXaxis()->SetTitle("MotherID");
        fh_califa_motherid_clustertype->GetYaxis()->SetTitle("ClusterType");
    }

    // Output spectra:
    for (int i = 0; i < 2; i++)
    {
        fh_check_QvsX[i] = new TH2F(Form("fhCheckQvsX%d", i), Form("Check Q vs X det%d", i), 1000, -5, 5, 20, 0, 10);
        fh_check_TvsX[i] =
            new TH2F(Form("fhCheckTvsX%d", i), Form("Check T vs X det%d", i), 500, -5, 5, 500, -100, 100);
        fh_check_XvsY[i] = new TH2F(Form("fhCheckXvsY%d", i), Form("Check X vs Y det%d", i), 500, -5, 5, 500, -5, 5);
    }
    for (int i = 2; i < 6; i++)
    {
        fh_check_QvsX[i] = new TH2F(Form("fhCheckQvsX%d", i), Form("Check Q vs X det%d", i), 300, -30, 30, 20, 0, 10);
        fh_check_TvsX[i] =
            new TH2F(Form("fhCheckTvsX%d", i), Form("Check T vs X det%d", i), 300, -30, 30, 500, -100, 100);
        fh_check_XvsY[i] =
            new TH2F(Form("fhCheckXvsY%d", i), Form("Check X vs Y det%d", i), 300, -30, 30, 300, -30, 30);
    }
    fh_check_QvsX[6] = new TH2F(Form("fhCheckQvsX%d", 6), Form("Check Q vs X det%d", 6), 100, -135., 135., 20, 0, 10);
    fh_check_TvsX[6] =
        new TH2F(Form("fhCheckTvsX%d", 6), Form("Check T vs X det%d", 6), 240, -60., 60., 500, -100, 100);
    fh_check_XvsY[6] =
        new TH2F(Form("fhCheckXvsY%d", 6), Form("Check X vs Y det%d", 6), 240, -60., 60., 200, -100., 100.);

    /*
           TCanvas* check = new TCanvas("CheckingMom", "CheckingMom", 10, 10, 1100, 700);
           check->Divide(7,3);
           for(int i = 0; i < 7; i++)
           {
                check->cd(i+1);
                fh_check_QvsX[i]->Draw("colz");
            }
           for(int i = 0; i < 7; i++)
           {
                check->cd(i+8);
                fh_check_TvsX[i]->Draw("colz");
            }
           for(int i = 0; i < 7; i++)
           {
                check->cd(i+15);
                fh_check_XvsY[i]->Draw("colz");
            }

    */

    return kSUCCESS;
}

void R3BPreTrackS494::Exec(Option_t* option)
{
    /*
        if (fNEvents / 10000. == (int)fNEvents / 10000)
            std::cout << "\rEvents: " << fNEvents << " / " << maxevent << " (" << (int)(fNEvents * 100. / maxevent)
                      << " %) " << std::flush;
    */
    fNEvents += 1;

    FairRootManager* mgr = FairRootManager::Instance();
    if (NULL == mgr)
        LOG(fatal) << "FairRootManager not found";

    if (header && !fSimu)
    {
        time = header->GetTimeStamp();

        // if (time > 0) cout << "header time: " << time << endl;
        if (time_start == 0 && time > 0)
        {
            time_start = time;
            fNEvents_start = fNEvents;
            // cout << "Start event number " << fNEvents_start << endl;
        }

        if (header->GetTrigger() == 12)
        {
            cout << "spill start" << endl;
            num_spills++;
            // cout << "Spill: " << num_spills << endl;
        }
        if (header->GetTrigger() == 13)
        {
            //  cout << "spill stop" << endl;
        }

        fh_Trigger->Fill(header->GetTrigger());
        //   check for requested trigger (Todo: should be done globablly / somewhere else)
        if ((fTrigger >= 0) && (header) && (header->GetTrigger() != fTrigger) && !fSimu)
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
        Int_t fTpat_bit1 = fTpat1 - 1;
        Int_t fTpat_bit2 = fTpat2 - 1;
        for (int i = 0; i < 16; i++)
        {
            tpatbin = (header->GetTpat() & (1 << i));
            if (tpatbin != 0 && (i < fTpat_bit1 || i > fTpat_bit2) && !fSimu)
            {
                counterWrongTpat++;
                return;
            }
        }
    }

    if (fMappedItems.at(DET_BMON) && !fSimu)
    {
        unsigned long IC = 0;
        unsigned long SEETRAM_raw = 0;
        Double_t SEETRAM = 0;
        unsigned long TOFDOR = 0;

        auto detBmon = fMappedItems.at(DET_BMON);
        Int_t nHitsbm = detBmon->GetEntriesFast();
        // cout<<"Bmon hits: "<<nHitsbm<<endl;

        for (Int_t ihit = 0; ihit < nHitsbm; ihit++)
        {
            R3BBeamMonitorMappedData* hit = (R3BBeamMonitorMappedData*)detBmon->At(ihit);
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
    if (fMappedItems.at(DET_ROLU) && !fSimu)
    {
        // rolu
        auto detRolu = fMappedItems.at(DET_ROLU);
        Int_t nHitsRolu = detRolu->GetEntriesFast();
        // cout<<"ROLU hits: "<<nHitsRolu<<endl;

        for (Int_t ihit = 0; ihit < nHitsRolu; ihit++)
        {
            R3BRoluMappedData* hitRolu = (R3BRoluMappedData*)detRolu->At(ihit);
            if (!hitRolu)
                continue;

            // channel numbers are stored 1-based (1..n)
            Int_t iDet = hitRolu->GetDetector(); // 1..
            Int_t iCha = hitRolu->GetChannel();  // 1..
            RoluCut = true;
        }
    }

    if (RoluCut && !fSimu)
    {

        if (fHitItems.at(DET_CALIFA))
        {
            // CALIFA
            auto detCalifa = fHitItems.at(DET_CALIFA);
            Int_t nHitsCalifa = detCalifa->GetEntriesFast();
            // cout<<"Califa hits: "<<nHitsCalifa<<endl;
            for (Int_t ihit = 0; ihit < nHitsCalifa; ihit++)
            {
                R3BCalifaClusterData* hitCalifa = (R3BCalifaClusterData*)detCalifa->At(ihit);
                if (!hitCalifa)
                    continue;

                if (IS_NAN(hitCalifa->GetEnergy()))
                    fh_califa_overflow_rolu->Fill(hitCalifa->GetClusterType(), hitCalifa->GetMotherCrystal());
                else
                    fh_califa_good_rolu->Fill(hitCalifa->GetClusterType(), hitCalifa->GetMotherCrystal());

                if (hitCalifa->GetClusterType() == 1) // gammas
                {
                    ULong64_t timeCalifa = hitCalifa->GetTime();
                    Double_t timerelCalifa = (double)(timeCalifa - time);
                    Int_t CrystalNb = hitCalifa->GetCrystalList().size();
                    Double_t Energy = hitCalifa->GetEnergy();

                    if (!(IS_NAN(Energy))) //
                    {
                        fh_califa_cluster_rolu->Fill(hitCalifa->GetCrystalList().size(), hitCalifa->GetMotherCrystal());
                    }
                }
            }
        }

        counterRolu++;
        return;
    }

    Bool_t CalifaHit = false;
    Int_t nHitsCalifa = 0;
    vector<double> eCalifa;
    vector<double> tCalifa;
    Double_t E_0 = 15.99065084 * 931.494028;
    Double_t E_Beam = E_0 + 16. * 500.;
    Double_t BETA = sqrt(1. - (E_0 * E_0) / (E_Beam * E_Beam));
    Double_t GAMMA = 1. / (sqrt(1. - BETA * BETA));
    if (fHitItems.at(DET_CALIFA))
    {
        CalifaHit = true;

        // CALIFA
        auto detCalifa = fHitItems.at(DET_CALIFA);
        nHitsCalifa = detCalifa->GetEntriesFast();
        // cout<<"Califa hits: "<<nHitsCalifa<<endl;
        for (Int_t ihit = 0; ihit < nHitsCalifa; ihit++)
        {
            R3BCalifaClusterData* hitCalifa = (R3BCalifaClusterData*)detCalifa->At(ihit);
            if (!hitCalifa)
                continue;

            ULong64_t timeCalifa = hitCalifa->GetTime();
            Double_t timerelCalifa = (double)(timeCalifa - time);
            // cout<<"califa time: "<<time <<", "<< timeCalifa<<", "<<diff<<endl;
            Int_t CrystalNb = hitCalifa->GetCrystalList().size();
            Double_t Energy = hitCalifa->GetEnergy();
            Double_t theta = hitCalifa->GetTheta();
            Double_t Energy_dc = Energy * GAMMA * (1 - BETA * TMath::Cos(theta));

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
            */
            if (IS_NAN(Energy))
                fh_califa_overflow_antirolu->Fill(hitCalifa->GetClusterType(), hitCalifa->GetMotherCrystal());
            else
                fh_califa_good_antirolu->Fill(hitCalifa->GetClusterType(), hitCalifa->GetMotherCrystal());

            fh_califa_motherid_clustertype->Fill(hitCalifa->GetMotherCrystal(), hitCalifa->GetClusterType());

            if (hitCalifa->GetClusterType() == 1 && !(IS_NAN(Energy))) // gammas
            {
                fh_califa_cluster->Fill(hitCalifa->GetCrystalList().size(), hitCalifa->GetMotherCrystal());

                if (hitCalifa->GetMotherCrystal() > 927 &&
                    hitCalifa->GetMotherCrystal() < 1953) // barrel gamma-channels
                {
                    fh_califa_energy->Fill(hitCalifa->GetCrystalList().size(), Energy);
                    fh_califa_energy_dc->Fill(hitCalifa->GetCrystalList().size(), Energy_dc);
                }
            }
        }
        if (CalifaHit)
            counterCalifa++;
    }

    if (fMCTrack) // for simulated data
    {
        // read in Monte Carlo Track parameter

        Int_t nHitsMCTrack = fMCTrack->GetEntriesFast();
        // cout << "MCTrack hits: " << nHitsMCTrack << endl;

        for (Int_t l = 0; l < nHitsMCTrack; l++)
        {
            R3BMCTrack* aTrack = (R3BMCTrack*)fMCTrack->At(l);

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

    Bool_t SelectedPairs = false;
    Double_t xTest = 0.;
    Double_t yTest = 0.;

    Int_t max = 500;
    Int_t detector[max];
    Double_t xdet[max];
    Double_t ydet[max];
    Double_t zdet[max];
    Double_t zdet_s[max];
    Double_t tdet[max];
    Int_t fibdet[max];
    Int_t fibdet_s[max];
    Int_t qdet[max];
    Double_t qdet_s[max];
    Double_t xdet_s[max];
    Double_t ydet_s[max];
    Double_t tdet_s[max];

    Double_t qTofd1ls[max];
    Double_t xTofd1ls[max];
    Double_t yTofd1ls[max];
    Double_t tTofd1ls[max];

    Double_t qTofd1rs[max];
    Double_t xTofd1rs[max];
    Double_t yTofd1rs[max];
    Double_t tTofd1rs[max];

    Double_t qTofd2ls[max];
    Double_t xTofd2ls[max];
    Double_t yTofd2ls[max];
    Double_t tTofd2ls[max];

    Double_t qTofd2rs[max];
    Double_t xTofd2rs[max];
    Double_t yTofd2rs[max];
    Double_t tTofd2rs[max];

    Double_t qTofd3ls[max];
    Double_t xTofd3ls[max];
    Double_t yTofd3ls[max];
    Double_t tTofd3ls[max];

    Double_t qTofd3rs[max];
    Double_t xTofd3rs[max];
    Double_t yTofd3rs[max];
    Double_t tTofd3rs[max];

    Int_t detector_s[max];
    Double_t xFi33[max];
    Double_t yFi33[max];
    Double_t qFi33[max];
    Double_t tFi33[max];
    Double_t timeFi33[max];
    Bool_t fFi33[max];
    Double_t xFi32[max];
    Double_t yFi32[max];
    Double_t qFi32[max];
    Double_t tFi32[max];
    Double_t timeFi32[max];
    Bool_t fFi32[max];
    Double_t xFi31[max];
    Double_t yFi31[max];
    Double_t qFi31[max];
    Double_t tFi31[max];
    Double_t timeFi31[max];
    Bool_t fFi31[max];
    Double_t xFi30[max];
    Double_t yFi30[max];
    Double_t qFi30[max];
    Double_t tFi30[max];
    Double_t timeFi30[max];
    Bool_t fFi30[max];
    Double_t xFi23a[max];
    Double_t yFi23a[max];
    Double_t qFi23a[max];
    Double_t tFi23a[max];
    Double_t timeFi23a[max];
    Bool_t fFi23a[max];
    Double_t xFi23b[max];
    Double_t yFi23b[max];
    Double_t qFi23b[max];
    Double_t tFi23b[max];
    Double_t timeFi23b[max];
    Bool_t fFi23b[max];
    Bool_t cTofd1ls[max];
    Bool_t cTofd1rs[max];
    Bool_t cTofd2ls[max];
    Bool_t cTofd2rs[max];
    Bool_t cTofd3ls[max];
    Bool_t cTofd3rs[max];
    Bool_t maxWertFi30;
    Bool_t maxWertFi31;
    Bool_t maxWertFi32;
    Bool_t maxWertFi33;
    Int_t detTofd1l[max] = { 0 };
    Int_t detTofd2l[max] = { 0 };
    Int_t detTofd3l[max] = { 0 };
    Int_t detTofd1r[max] = { 0 };
    Int_t detTofd2r[max] = { 0 };
    Int_t detTofd3r[max] = { 0 };
    Bool_t cTofddet[max];
    Double_t xtofdtemp6[40], xtofdtemp7[40], xtofdtemp8[40], xtofdtemp9[40], xtofdtemp10[40], xtofdtemp11[40];
    Double_t ytofdtemp6[40], ytofdtemp7[40], ytofdtemp8[40], ytofdtemp9[40], ytofdtemp10[40], ytofdtemp11[40];
    Double_t qtofdtemp6[40], qtofdtemp7[40], qtofdtemp8[40], qtofdtemp9[40], qtofdtemp10[40], qtofdtemp11[40];
    Double_t ttofdtemp6[40], ttofdtemp7[40], ttofdtemp8[40], ttofdtemp9[40], ttofdtemp10[40], ttofdtemp11[40];

    for (int i = 0; i < max; i++)
    {
        fFi33[i] = false;
        fFi32[i] = false;
        fFi31[i] = false;
        fFi30[i] = false;
        fFi23a[i] = false;
        fFi23b[i] = false;
        cTofd1ls[i] = false;
        cTofd1rs[i] = false;
        cTofd2ls[i] = false;
        cTofd2rs[i] = false;
        cTofd3ls[i] = false;
        cTofd3rs[i] = false;
        xTofd1ls[i] = -1000;
        yTofd1ls[i] = -1000;
        qTofd1ls[i] = -1;
        tTofd1ls[i] = -1000;
        xTofd1rs[i] = -1000;
        yTofd1rs[i] = -1000;
        qTofd1rs[i] = -1;
        tTofd1rs[i] = -1000;
        xTofd2ls[i] = -1000;
        yTofd2ls[i] = -1000;
        qTofd2ls[i] = -1;
        tTofd2ls[i] = -1000;
        xTofd2rs[i] = -1000;
        yTofd2rs[i] = -1000;
        qTofd2rs[i] = -1;
        tTofd2rs[i] = -1000;
        xTofd3ls[i] = -1000;
        yTofd3ls[i] = -1000;
        qTofd3ls[i] = -1;
        tTofd3ls[i] = -1000;
        xTofd3rs[i] = -1000;
        yTofd3rs[i] = -1000;
        qTofd3rs[i] = -1;
        tTofd3rs[i] = -1000;
        timeFi23b[i] = -1000;
        xFi23b[i] = -1000;
        yFi23b[i] = -1000;
        qFi23b[i] = -1;
        tFi23b[i] = -1000;
        timeFi23b[i] = -1000;
        xFi23a[i] = -1000;
        yFi23a[i] = -1000;
        qFi23a[i] = -1;
        tFi23a[i] = -1000;
        timeFi23a[i] = -1000;
        xFi30[i] = -1000;
        yFi30[i] = -1000;
        qFi30[i] = -1;
        tFi30[i] = -1000;
        timeFi30[i] = -1000;
        xFi31[i] = -1000;
        yFi31[i] = -1000;
        qFi31[i] = -1;
        tFi31[i] = -1000;
        timeFi31[i] = -1000;
        xFi32[i] = -1000;
        yFi32[i] = -1000;
        qFi32[i] = -1;
        tFi32[i] = -1000;
        timeFi32[i] = -1000;
        xFi33[i] = -1000;
        yFi33[i] = -1000;
        qFi33[i] = -1;
        tFi33[i] = -1000;
        timeFi33[i] = -1000;
        detector[i] = -1;
        xdet[i] = -1000;
        ydet[i] = -1000;
        zdet[i] = -1000;
        tdet[i] = -1000;
        qdet[i] = -1;
        cTofddet[i] = false;
        qdet_s[i] = -1;
        xdet_s[i] = -1000;
        ydet_s[i] = -1000;
        tdet_s[i] = -1000;
        zdet_s[i] = -1000;
        detector_s[i] = -1;
        fibdet[i] = 0;
        fibdet_s[i] = 0;
    }

    countdet = 0;

    Int_t n_det = 100;

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
    Int_t id, id1, id2;

    Int_t det = 0;
    Int_t det1 = 0;
    Int_t det2 = 0;

    // SET THE OPTIONS ***************
    Bool_t debug_tofd = false;
    Bool_t debug_fib = false;
    Bool_t debug_fibcut = false;
    Bool_t debug_write = false;

    /*
        if (fNEvents == 22609 || fNEvents == 111385 || fNEvents == 297730 || fNEvents == 325745 || fNEvents == 336201)
        {
            debug_tofd = true;
            debug_fib = false;
            debug_fibcut = false;
            debug_write = true;
        }
    */
    Bool_t fibCuts = true;

    // only consider fiber with maximum energy deposit, only for sweep runs with beam
    Bool_t maxWerte = false;
    if (!fPairs)
    {
        maxWerte = true;
        fibCuts = false;
    }
    // if fibCuts true, dx1 (Fib3xvs3x), dx2(Fib3xvs2x), dx3(Fib3xvsTofd) used for cuts in xfib correlations

    Double_t dx1 = 2;            // cm
    Double_t dx2 = 5.;           // cm
    Double_t dx3 = 10.;          // cm
    Double_t dx4 = 0.7;          // cm
    Double_t dtft = 11.;         // ns
    Double_t dTofFib = 250.;     // ns
    Double_t dytofdfib23b = 10.; // 10 cm
                                 // END CHOOSING OPTIONS **********

    for (int i = 0; i < n_det; i++)
    {
        x[i] = -1000.;
        y[i] = -1000.;
        z[i] = -1000.;
        q[i] = -1.;
        t[i] = -1000.;

        x1[i] = -1000.;
        y1[i] = -1000.;
        z1[i] = -1000.;
        q1[i] = -1.;
        t1[i] = -1000.;

        x2[i] = -1000.;
        y2[i] = -1000.;
        z2[i] = -1000.;
        q2[i] = -1.;
        t2[i] = -1000.;

        xMax[i] = -1000.;
        yMax[i] = -1000.;
        zMax[i] = -1000.;
        qMax[i] = -1000.;
        tMax[i] = -1000.;
    }

    // is also number of ifibcount
    Int_t fi23a = 0;
    Int_t fi23b = 1;
    Int_t fi30 = 2;
    Int_t fi31 = 3;
    Int_t fi32 = 4;
    Int_t fi33 = 5;
    Int_t tofd1rc = 6;
    Int_t tofd1lc = 7;
    Int_t tofd2rc = 8;
    Int_t tofd2lc = 9;
    Int_t tofd3rc = 10;
    Int_t tofd3lc = 11;
    Int_t tofd1r = 12;
    Int_t tofd1l = 13;
    Int_t tofd2r = 14;
    Int_t tofd2l = 15;
    Int_t tofd3r = 16;
    Int_t tofd3l = 17;
    Int_t fi23ar = 18;
    Int_t fi23al = 19;

    Double_t tof = 0.;
    Bool_t pair = false;
    Bool_t single = false;
    Double_t tStart = -1000, xStart = -1000, yStart = -1000, qStart = 0;
    Bool_t first = true;
    Bool_t alpha = false;
    Bool_t carbon = false;
    Double_t x_4He = 0.;
    Double_t y_4He = 0.;
    Double_t z_4He = 0.;
    Double_t x_12C = 0.;
    Double_t y_12C = 0.;
    Double_t z_12C = 0.;

    //   cout<<"Before fib23 cut"<<endl;
    auto detHit23aa = fHitItems.at(DET_FI23A);
    Int_t nHits23aa = detHit23aa->GetEntriesFast();
    if (!fSimu && !fPairs && nHits23aa > 0)
        return;
    //   if (fPairs && nHits23aa < 1)
    //     return;
    //  cout<<"After fib23a cut"<<endl;
    auto detHit23bb = fHitItems.at(DET_FI23B);
    Int_t nHits23bb = detHit23bb->GetEntriesFast();
    if (!fSimu && !fPairs && nHits23bb > 0)
        return;
    //  if (fPairs && nHits23bb < 1)
    //    return;
    //  cout<<"After fib23b cut"<<endl;

    //   cout<<"*** Entering analysis ***"<<endl;

    // **********************************

    auto detTofd = fHitItems.at(DET_TOFD);
    Int_t nHits = detTofd->GetEntriesFast();
    if (nHits > 0)
    {
        fh_tofd_mult->Fill(nHits);
        counterTofd++;
        if (debug_tofd)
        {
            cout << "********************************" << endl;
            cout << "ToFD hits: " << nHits << endl;
        }
    }

    if (nHits > 100)
        return;

    // Trying to restore Q of tofi bar 14
    Double_t qplane2[10], qplane3[10], tplane3[10];
    for (Int_t i = 0; i < 10; i++)
    {
        qplane3[i] = -1000;
        qplane2[i] = -1000;
        tplane3[i] = -1000;
    }
    Int_t nplane2 = 0, nplane3 = 0;
    for (Int_t ihit = 0; ihit < nHits; ihit++)
    {

        R3BTofdHitData* hitTofd = (R3BTofdHitData*)detTofd->At(ihit);

        if (IS_NAN(hitTofd->GetTime()))
            continue;

        Int_t iplane = hitTofd->GetDetId();
        Int_t ibar = hitTofd->GetBarId();

        if (iplane == 2 && ibar == 24)
        {
            qplane2[nplane2] = hitTofd->GetEloss();
            nplane2++;
        }

        if (iplane == 3 && ibar == 14)
        {
            qplane3[nplane3] = hitTofd->GetEloss();
            tplane3[nplane3] = hitTofd->GetTime();
            nplane3++;
        }
    }

    Bool_t qtofi_bar14_z6[nplane3];
    for (Int_t i = 0; i < nplane3; i++)
    {
        qtofi_bar14_z6[i] = false;
    }
    for (Int_t i = 0; i < nplane2; i++)
    {
        for (Int_t j = 0; j < nplane3; j++)
        {
            fh_qtofi_vs_qtofd->Fill(qplane2[i], qplane3[j]);

            if (fGraphCuts)
            {
                if (cut_qtofi_bar14->IsInside(qplane2[i], qplane3[j]))
                    qtofi_bar14_z6[j] = true;
            }
        }
    }

    Int_t multTofd = 0;
    Int_t i1l = 0, i1r = 0, i2l = 0, i2r = 0, i3l = 0, i3r = 0;
    Int_t i1lc = 0, i1rc = 0, i2lc = 0, i2rc = 0, i3lc = 0, i3rc = 0;
    Int_t icountplane3bar14 = 0;
    // loop over ToFD

    if (debug_tofd)
        cout << "*** NEW EVENT *** " << fNEvents << endl;
    for (Int_t ihit = 0; ihit < nHits; ihit++)
    {

        R3BTofdHitData* hitTofd = (R3BTofdHitData*)detTofd->At(ihit);
        pair = false;

        if (IS_NAN(hitTofd->GetTime()))
            continue;
        /*
                if (debug_tofd)
                {
                    cout << "Hit " << ihit << " of " << nHits << " charge " << hitTofd->GetEloss() << " time "
                         << hitTofd->GetTime() << endl;
                }
        */

        Int_t iplane = hitTofd->GetDetId();
        Int_t ibar = hitTofd->GetBarId();

        //  if (iplane == 2 && ibar == 24)
        //    continue;

        Double_t qqq = hitTofd->GetEloss(); // / 1.132;

        if (iplane == 3 && ibar == 14)
        {
            if (qqq > 5 && qqq < 7 && !qtofi_bar14_z6[icountplane3bar14])
                qqq = 100; // we will use q=6 only if cut true; not to get values outside of cut we set 100

            if (qtofi_bar14_z6[icountplane3bar14])
                qqq = 6;

            if (nplane3 > 1)
                icountplane3bar14++;
        }

        if (iplane == 3 && ibar == 13)
            continue; // no charge resolution

        Double_t ttt = hitTofd->GetTime();
        fh_tofd_time->Fill(ttt);

        if (fCuts && (ttt < -150. || ttt > 100.)) // change time cuts  -150,100
        {                                         // trigger window -1500, 1500
            if (debug_tofd)
                cout << "No trigger particle!" << endl;
            continue;
        }

        // if (fB == -1710)
        //{
        //    qqq = qqq * 1.11; // / 1.132;
        //}

        if (fSimu)
        {
            // qqq = qqq * 3.2706460;
            //  if (qqq > 0.01)
            //    qqq = qqq * 26.76 * 4.4085231;
            // if (qqq < 0.01)
            //     qqq = qqq * 26.76 * 4.4085231 * 3.0599755;
            // for digiHit_s494_500AMeV_breakup_all_geant3_ideal.root
            //   qqq = (int)(qqq * 88.196293 + 1.5024036 + 0.5);
        }

        Double_t xxx = hitTofd->GetX();
        Double_t yyy = hitTofd->GetY();

        if (!fSimu)
        {
            if (iplane < 3)
            {
                yyy = yyy * (-1.); // -1 until we solve problem with y direction
                if (qqq > 1.4 && qqq < 2.6)
                {
                    yyy = (yyy + 0.7495721) * 7.6358668 / 4.7541388 - 0.0817257;
                }
            }

            if (iplane == 1)
            {
                if (qqq > 1.4 && qqq < 2.6)
                    fh_tofd_x_vs_y_z[0]->Fill(ibar, yyy);
                if (qqq > 5.5 && qqq < 6.6)
                    fh_tofd_x_vs_y_z[1]->Fill(ibar, yyy);
                if (qqq > 7.2 && qqq < 8.8)
                    fh_tofd_x_vs_y_z[2]->Fill(ibar, yyy);
            }
            if (iplane == 2)
            {
                if (qqq > 1.4 && qqq < 2.6)
                    fh_tofd_x_vs_y_z[3]->Fill(ibar, yyy);
                if (qqq > 5.5 && qqq < 6.6)
                    fh_tofd_x_vs_y_z[4]->Fill(ibar, yyy);
                if (qqq > 7.2 && qqq < 8.8)
                    fh_tofd_x_vs_y_z[5]->Fill(ibar, yyy);
            }
        }

        fh_tofd_charge->Fill(yyy, qqq);

        Double_t y_corr = 0.;
        Double_t randx;
        Double_t randy;
        randx = (std::rand() / (float)RAND_MAX) * 2.8 - 1.4;
        randy = (std::rand() / (float)RAND_MAX) * 2.8 - 1.4;
        // randx = (std::rand() / (float)RAND_MAX) - 0.5;

        fh_xy_tofd->Fill(xxx + 0. * randx, yyy);

        // first looking for the right charge
        if (1 == 1) // fB == -1250 || fB == -1710)
        {
            fh_tofd_q2_vs_q1->Fill(q1[det1], qqq);

            if (!fPairs) // fPairs given in macro, true 2 particles, false 16O tracking
            {
                if (qqq < 7.2 || qqq > 8.8)
                {
                    if (debug_tofd)
                        cout << "Not the right charge! Charge <> 8 " << qqq << endl;

                    continue;
                }
                else
                {
                    qqq = 8.; // this istemp, as hitpar from He run are not quite good any more.
                }
            }

            if (fPairs && !(qqq > 1.75 && qqq < 2.4) && !(qqq > 5.5 && qqq < 6.5))
            {
                if (debug_tofd)
                    cout << "Not the right charge! Charge = " << qqq << " Charge <> 2 and Charge <> 6" << endl;
                Nevent_wrongcharge++;
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

        if (xxx <= 0.)
        {
            // tof rechts
            if (id2 == 1)
            {
                if (maxWerte || !fibCuts)
                {
                    det2 = tofd1rc;
                }
                else
                {
                    det2 = tofd1r;
                }
            }
            else if (id2 == 2)
            {
                if (maxWerte || !fibCuts)
                {
                    det2 = tofd2rc;
                }
                else
                {
                    det2 = tofd2r;
                }
            }
            else if (id2 == 3)
            {
                if (maxWerte || !fibCuts)
                {
                    det2 = tofd3rc;
                }
                else
                {
                    det2 = tofd3r;
                }
            }
        }
        else
        {
            // tof links
            if (id2 == 1)
            {
                if (maxWerte || !fibCuts)
                {
                    det2 = tofd1lc;
                }
                else
                {
                    det2 = tofd1l;
                }
            }
            else if (id2 == 2)
            {
                if (maxWerte || !fibCuts)
                {
                    det2 = tofd2lc;
                }
                else
                {
                    det2 = tofd2l;
                }
            }
            else if (id2 == 3)
            {
                if (maxWerte || !fibCuts)
                {
                    det2 = tofd3lc;
                }
                else
                {
                    det2 = tofd3l;
                }
            }
        }

        //  cout<<"NEW HIT: "<<nHits<<", "<<ihit<<"; "<<id2<<", "<<qqq<<", "<<xxx<<endl;

        x2[det2] = xxx / 100.; // for Rene's tracker everything in meters; Dima's tracker in cm
        y2[det2] = yyy / 100. + y_corr;
        if (y2[det2] < -0.8 || y2[det2] > 0.8) // in meters!
        {
            if (debug_tofd)
                cout << "Wrong ytoFD-position! " << y2[det2] << " for charge " << q2[det2] << endl;
            continue;
        }

        z2[det2] = 0.;
        q2[det2] = qqq;
        t2[det2] = hitTofd->GetTime();

        fh_tofd_time_ac->Fill(t2[det2]);
        fh_tofd_charge_ac->Fill(y2[det2] * 100., q2[det2]);

        // register hits for tracker as long a time is in the coincidence window
        Double_t tDiffTofd = t2[det2] - t1[det1];
        fh_TimePreviousEvent->Fill(t2[det2] - t1[det1], q2[det2] + q1[det1]);

        // if (((q2[det2] + q1[det1])> 7.5 && (q2[det2] + q1[det1])< 8.5 && abs(t2[det2] - t1[det1]) < 50.) || first) //
        // change back to 5.
        Int_t detHe = -1, det12C = -1;
        if ((abs(t2[det2] - t1[det1]) < 50.) || first) // change back to 5.
        {
            // if (debug_tofd)
            //   cout << "2 particle within 5 ns   " << first << endl;
            // register point for tracker

            if (!first)
            {
                fh_tofd_q2_vs_q1_ac->Fill(q1[det1], q2[det2]);
            }

            if (debug_tofd)
            {
                cout << "registered"
                     << " x: " << x2[det2] * 100 << " y: " << y2[det2] * 100 << " q: " << q2[det2] << " t: " << t2[det2]
                     << " ID " << id2 << ", Det: " << det2 << endl;
            }

            detector[countdet] = det2;
            xdet[countdet] = x2[det2];
            ydet[countdet] = y2[det2];
            zdet[countdet] = z2[det2];
            tdet[countdet] = t2[det2];
            qdet[countdet] = (int)(q2[det2] + 0.5); // q for tracker must be integer
            cTofddet[countdet] = false;
            countdet++;

            single = true;
            first = false;
            tStart = t2[det2];
            xStart = x2[det2];
            yStart = y2[det2];
            qStart = q2[det2];

            det1 = det2;
            x1[det1] = x2[det2];
            y1[det1] = y2[det2];
            z1[det1] = 0.; // this is coordinate
            q1[det1] = q2[det2];
            t1[det1] = t2[det2];
            id1 = id2;

            // since we had a coincidence, continue with next event, if not last event.
            if (ihit < nHits - 1)
                continue;
        }
        else
        {
            x2[det2] = -1000.;
            y2[det2] = -1000.;
            z2[det2] = 0.;
            q2[det2] = 0.;
            t2[det2] = -1000.;
        }

        for (Int_t ic = 0; ic < countdet; ic++)
        {
            if (abs(qdet[ic] - 2.) < 0.5)
            {
                alpha = true;
                x_4He = xdet[ic];
                y_4He = ydet[ic];
                z_4He = zdet[ic];
                detHe = detector[ic];
            }
            if (abs(qdet[ic] - 6.) < 0.5)
            {
                carbon = true;
                x_12C = xdet[ic];
                y_12C = ydet[ic];
                z_12C = zdet[ic];
                det12C = detector[ic];
            }
            if (detector[ic] == tofd1l || detector[ic] == tofd1lc)
            {
                xTofd1ls[i1lc] = xdet[ic];
                yTofd1ls[i1lc] = ydet[ic];
                tTofd1ls[i1lc] = tdet[ic];
                qTofd1ls[i1lc] = qdet[ic];
                cTofd1ls[i1lc] = cTofddet[ic];
                i1lc++;
            }
            if (detector[ic] == tofd2l || detector[ic] == tofd2lc)
            {
                xTofd2ls[i2lc] = xdet[ic];
                yTofd2ls[i2lc] = ydet[ic];
                tTofd2ls[i2lc] = tdet[ic];
                qTofd2ls[i2lc] = qdet[ic];
                cTofd2ls[i2lc] = cTofddet[ic];
                i2lc++;
            }
            if (detector[ic] == tofd3l || detector[ic] == tofd3lc)
            {
                xTofd3ls[i3lc] = xdet[ic];
                yTofd3ls[i3lc] = ydet[ic];
                tTofd3ls[i3lc] = tdet[ic];
                qTofd3ls[i3lc] = qdet[ic];
                cTofd3ls[i3lc] = cTofddet[ic];
                i3lc++;
            }
            if (detector[ic] == tofd1r || detector[ic] == tofd1rc)
            {
                xTofd1rs[i1rc] = xdet[ic];
                yTofd1rs[i1rc] = ydet[ic];
                tTofd1rs[i1rc] = tdet[ic];
                qTofd1rs[i1rc] = qdet[ic];
                cTofd1rs[i1rc] = cTofddet[ic];
                i1rc++;
            }
            if (detector[ic] == tofd2r || detector[ic] == tofd2rc)
            {
                xTofd2rs[i2rc] = xdet[ic];
                yTofd2rs[i2rc] = ydet[ic];
                tTofd2rs[i2rc] = tdet[ic];
                qTofd2rs[i2rc] = qdet[ic];
                cTofd2rs[i2rc] = cTofddet[ic];
                i2rc++;
            }
            if (detector[ic] == tofd3r || detector[ic] == tofd3rc)
            {
                xTofd3rs[i3rc] = xdet[ic];
                yTofd3rs[i3rc] = ydet[ic];
                tTofd3rs[i3rc] = tdet[ic];
                qTofd3rs[i3rc] = qdet[ic];
                cTofd3rs[i3rc] = cTofddet[ic];
                i3rc++;
            }
        }

        if (!single)
            continue;

        if (fPairs && !(alpha && carbon))
        {
            Nevent_singles++;
            continue;
        }
        alpha = false;
        carbon = false;

        if (debug_tofd)
        {
            cout << "Number of sorted ToFD hits:" << countdet << endl;

            for (Int_t ic = 0; ic < countdet; ic++)
            {
                cout << "Det: " << detector[ic] << ", x: " << xdet[ic] * 100 << ", y: " << ydet[ic]
                     << ", t: " << tdet[ic] << ", q: " << qdet[ic] << endl;
            }
        }

        if (debug_tofd)
        {
            delta = sqrt((x_12C - x_4He) * (x_12C - x_4He) + (y_12C - y_4He) * (y_12C - y_4He));
            cout << "found good event, x_12C " << x_12C * 100. << " y_12C " << y_12C * 100. << "det: " << det12C
                 << endl;
            cout << "found good event, x_4He " << x_4He * 100. << " y_4He " << y_4He * 100. << "det: " << detHe << endl;
            //  cout << "Delta: " << delta << endl;
        }

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

        Int_t multsortedTofd = countdet;
        multTofd++;

        counterTofdMulti++;

        // cut in ToT for Fibers
        Double_t cutQ = 0.;
        if (!fPairs || fB != -1710)
            // check if cut can be applied
            cutQ = -10.;

        if (debug_fib)
            cout << "start fiber analysis" << endl;

        // loop over fiber 33
        hits33 = 0;
        hits33bc = 0;
        auto detHit33 = fHitItems.at(DET_FI33);
        Int_t nHits33 = detHit33->GetEntriesFast();
        if (debug_fib)
            cout << "Fi33 hits: " << nHits33 << endl;

        Int_t mult33 = 0, imax33 = -1;
        for (Int_t ihit33 = 0; ihit33 < nHits33; ihit33++)
        {
            det = fi33;
            R3BBunchedFiberHitData* hit33 = (R3BBunchedFiberHitData*)detHit33->At(ihit33);
            x1[det] = hit33->GetX() / 100.;
            y1[det] = hit33->GetY() / 100.;
            z1[det] = 0.;
            q1[det] = hit33->GetEloss();

            t1[det] = hit33->GetTime();
            tof = tStart - t1[det];

            // Fill histograms before cuts
            fh_Fib_ToF[det]->Fill(x1[det] * 100., tof);
            fh_xy_Fib[det]->Fill(x1[det] * 100., y1[det] * 100.);
            fh_mult_Fib[det]->Fill(nHits33);
            fh_ToT_Fib[det]->Fill(x1[det] * 100., q1[det]);
            fh_Fibs_vs_Tofd[det]->Fill(x1[tofd1r] * 100. + randx, x1[det] * 100.);
            fh_Fib_vs_Events[det]->Fill(fNEvents, x1[det] * 100.);
            fh_ToF_vs_Events[det]->Fill(fNEvents, tof);
            fh_Fib_Time[det]->Fill(x1[det] * 100., t1[det]);

            if (debug_fib)
                cout << "Fi33 bc: " << ihit33 << " x1: " << x1[det] * 100 << " y1: " << y1[det] * 100
                     << " q1: " << q1[det] << " t1: " << t1[det] << ", det: " << det << endl;

            hits33bc++;

            // Cuts on Fi33
            if (fCuts && (t1[det] < -60 || t1[det] > 50) && fPairs)
                continue;
            if (fCuts && (x1[det] < -0.3 || x1[det] > 0.3) && fPairs)
                continue;
            //  if (fCuts && (y1[det] < -0.3 || y1[det] > 0.3) && fPairs)
            //    continue;

            if (fCuts && (t1[det] < -15 || t1[det] > 15) && !fSimu && !fPairs)
                continue;
            if (fCuts && (tof < 30 || tof > 80) && !fSimu && !fPairs)
                continue;
            if (fCuts && (x1[det] * 100. < fX2min || x1[det] * 100. > fX2max) && !fPairs)
                continue;
            //            if (fGraphCuts && !cut_Fi33vsTofd->IsInside(x1[tofd1r] * 100., x1[det] * 100.))
            //                continue;

            if (debug_fib)
                cout << "Fi33 ac: " << ihit33 << " x1: " << x1[det] * 100 << " y1: " << y1[det] * 100
                     << " q1: " << q1[det] << " t1: " << t1[det] << ", det: " << det << endl;
            hits33++;

            if (q1[det] > qMax[det])
            {
                qMax[det] = q1[det];
                xMax[det] = x1[det];
                yMax[det] = y1[det];
                zMax[det] = z1[det];
                tMax[det] = t1[det];
            }

            xFi33[mult33] = x1[det];
            yFi33[mult33] = y1[det];
            /*  if (q1[det] > 7.65 && !fSimu)
                  qFi33[mult33] = 6.;
              else
                  qFi33[mult33] = 2;
              if (fSimu)*/
            qFi33[mult33] = q1[det];
            tFi33[mult33] = tof;
            timeFi33[mult33] = t1[det];
            mult33++;
            if (mult33 > 100)
                continue;

            // Fill histograms
            if (!fibCuts && !maxWerte)
            {
                fh_Fib_ToF_ac[det]->Fill(x1[det] * 100., tof);
                fh_xy_Fib_ac[det]->Fill(x1[det] * 100., y1[det] * 100.);
                fh_mult_Fib_ac[det]->Fill(mult33);
                fh_ToT_Fib_ac[det]->Fill(x1[det] * 100., q1[det]);
                fh_Fibs_vs_Tofd_ac[det]->Fill(x1[tofd1r] * 100. + randx, x1[det] * 100.);
                fh_Fib_vs_Events_ac[det]->Fill(fNEvents, x1[det] * 100.);
                fh_ToF_vs_Events_ac[det]->Fill(fNEvents, tof);
                fh_Fib_Time_ac[det]->Fill(x1[det] * 100., t1[det]);

                detector[countdet] = det;
                xdet[countdet] = x1[det];
                ydet[countdet] = y1[det];
                zdet[countdet] = z1[det];
                qdet[countdet] = q1[det];
                tdet[countdet] = t1[det];

                if (debug_fib)
                    cout << "Fi33: " << detector[countdet] << " x1: " << xdet[countdet] * 100
                         << " y1: " << ydet[countdet] * 100 << " q1: " << qdet[countdet] << " t1: " << tdet[countdet]
                         << endl;

                countdet++;
            }
        }

        if (hits33 > 0 && maxWerte)
        {
            detector[countdet] = fi33;
            xdet[countdet] = xMax[fi33];
            ydet[countdet] = yMax[fi33];
            zdet[countdet] = zMax[fi33];
            qdet[countdet] = qMax[fi33];
            tdet[countdet] = tMax[fi33];

            fh_Fib_ToF_ac[det]->Fill(xdet[countdet] * 100., tStart - tdet[countdet]);
            fh_xy_Fib_ac[det]->Fill(xdet[countdet] * 100., ydet[countdet] * 100.);
            fh_mult_Fib_ac[det]->Fill(mult33);
            fh_ToT_Fib_ac[det]->Fill(xdet[countdet] * 100., qdet[countdet]);
            fh_Fibs_vs_Tofd_ac[det]->Fill(x1[tofd1r] * 100. + randx, xdet[countdet] * 100.);
            fh_Fib_vs_Events_ac[det]->Fill(fNEvents, xdet[countdet] * 100.);
            fh_ToF_vs_Events_ac[det]->Fill(fNEvents, tStart - tdet[countdet]);
            fh_Fib_Time_ac[det]->Fill(xdet[countdet] * 100., tdet[countdet]);

            if (debug_fib)
                cout << "Fi33: " << detector[countdet] << " x1: " << xdet[countdet] * 100
                     << " y1: " << ydet[countdet] * 100 << " q1: " << qdet[countdet] << " t1: " << tdet[countdet]
                     << endl;

            countdet++;
            if (hits33 > 0)
                maxWertFi33 = true;
        }

        // loop over fiber 31
        hits31 = 0;
        hits31bc = 0;
        auto detHit31 = fHitItems.at(DET_FI31);
        Int_t nHits31 = detHit31->GetEntriesFast();
        if (debug_fib)
            cout << "Fi31 hits: " << nHits31 << endl;
        Int_t mult31 = 0;
        for (Int_t ihit31 = 0; ihit31 < nHits31; ihit31++)
        {
            det = fi31;
            R3BBunchedFiberHitData* hit31 = (R3BBunchedFiberHitData*)detHit31->At(ihit31);
            x1[det] = hit31->GetX() / 100.;
            y1[det] = hit31->GetY() / 100.;
            z1[det] = 0.;
            q1[det] = hit31->GetEloss();

            t1[det] = hit31->GetTime();
            tof = tStart - t1[det];

            // Fill histograms before cuts
            fh_Fib_ToF[det]->Fill(x1[det] * 100., tof);
            fh_xy_Fib[det]->Fill(x1[det] * 100., y1[det] * 100);
            fh_mult_Fib[det]->Fill(nHits31);
            fh_ToT_Fib[det]->Fill(x1[det] * 100., q1[det]);
            fh_Fibs_vs_Tofd[det]->Fill(x1[tofd1r] * 100. + randx, x1[det] * 100.);
            fh_Fib_vs_Events[det]->Fill(fNEvents, x1[det] * 100.);
            fh_ToF_vs_Events[det]->Fill(fNEvents, tof);
            fh_Fib_Time[det]->Fill(x1[det] * 100., t1[det]);

            if (debug_fib)
                cout << "Fi31 bc: " << ihit31 << " x1: " << x1[det] * 100 << " y1: " << y1[det] * 100
                     << " q1: " << q1[det] << " t1: " << t1[det] << ", det: " << det << endl;
            hits31bc++;

            // Cuts on Fi31
            if (fCuts && (t1[det] < -60 || t1[det] > 50) && fPairs)
                continue;
            if (fCuts && (x1[det] < -0.3 || x1[det] > 0.3) && fPairs)
                continue;
            //   if (fCuts && (y1[det] < -0.3 || y1[det] > 0.3) && fPairs)
            //     continue;

            if (fCuts && (t1[det] < -15 || t1[det] > 15) && !fSimu && !fPairs)
                continue;
            if (fCuts && (tof < 30 || tof > 80) && !fSimu && !fPairs)
                continue;
            if (fCuts && (x1[det] * 100. < fX1min || x1[det] * 100. > fX1max) && !fPairs)
                continue;
            //            if (fGraphCuts && !cut_Fi33vsTofd->IsInside(x1[tofd1r] * 100., x1[det] * 100.))
            //                continue;

            if (debug_fib)
                cout << "Fi31 ac: " << ihit31 << " x1: " << x1[det] * 100 << " y1: " << y1[det] * 100
                     << " q1: " << q1[det] << " t1: " << t1[det] << ", det: " << det << endl;

            hits31++;

            xFi31[mult31] = x1[det];
            yFi31[mult31] = y1[det];
            /*  if (q1[det] > 7.6 && !fSimu)
              {
                  qFi31[mult31] = 6.;
              }
              else
              {
                  qFi31[mult31] = 2;
              }
              if (fSimu) */
            qFi31[mult31] = q1[det];
            tFi31[mult31] = tof;
            timeFi31[mult31] = t1[det];

            //  cout << "Fi31 ac: " << mult31 << " x1: " << xFi31[mult31] << " y1: " << yFi31[mult31]  << " q1: " <<
            //  qFi31[mult31]
            //     << " t1: " << timeFi31[mult31] <<", "<<q1[tofd1r]<<", "<<q1[tofd2r]<< endl;

            mult31++;

            if (mult31 > 100)
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
            if (!fibCuts && !maxWerte)
            {
                fh_Fib_ToF_ac[det]->Fill(x1[det] * 100., tof);
                fh_xy_Fib_ac[det]->Fill(x1[det] * 100., y1[det] * 100.);
                fh_mult_Fib_ac[det]->Fill(mult31);
                fh_ToT_Fib_ac[det]->Fill(x1[det] * 100., q1[det]);
                fh_Fibs_vs_Tofd_ac[det]->Fill(x1[tofd1r] * 100. + randx, x1[det] * 100.);
                fh_Fib_vs_Events_ac[det]->Fill(fNEvents, x1[det] * 100.);
                fh_ToF_vs_Events_ac[det]->Fill(fNEvents, tof);
                fh_Fib_Time_ac[det]->Fill(x1[det] * 100., t1[det]);

                detector[countdet] = det;
                xdet[countdet] = x1[det];
                ydet[countdet] = y1[det];
                zdet[countdet] = z1[det];
                qdet[countdet] = q1[det];
                tdet[countdet] = t1[det];

                if (debug_fib)
                    cout << "Fi31: " << detector[countdet] << " x1: " << xdet[countdet] << " y1: " << ydet[countdet]
                         << " q1: " << qdet[countdet] << " t1: " << tdet[countdet] << endl;

                countdet++;
            }
        }

        if (hits31 > 0 && maxWerte)
        {
            detector[countdet] = fi31;
            xdet[countdet] = xMax[fi31];
            ydet[countdet] = yMax[fi31];
            zdet[countdet] = zMax[fi31];
            qdet[countdet] = qMax[fi31];
            tdet[countdet] = tMax[fi31];

            fh_Fib_ToF_ac[det]->Fill(xdet[countdet] * 100., tStart - tdet[countdet]);
            fh_xy_Fib_ac[det]->Fill(xdet[countdet] * 100., ydet[countdet] * 100.);
            fh_mult_Fib_ac[det]->Fill(mult33);
            fh_ToT_Fib_ac[det]->Fill(xdet[countdet] * 100., qdet[countdet]);
            fh_Fibs_vs_Tofd_ac[det]->Fill(x1[tofd1r] * 100. + randx, xdet[countdet] * 100.);
            fh_Fib_vs_Events_ac[det]->Fill(fNEvents, xdet[countdet] * 100.);
            fh_ToF_vs_Events_ac[det]->Fill(fNEvents, tStart - tdet[countdet]);
            fh_Fib_Time_ac[det]->Fill(xdet[countdet] * 100., tdet[countdet]);

            if (debug_fib)
                cout << "Fi31: " << detector[countdet] << " x1: " << xdet[countdet] << " y1: " << ydet[countdet]
                     << " q1: " << qdet[countdet] << " t1: " << tdet[countdet] << endl;

            countdet++;
            if (hits31 > 0)
                maxWertFi31 = true;
        }

        // loop over fiber 32
        hits32 = 0;
        hits32bc = 0;
        auto detHit32 = fHitItems.at(DET_FI32);
        Int_t nHits32 = detHit32->GetEntriesFast();
        if (debug_fib)
            cout << "Fi32 hits: " << nHits32 << endl;
        Int_t mult32 = 0;
        for (Int_t ihit32 = 0; ihit32 < nHits32; ihit32++)
        {
            det = fi32;
            R3BBunchedFiberHitData* hit32 = (R3BBunchedFiberHitData*)detHit32->At(ihit32);
            x1[det] = hit32->GetX() / 100.;
            y1[det] = hit32->GetY() / 100.;
            z1[det] = 0.;
            q1[det] = hit32->GetEloss();

            t1[det] = hit32->GetTime();
            tof = tStart - t1[det];

            // Fill histograms before cuts
            fh_Fib_ToF[det]->Fill(x1[det] * 100., tof);
            fh_xy_Fib[det]->Fill(x1[det] * 100., y1[det] * 100.);
            fh_mult_Fib[det]->Fill(nHits32);
            fh_ToT_Fib[det]->Fill(x1[det] * 100., q1[det]);
            fh_Fibs_vs_Tofd[det]->Fill(x1[tofd1l] * 100. + randx, x1[det] * 100.);
            fh_Fib_vs_Events[det]->Fill(fNEvents, x1[det] * 100.);
            fh_ToF_vs_Events[det]->Fill(fNEvents, tof);
            fh_Fib_Time[det]->Fill(x1[det] * 100., t1[det]);

            hits32bc++;

            if (debug_fib)
                cout << "Fi32 bc: " << ihit32 << " x1: " << x1[det] * 100 << " y1: " << y1[det] * 100
                     << " q1: " << q1[det] << " t1: " << t1[det] << ", det: " << det << endl;

            // Cuts on Fi32
            if (fCuts && (t1[det] < -60 || t1[det] > 50) && fPairs)
                continue;
            if (fCuts && (x1[det] < -0.3 || x1[det] > 0.3) && fPairs)
                continue;
            //  if (fCuts && (y1[det] < -0.3 || y1[det] > 0.3) && fPairs)
            //    continue;

            if (fCuts && (t1[det] < -20 || t1[det] > 30) && !fSimu && !fPairs)
                continue;
            if (fCuts && (tof < 20 || tof > 50) && !fSimu && !fPairs)
                continue;
            if (fCuts && (x1[det] * 100. < fX2min || x1[det] * 100. > fX2max) && !fPairs)
                continue;

            hits32++;

            if (debug_fib)
                cout << "Fi32 ac: " << ihit32 << " x1: " << x1[det] * 100 << " y1: " << y1[det] * 100
                     << " q1: " << q1[det] << " t1: " << t1[det] << ", det: " << det << endl;

            xFi32[mult32] = x1[det];
            yFi32[mult32] = y1[det];
            /*  if (q1[det] > 7.6 && !fSimu)
                  qFi32[mult32] = 6.;
              else
                  qFi32[mult32] = 2.;
              if (fSimu)*/
            qFi32[mult32] = q1[det];
            tFi32[mult32] = tof;
            timeFi32[mult32] = t1[det];
            mult32++;

            if (mult32 > 100)
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
            if (!fibCuts && !maxWerte)
            {
                fh_Fib_ToF_ac[det]->Fill(x1[det] * 100., tof);
                fh_xy_Fib_ac[det]->Fill(x1[det] * 100., y1[det] * 100.);
                fh_mult_Fib_ac[det]->Fill(mult32);
                fh_ToT_Fib_ac[det]->Fill(x1[det] * 100., q1[det]);
                fh_Fibs_vs_Tofd_ac[det]->Fill(x1[tofd1r] * 100. + randx, x1[det] * 100.);
                fh_Fib_vs_Events_ac[det]->Fill(fNEvents, x1[det] * 100.);
                fh_ToF_vs_Events_ac[det]->Fill(fNEvents, tof);
                fh_Fib_Time_ac[det]->Fill(x1[det] * 100., t1[det]);

                detector[countdet] = det;
                xdet[countdet] = x1[det];
                ydet[countdet] = y1[det];
                zdet[countdet] = z1[det];
                qdet[countdet] = q1[det];
                tdet[countdet] = t1[det];

                if (debug_fib)
                    cout << "Fi32: " << detector[countdet] << " x1: " << xdet[countdet] << " y1: " << ydet[countdet]
                         << " q1: " << qdet[countdet] << " t1: " << tdet[countdet] << endl;

                countdet++;
            }
        }

        if (hits32 > 0 && maxWerte)
        {
            detector[countdet] = fi32;
            xdet[countdet] = xMax[fi32];
            ydet[countdet] = yMax[fi32];
            zdet[countdet] = zMax[fi32];
            qdet[countdet] = qMax[fi32];
            tdet[countdet] = tMax[fi32];

            fh_Fib_ToF_ac[det]->Fill(xdet[countdet] * 100., tStart - tdet[countdet]);
            fh_xy_Fib_ac[det]->Fill(xdet[countdet] * 100., ydet[countdet] * 100.);
            fh_mult_Fib_ac[det]->Fill(mult33);
            fh_ToT_Fib_ac[det]->Fill(xdet[countdet] * 100., qdet[countdet]);
            fh_Fibs_vs_Tofd_ac[det]->Fill(x1[tofd1r] * 100. + randx, xdet[countdet] * 100.);
            fh_Fib_vs_Events_ac[det]->Fill(fNEvents, xdet[countdet] * 100.);
            fh_ToF_vs_Events_ac[det]->Fill(fNEvents, tStart - tdet[countdet]);
            fh_Fib_Time_ac[det]->Fill(xdet[countdet] * 100., tdet[countdet]);

            if (debug_fib)
                cout << "Fi32: " << detector[countdet] << " x1: " << xdet[countdet] << " y1: " << ydet[countdet]
                     << " q1: " << qdet[countdet] << " t1: " << tdet[countdet] << endl;

            countdet++;
            if (hits32 > 0)
                maxWertFi32 = true;
        }

        // loop over fiber 30
        hits30 = 0;
        hits30bc = 0;
        auto detHit30 = fHitItems.at(DET_FI30);
        Int_t nHits30 = detHit30->GetEntriesFast();
        if (debug_fib)
            cout << "Fi30 hits: " << nHits30 << endl;
        Int_t mult30 = 0;
        for (Int_t ihit30 = 0; ihit30 < nHits30; ihit30++)
        {
            det = fi30;
            R3BBunchedFiberHitData* hit30 = (R3BBunchedFiberHitData*)detHit30->At(ihit30);
            x1[det] = hit30->GetX() / 100.;
            y1[det] = hit30->GetY() / 100.;
            z1[det] = 0.;
            q1[det] = hit30->GetEloss();

            t1[det] = hit30->GetTime();
            tof = tStart - t1[det];

            // Fill histograms before cuts
            fh_Fib_ToF[det]->Fill(x1[det] * 100., tof);
            fh_xy_Fib[det]->Fill(x1[det] * 100., y1[det] * 100.);
            fh_mult_Fib[det]->Fill(nHits30);
            fh_ToT_Fib[det]->Fill(x1[det] * 100., q1[det]);
            fh_Fibs_vs_Tofd[det]->Fill(x1[tofd1l] * 100. + randx, x1[det] * 100.);
            fh_Fib_vs_Events[det]->Fill(fNEvents, x1[det] * 100.);
            fh_ToF_vs_Events[det]->Fill(fNEvents, tof);
            fh_Fib_Time[det]->Fill(x1[det] * 100., t1[det]);

            hits30bc++;

            if (debug_fib)
                cout << "Fi30 bc: " << ihit30 << " x1: " << x1[det] * 100 << " y1: " << y1[det] * 100
                     << " q1: " << q1[det] << " t1: " << t1[det] << ", det: " << det << endl;

            // Cuts on Fi30
            if (fCuts && (t1[det] < -50 || t1[det] > 40) && fPairs)
                continue;
            if (fCuts && (x1[det] < -0.3 || x1[det] > 0.3) && fPairs)
                continue;
            //  if (fCuts && (y1[det] < -0.3 || y1[det] > 0.3) && fPairs)
            //    continue;

            if (fCuts && (t1[det] < -20 || t1[det] > 30) && !fSimu && !fPairs)
                continue;
            if (fCuts && (tof < 20 || tof > 50) && !fSimu && !fPairs)
                continue;
            if (fCuts && (x1[det] * 100. < fX1min || x1[det] * 100. > fX1max) && !fPairs)
                continue;

            hits30++;
            if (debug_fib)
                cout << "Fi30 ac: " << ihit30 << " x1: " << x1[det] * 100 << " y1: " << y1[det] * 100
                     << " q1: " << q1[det] << " t1: " << t1[det] << ", det: " << det << endl;

            xFi30[mult30] = x1[det];
            yFi30[mult30] = y1[det];
            qFi30[mult30] = q1[det];
            tFi30[mult30] = tof;
            timeFi30[mult30] = t1[det];
            mult30++;
            if (mult30 > 100)
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
            if (!fibCuts && !maxWerte)
            {
                fh_Fib_ToF_ac[det]->Fill(x1[det] * 100., tof);
                fh_xy_Fib_ac[det]->Fill(x1[det] * 100., y1[det] * 100.);
                fh_mult_Fib_ac[det]->Fill(mult30);
                fh_ToT_Fib_ac[det]->Fill(x1[det] * 100., q1[det]);
                fh_Fibs_vs_Tofd_ac[det]->Fill(x1[tofd1r] * 100. + randx, x1[det] * 100.);
                fh_Fib_vs_Events_ac[det]->Fill(fNEvents, x1[det] * 100.);
                fh_ToF_vs_Events_ac[det]->Fill(fNEvents, tof);
                fh_Fib_Time_ac[det]->Fill(x1[det] * 100., t1[det]);

                detector[countdet] = det;
                xdet[countdet] = x1[det];
                ydet[countdet] = y1[det];
                zdet[countdet] = z1[det];
                qdet[countdet] = q1[det];
                tdet[countdet] = t1[det];

                if (debug_fib)
                    cout << "Fi30: " << detector[countdet] << " x1: " << xdet[countdet] << " y1: " << ydet[countdet]
                         << " q1: " << qdet[countdet] << " t1: " << tdet[countdet] << endl;
                countdet++;
            }
        }

        if (hits30 > 0 && maxWerte)
        {
            detector[countdet] = fi30;
            xdet[countdet] = xMax[fi30];
            ydet[countdet] = yMax[fi30];
            zdet[countdet] = zMax[fi30];
            qdet[countdet] = qMax[fi30];
            qdet[countdet] = qdet[fi30];
            tdet[countdet] = tMax[fi30];

            fh_Fib_ToF_ac[det]->Fill(xdet[countdet] * 100., tStart - tdet[countdet]);
            fh_xy_Fib_ac[det]->Fill(xdet[countdet] * 100., ydet[countdet] * 100.);
            fh_mult_Fib_ac[det]->Fill(mult33);
            fh_ToT_Fib_ac[det]->Fill(xdet[countdet] * 100., qdet[countdet]);
            fh_Fibs_vs_Tofd_ac[det]->Fill(x1[tofd1r] * 100. + randx, xdet[countdet] * 100.);
            fh_Fib_vs_Events_ac[det]->Fill(fNEvents, xdet[countdet] * 100.);
            fh_ToF_vs_Events_ac[det]->Fill(fNEvents, tStart - tdet[countdet]);
            fh_Fib_Time_ac[det]->Fill(xdet[countdet] * 100., tdet[countdet]);

            if (debug_fib)
                cout << "Fi30: " << detector[countdet] << " x1: " << xdet[countdet] << " y1: " << ydet[countdet]
                     << " q1: " << qdet[countdet] << " t1: " << tdet[countdet] << endl;

            countdet++;
            maxWertFi30 = true;
        }

        // loop over fiber 23a
        Int_t mult23a = 0, imax23a = -1;
        auto detHit23a = fHitItems.at(DET_FI23A);
        Int_t nHits23a = detHit23a->GetEntriesFast();
        if (debug_fib)
            cout << "Fi23a hits: " << nHits23a << endl;
        // fh_mult_Fib[fi23a]->Fill(nHits23a);
        for (Int_t ihit23a = 0; ihit23a < nHits23a; ihit23a++)
        {
            det = fi23a;
            R3BBunchedFiberHitData* hit23a = (R3BBunchedFiberHitData*)detHit23a->At(ihit23a);
            x1[det] = hit23a->GetX() / 100.;
            y1[det] = hit23a->GetY() / 100.;
            z1[det] = 0.;
            q1[det] = hit23a->GetEloss();

            // cout << "Fib23a x: " << hit23a->GetX() << " y: " << hit23a->GetY() << endl;

            /*  if (fSimu)
              {
                  if (q1[det] > 9.)
                      q1[det] = 8.;
                  else if (q1[det] > 3.)
                      q1[det] = 6.;
                  else if (q1[det] > 0.)
                      q1[det] = 2.;
              }*/

            t1[det] = hit23a->GetTime();
            tof = tStart - t1[det];

            //     cout << "Fi23 bc: " << ihit23a << " x1: " << x1[det]*100 << " y1: " << y1[det]*100 << " q1: " <<
            //     q1[det]
            //            << " t1: " << t1[det] << endl;

            // Fill histograms before cuts
            fh_Fib_ToF[det]->Fill(x1[det] * 100., tof);
            fh_xy_Fib[det]->Fill(x1[det] * 100., y1[det] * 100.);
            fh_mult_Fib[det]->Fill(nHits23a);
            fh_ToT_Fib[det]->Fill(x1[det] * 100., q1[det]);
            fh_Fibs_vs_Tofd[det]->Fill(xStart * 100. + randx, x1[det] * 100.);
            fh_Fib_vs_Events[det]->Fill(fNEvents, x1[det] * 100.);
            fh_ToF_vs_Events[det]->Fill(fNEvents, tof);
            fh_Fib_Time[det]->Fill(x1[det] * 100., t1[det]);

            if (debug_fib)
                cout << "Fi23a bc: " << ihit23a << " x1: " << x1[det] * 100 << " y1: " << y1[det] * 100
                     << " q1: " << q1[det] << " t1: " << t1[det] << ", det: " << det << endl;

            Double_t tmem = 0;
            if (ihit23a > 0)
            {
                fh_TimePreviousEventFi23a->Fill(tDiffTofd, tmem - t1[det1]);
            }
            tmem = t1[det1];

            // Cuts on Fi23a
            // if (fCuts && (y1[det] * 100. < -7.5 || y1[det] * 100. > 7.5) && !fSimu)
            //   continue;
            //   if (fCuts && (t1[det] < -40 || t1[det] > 40) && !fSimu && fB == -1710)
            //       continue;
            if (fCuts && (abs(x1[det] * 100.) < 0.1512))
                continue;
            //  if (fCuts && (x1[det] < -0.06 || x1[det] > 0.06))
            //    continue;
            // if (fCuts && (y1[det] < -0.1 || y1[det] > 0.1) && fPairs)
            //   continue;

            if (fCuts && (y1[det] < -0.03 || y1[det] > 0.03) && !fSimu && !fPairs)
                continue;
            if (fCuts && (tof < 10 || tof > 90) && !fSimu && !fPairs)
                continue;

            if (debug_fib)
                cout << "Fi23a ac: " << ihit23a << " x1: " << x1[det] * 100 << " y1: " << y1[det] * 100
                     << " q1: " << q1[det] << " t1: " << t1[det] << ", det: " << det << endl;

            if (q1[det] > qMax[det])
            {
                qMax[det] = q1[det];
                xMax[det] = x1[det];
                yMax[det] = y1[det];
                zMax[det] = z1[det];
                tMax[det] = t1[det];
                imax23a = mult23a;
            }

            xFi23a[mult23a] = x1[det];
            yFi23a[mult23a] = y1[det];
            qFi23a[mult23a] = q1[det];
            tFi23a[mult23a] = tof;
            timeFi23a[mult23a] = t1[det];

            mult23a++;
            if (mult23a > 100)
                continue;

            // Fill histograms
            if (!fibCuts && !maxWerte)
            {
                fh_Fib_ToF_ac[det]->Fill(x1[det] * 100., tof);
                fh_xy_Fib_ac[det]->Fill(x1[det] * 100., y1[det] * 100.);
                fh_mult_Fib_ac[det]->Fill(mult23a);
                fh_ToT_Fib_ac[det]->Fill(x1[det] * 100., q1[det]);
                fh_Fibs_vs_Tofd_ac[det]->Fill(x1[tofd1r] * 100. + randx, x1[det] * 100.);
                fh_Fib_vs_Events_ac[det]->Fill(fNEvents, x1[det] * 100.);
                fh_ToF_vs_Events_ac[det]->Fill(fNEvents, tof);
                fh_Fib_Time_ac[det]->Fill(x1[det] * 100., t1[det]);

                detector[countdet] = det;
                xdet[countdet] = x1[det];
                ydet[countdet] = y1[det];
                zdet[countdet] = z1[det];
                qdet[countdet] = q1[det];
                tdet[countdet] = t1[det];

                if (debug_fib)
                    cout << "Fi23a: " << detector[countdet] << " x1: " << xdet[countdet] << " y1: " << ydet[countdet]
                         << " q1: " << qdet[countdet] << " t1: " << tdet[countdet] << endl;

                countdet++;
            }
        }
        /*
                for(Int_t i23a=0; i23a<mult23a;i23a++)
                {
                    if(i23a == imax23a) qFi23a[i23a] = 6;
                    else qFi23a[i23a] = 2;
                }
        */
        if (mult23a > 0 && maxWerte)
        {
            detector[countdet] = fi23a;
            xdet[countdet] = xMax[fi23a];
            ydet[countdet] = yMax[fi23a];
            zdet[countdet] = zMax[fi23a];
            qdet[countdet] = 8;
            tdet[countdet] = tMax[fi23a];

            fh_Fib_ToF_ac[det]->Fill(xdet[countdet] * 100., tStart - tdet[countdet]);
            fh_xy_Fib_ac[det]->Fill(xdet[countdet] * 100., ydet[countdet] * 100.);
            fh_mult_Fib_ac[det]->Fill(mult33);
            fh_ToT_Fib_ac[det]->Fill(xdet[countdet] * 100., qdet[countdet]);
            fh_Fibs_vs_Tofd_ac[det]->Fill(x1[tofd1r] * 100. + randx, xdet[countdet] * 100.);
            fh_Fib_vs_Events_ac[det]->Fill(fNEvents, xdet[countdet] * 100.);
            fh_ToF_vs_Events_ac[det]->Fill(fNEvents, tStart - tdet[countdet]);
            fh_Fib_Time_ac[det]->Fill(xdet[countdet] * 100., tdet[countdet]);

            if (debug_fib)
                cout << "Fi23a: " << detector[countdet] << " x1: " << xdet[countdet] << " y1: " << ydet[countdet]
                     << " q1: " << qdet[countdet] << " t1: " << tdet[countdet] << endl;

            countdet++;
        }

        // loop over fiber 23b
        Int_t mult23b = 0, imax23b = -1;
        auto detHit23b = fHitItems.at(DET_FI23B);
        Int_t nHits23b = detHit23b->GetEntriesFast();
        if (debug_fib)
            cout << "Fi23b hits: " << nHits23b << endl;
        //   fh_mult_Fib[fi23b]->Fill(nHits23b);
        for (Int_t ihit23b = 0; ihit23b < nHits23b; ihit23b++)
        {
            det = fi23b;
            R3BBunchedFiberHitData* hit23b = (R3BBunchedFiberHitData*)detHit23b->At(ihit23b);
            x1[det] = hit23b->GetX() / 100.;
            y1[det] = hit23b->GetY() / 100.;
            z1[det] = 0.;
            q1[det] = hit23b->GetEloss();

            // cout << "Fib23b x: " << hit23b->GetX() << " y: " << hit23b->GetY() << endl;

            /*  if (fSimu)
              {
                  if (q1[det] > 9.)
                      q1[det] = 8.;
                  else if (q1[det] > 3.)
                      q1[det] = 6.;
                  else if (q1[det] > 0.)
                      q1[det] = 2.;
              }*/

            t1[det] = hit23b->GetTime();
            tof = tStart - t1[det];

            // Fill histograms before cuts
            fh_Fib_ToF[det]->Fill(y1[det] * 100., tof);
            fh_xy_Fib[det]->Fill(x1[det] * 100., y1[det] * 100.);
            fh_mult_Fib[det]->Fill(nHits23b);
            fh_ToT_Fib[det]->Fill(y1[det] * 100., q1[det]);
            fh_Fibs_vs_Tofd[det]->Fill(yStart * 100., y1[det] * 100.);
            fh_Fib_vs_Events[det]->Fill(fNEvents, y1[det] * 100.);
            fh_ToF_vs_Events[det]->Fill(fNEvents, tof);
            fh_Fib_Time[det]->Fill(y1[det] * 100., t1[det]);
            Double_t tmem = 0;
            if (ihit23b > 0)
            {
                fh_TimePreviousEventFi23b->Fill(tDiffTofd, tmem - t1[det1]);
            }
            tmem = t1[det1];
            if (debug_fib)
                cout << "Fi23b bc: " << ihit23b << " x1: " << x1[det] * 100. << " y1: " << y1[det] * 100.
                     << " q1: " << q1[det] << " t1: " << t1[det] << ", det: " << det << endl;

            // Cuts on Fi23b
            // if (fCuts && (x1[det] * 100. < -7.5 || x1[det] * 100. > 7.5) && !fSimu)
            //   continue;
            //   if (fCuts && (t1[det] < -40 || t1[det] > 40) && !fSimu && fB == -1710)
            //      continue;
            if (fCuts && (abs(y1[det] * 100.) < 0.1512))
                continue;

            if (fCuts && (y1[det] < -0.03 || y1[det] > 0.03) && !fSimu && !fPairs)
                continue;
            if (fCuts && (tof < 10 || tof > 90) && !fSimu && !fPairs)
                continue;

            if (q1[det] > qMax[det])
            {
                qMax[det] = q1[det];
                xMax[det] = x1[det];
                yMax[det] = y1[det];
                zMax[det] = z1[det];
                tMax[det] = t1[det];
                imax23b = mult23b;
            }

            xFi23b[mult23b] = x1[det];
            yFi23b[mult23b] = y1[det];
            qFi23b[mult23b] = q1[det];
            tFi23b[mult23b] = tof;
            timeFi23b[mult23b] = t1[det];
            mult23b++;

            if (mult23b > 100)
                continue;

            // Fill histograms
            if (!fibCuts && !maxWerte && mult23b > 0)
            {
                fh_Fib_ToF_ac[det]->Fill(x1[det] * 100., tof);
                fh_xy_Fib_ac[det]->Fill(x1[det] * 100., y1[det] * 100.);
                fh_mult_Fib_ac[det]->Fill(mult23b);
                fh_ToT_Fib_ac[det]->Fill(x1[det] * 100., q1[det]);
                fh_Fibs_vs_Tofd_ac[det]->Fill(x1[tofd1r] * 100. + randx, x1[det] * 100.);
                fh_Fib_vs_Events_ac[det]->Fill(fNEvents, x1[det] * 100.);
                fh_ToF_vs_Events_ac[det]->Fill(fNEvents, tof);
                fh_Fib_Time_ac[det]->Fill(x1[det] * 100., t1[det]);

                detector[countdet] = det;
                xdet[countdet] = x1[det];
                ydet[countdet] = y1[det];
                zdet[countdet] = z1[det];
                qdet[countdet] = 0;
                tdet[countdet] = t1[det];

                if (debug_fib)
                    cout << "Fi23b: " << detector[countdet] << " x1: " << xdet[countdet] << " y1: " << ydet[countdet]
                         << " q1: " << qdet[countdet] << " t1: " << tdet[countdet] << endl;
                countdet++;
            }

            if (debug_fib)
                cout << "Fi23b ac: " << ihit23b << " x1: " << x1[det] * 100. << " y1: " << y1[det] * 100.
                     << " q1: " << q1[det] << " t1: " << t1[det] << ", det: " << det << endl;
        }
        /*
               for(Int_t i23b=0; i23b<mult23b;i23b++)
               {
                   if(i23b == imax23b) qFi23b[i23b] = 6;
                   else qFi23b[i23b] = 2;
               }
        */

        if (mult23b > 0 && maxWerte)
        {
            detector[countdet] = fi23b;
            xdet[countdet] = xMax[fi23b];
            ydet[countdet] = yMax[fi23b];
            zdet[countdet] = zMax[fi23b];
            qdet[countdet] = 8;
            tdet[countdet] = tMax[fi23b];

            fh_Fib_ToF_ac[det]->Fill(xdet[countdet] * 100., tStart - tdet[countdet]);
            fh_xy_Fib_ac[det]->Fill(xdet[countdet] * 100., ydet[countdet] * 100.);
            fh_mult_Fib_ac[det]->Fill(mult33);
            fh_ToT_Fib_ac[det]->Fill(xdet[countdet] * 100., qdet[countdet]);
            fh_Fibs_vs_Tofd_ac[det]->Fill(x1[tofd1r] * 100. + randx, xdet[countdet] * 100.);
            fh_Fib_vs_Events_ac[det]->Fill(fNEvents, xdet[countdet] * 100.);
            fh_ToF_vs_Events_ac[det]->Fill(fNEvents, tStart - tdet[countdet]);
            fh_Fib_Time_ac[det]->Fill(xdet[countdet] * 100., tdet[countdet]);

            if (debug_fib)
                cout << "Fi23b: " << detector[countdet] << " x1: " << xdet[countdet] << " y1: " << ydet[countdet]
                     << " q1: " << qdet[countdet] << " t1: " << tdet[countdet] << endl;

            countdet++;
        }

        // if we have magnetic field runs we do not have hits in Fib23
        if (!fSimu && nHits23a < 1 && nHits23b < 1 && !fPairs && fB != -1710)
        {
            det = fi23a;
            detector[countdet] = det;
            xdet[countdet] = 0.;
            ydet[countdet] = 0.;
            zdet[countdet] = 0.;
            qdet[countdet] = 8.;
            tdet[countdet] = 0.;

            fh_Fib_ToF_ac[det]->Fill(xdet[countdet] * 100., tStart - tdet[countdet]);
            fh_xy_Fib_ac[det]->Fill(xdet[countdet] * 100., ydet[countdet] * 100.);
            fh_mult_Fib_ac[det]->Fill(mult23a);
            fh_ToT_Fib_ac[det]->Fill(xdet[countdet] * 100., qdet[countdet]);
            fh_Fibs_vs_Tofd_ac[det]->Fill(x1[tofd1r] * 100. + randx, xdet[countdet] * 100.);
            fh_Fib_vs_Events_ac[det]->Fill(fNEvents, xdet[countdet] * 100.);
            fh_ToF_vs_Events_ac[det]->Fill(fNEvents, tStart - tdet[countdet]);
            fh_Fib_Time_ac[det]->Fill(xdet[countdet] * 100., tdet[countdet]);

            countdet++;

            det = fi23b;
            detector[countdet] = det;
            xdet[countdet] = 0.;
            ydet[countdet] = 0.;
            zdet[countdet] = 0.;
            qdet[countdet] = 8.;
            tdet[countdet] = 0.;

            fh_Fib_ToF_ac[det]->Fill(xdet[countdet] * 100., tStart - tdet[countdet]);
            fh_xy_Fib_ac[det]->Fill(xdet[countdet] * 100., ydet[countdet] * 100.);
            fh_mult_Fib_ac[det]->Fill(mult23a);
            fh_ToT_Fib_ac[det]->Fill(xdet[countdet] * 100., qdet[countdet]);
            fh_Fibs_vs_Tofd_ac[det]->Fill(x1[tofd1r] * 100. + randx, xdet[countdet] * 100.);
            fh_Fib_vs_Events_ac[det]->Fill(fNEvents, xdet[countdet] * 100.);
            fh_ToF_vs_Events_ac[det]->Fill(fNEvents, tStart - tdet[countdet]);
            fh_Fib_Time_ac[det]->Fill(xdet[countdet] * 100., tdet[countdet]);

            countdet++;
        }
        if (fSimu && nHits23a < 1 && !fPairs && fB != -1710)
        {
            det = fi23a;
            detector[countdet] = det;
            xdet[countdet] = 0.;
            ydet[countdet] = 0.;
            zdet[countdet] = 0.;
            qdet[countdet] = 8.;
            tdet[countdet] = 0.;

            fh_Fib_ToF_ac[det]->Fill(xdet[countdet] * 100., tStart - tdet[countdet]);
            fh_xy_Fib_ac[det]->Fill(xdet[countdet] * 100., ydet[countdet] * 100.);
            fh_mult_Fib_ac[det]->Fill(mult23a);
            fh_ToT_Fib_ac[det]->Fill(xdet[countdet] * 100., qdet[countdet]);
            fh_Fibs_vs_Tofd_ac[det]->Fill(x1[tofd1r] * 100. + randx, xdet[countdet] * 100.);
            fh_Fib_vs_Events_ac[det]->Fill(fNEvents, xdet[countdet] * 100.);
            fh_ToF_vs_Events_ac[det]->Fill(fNEvents, tStart - tdet[countdet]);
            fh_Fib_Time_ac[det]->Fill(xdet[countdet] * 100., tdet[countdet]);

            countdet++;
        }
        if (fSimu && nHits23b < 1 && !fPairs && fB != -1710)
        {
            det = fi23b;
            detector[countdet] = det;
            xdet[countdet] = 0.;
            ydet[countdet] = 0.;
            zdet[countdet] = 0.;
            qdet[countdet] = 8.;
            tdet[countdet] = 0.;

            fh_Fib_ToF_ac[det]->Fill(xdet[countdet] * 100., tStart - tdet[countdet]);
            fh_xy_Fib_ac[det]->Fill(xdet[countdet] * 100., ydet[countdet] * 100.);
            fh_mult_Fib_ac[det]->Fill(mult23a);
            fh_ToT_Fib_ac[det]->Fill(xdet[countdet] * 100., qdet[countdet]);
            fh_Fibs_vs_Tofd_ac[det]->Fill(x1[tofd1r] * 100. + randx, xdet[countdet] * 100.);
            fh_Fib_vs_Events_ac[det]->Fill(fNEvents, xdet[countdet] * 100.);
            fh_ToF_vs_Events_ac[det]->Fill(fNEvents, tStart - tdet[countdet]);
            fh_Fib_Time_ac[det]->Fill(xdet[countdet] * 100., tdet[countdet]);

            countdet++;
        }

        Bool_t cond1 = kFALSE;
        Bool_t cond2 = kFALSE;
        Bool_t cond3 = kFALSE;
        Bool_t cond4 = kFALSE;
        Bool_t cond5 = kFALSE;
        Bool_t cond6 = kFALSE;

        Bool_t tempFi33 = false;
        Bool_t tempFi32 = false;
        Bool_t tempFi23al = false;
        Bool_t tempFi23bl = false;
        Bool_t tempFi23ar = false;
        Bool_t tempFi23br = false;

        Nhits_before_fibcor =
            Nhits_before_fibcor + mult23a + mult23b + mult30 + mult32 + mult31 + mult33 + i1rc + i1lc + i2rc + i2lc;

        // FIBCORREL
        if (fibCuts && (mult30 + mult32) != 1 && (mult31 + mult33) != 1)
        //        if (fibCuts  && (mult30 + mult32) != 1 && mult23a > 0 && mult23b > 0 &&
        //            (mult31 + mult33) != 1 )
        //        if(fibCuts)
        {

            if (debug_fibcut)
                cout << "Entering fibcuts with multip fi23a/fi23b/fi30/fi31/fi32/fi33/tofd1r/tofd1l/tofd2r/tofd2l: "
                     << mult23a << ", " << mult23b << ", " << mult30 << ", " << mult31 << ", " << mult32 << ", "
                     << mult33 << ", " << i1rc << ", " << i1lc << ", " << i2rc << ", " << i2lc << endl;
            if (debug_fibcut)
                cout << "**** FibCuts ON THE RIGHT SIDE *****" << endl;

            Double_t offset_corr = 0.;
            if (abs(fB) == 2310. || abs(fB) == 2160. || abs(fB) == 1110)
                offset_corr = 3.;
            // Plots of correlations of Fiber detectors and register events for tracke
            if (debug_fibcut)
                cout << "**Loop over Tofd1r ****" << endl;

            countdet = 0;
            for (Int_t i1 = 0; i1 < i1rc; i1++)
            {

                if (qTofd1rs[i1] > 0 && !cTofd1rs[i1])
                {
                    detector[countdet] = tofd1rc;
                    xdet[countdet] = xTofd1rs[i1];
                    ydet[countdet] = yTofd1rs[i1];
                    zdet[countdet] = 0.;
                    tdet[countdet] = tTofd1rs[i1];
                    qdet[countdet] = (int)(qTofd1rs[i1] + 0.5); // q for tracker must be integer

                    if (debug_fibcut)
                        cout << "Tofd1r saved in fibcuts: " << xdet[countdet] * 100 << "; " << ydet[countdet] * 100
                             << "; " << qdet[countdet] << endl;

                    // fh_xy_tofd_ac->Fill(x1[tofd1r] * 100., y1[tofd1r] * 100.);

                    countdet++;
                    cTofd1rs[i1] = true;
                }

                for (Int_t i = 0; i < mult31; i++)
                {
                    fh_tofd_vs_Fib31->Fill(xFi31[i] * 100., xTofd1rs[i1] * 100.);
                    fh_ytofd_vs_yFib31->Fill(yFi31[i] * 100., yTofd1rs[i1] * 100.);
                    fh_qtofd_vs_qFib[fi31]->Fill(qFi31[i], qTofd1rs[i1]);

                    Double_t xtemp31;
                    xtemp31 = -0.732 * xTofd1rs[i1] * 100. - 29.364; // erste ToFD Ebene

                    if (debug_fibcut)
                    {
                        cout << "Fib31 before cut: " << i << " x: " << xFi31[i] * 100. << " q: " << qFi31[i] << endl;
                        cout << "tofd1r before fi31 cut: " << xTofd1rs[i1] * 100. << ", " << yTofd1rs[i1] * 100. << ", "
                             << qTofd1rs[i1] << endl;
                        cout << "diff: " << std::abs(xtemp31 - xFi31[i] * 100.) << ", " << dx3 << endl;
                    }

                    if (fPairs && std::abs(xtemp31 - xFi31[i] * 100.) > dx3)
                        continue;

                    if ((tTofd1rs[i1] - timeFi31[i]) > 30. || (tTofd1rs[i1] - timeFi31[i]) < -20.)
                        continue;

                    if (debug_fibcut)
                    {
                        cout << "Fib31 after cut: " << i << " x: " << xFi31[i] * 100. << " q: " << qFi31[i] << endl;
                        cout << "tofd1r after fi31 cut: " << xTofd1rs[i1] * 100. << ", " << yTofd1rs[i1] * 100. << ", "
                             << qTofd1rs[i1] << endl;
                    }

                    fh_tofd_vs_Fib31_ac->Fill(xFi31[i] * 100., xTofd1rs[i1] * 100.);
                    fh_Fib_Time_ac[fi31]->Fill(xFi31[i] * 100., timeFi31[i]);
                    fh_xy_Fib_ac[fi31]->Fill(xFi31[i] * 100., yFi31[i] * 100.);
                    fh_ToT_Fib_ac[fi31]->Fill(xFi31[i] * 100., qFi31[i]);
                    fh_Fibs_vs_Tofd_ac[fi31]->Fill(xTofd1rs[i1] * 100. + randx, xFi31[i] * 100.);
                    fh_Fib_ToF_ac[fi31]->Fill(xFi31[i] * 100., tTofd1rs[i1] - timeFi31[i]);
                    fh_Fib_ToF_vs_Qtofd[fi31]->Fill(qFi31[i], tTofd1rs[i1] - timeFi31[i]);
                    fh_qtofd_vs_qFib_ac[fi31]->Fill(qFi31[i], qTofd1rs[i1]);

                    for (Int_t j = 0; j < mult33; j++)
                    {
                        Double_t x31 =
                            0.866 * xFi33[j] * 100. - 4.637 - offset_corr; // 1.1165 * xFi31[i] * 100. + 2.53856;

                        fh_Fib33_vs_Fib31->Fill(xFi31[i] * 100., xFi33[j] * 100.);
                        fh_Fib33_vs_Fib31_dx->Fill(xFi31[i] * 100., xFi33[j] * 100. - xFi31[i] * 100.);
                        fh_Fib33_vs_Fib31_dt->Fill(xFi33[j] * 100., timeFi31[i] - timeFi33[j]);
                        fh_qtofd_vs_qFib[fi33]->Fill(qFi33[j], qTofd1rs[i1]);

                        tempFi33 = false;
                        // if (fSimu)
                        //	tempFi33 = true;
                        if (abs(timeFi31[i] - timeFi33[j]) < dtft && (tTofd1rs[i1] - timeFi33[j]) < 30 &&
                            (tTofd1rs[i1] - timeFi33[j]) > -25.)
                            tempFi33 = true;
                        if (abs(xFi31[i] * 100. - x31) < dx1 && xFi31[i] * 100. > -30. && xFi33[j] * 100. > -30. &&
                            tempFi33)
                        {
                            fh_Fib33_vs_Fib31_ac->Fill(xFi31[i] * 100., xFi33[j] * 100.);
                            fh_Fib33_vs_Fib31_dx_ac->Fill(xFi31[i] * 100., xFi33[j] * 100. - xFi31[i] * 100.);

                            fh_Fib_Time_ac[fi33]->Fill(xFi33[j] * 100., timeFi33[j]);
                            fh_xy_Fib_ac[fi33]->Fill(xFi33[j] * 100., yFi33[j] * 100.);
                            fh_ToT_Fib_ac[fi33]->Fill(xFi33[j] * 100., qFi33[j]);
                            fh_Fibs_vs_Tofd_ac[fi33]->Fill(xTofd1rs[i1] * 100. + randx, xFi33[j] * 100.);
                            fh_Fib_ToF_ac[fi33]->Fill(xFi33[j] * 100., tTofd1rs[i1] - timeFi33[j]);
                            fh_qtofd_vs_qFib_ac[fi33]->Fill(qFi33[j], qTofd1rs[i1]);
                            fh_Fib33_vs_Fib31_tot->Fill(qFi33[j], qFi31[i]);
                            fh_Fib_ToF_vs_Qtofd[fi33]->Fill(qFi33[j], tTofd1rs[i1] - timeFi33[j]);

                            if (debug_fibcut)
                                cout << "Fib33 after fi31 cuts: " << j << " x: " << xFi33[j] << " q: " << qFi33[j]
                                     << endl;

                            cond1 = kTRUE;
                            if (debug_fibcut)
                                cout << "cond1" << endl;

                            if (!maxWerte)
                            {
                                if (!fFi31[i])
                                {
                                    Double_t qtemp = qFi31[i];
                                    if (fGraphCuts && !fSimu)
                                    {
                                        if (cut_qfi31_qfi33->IsInside(qFi33[j], qFi31[i]))
                                            qtemp = 6.;
                                        else
                                            qtemp = 2.;
                                    }
                                    detector[countdet] = fi31;
                                    xdet[countdet] = xFi31[i];
                                    ydet[countdet] = yFi31[i];
                                    zdet[countdet] = 0.;
                                    qdet[countdet] = qtemp; // qFi31[i];
                                    tdet[countdet] = tTofd1rs[i1] - timeFi31[i];
                                    countdet++;
                                    fFi31[i] = true;
                                }
                                if (!fFi33[j])
                                {
                                    Double_t qtemp = qFi33[j];
                                    if (fGraphCuts && !fSimu)
                                    {
                                        if (cut_qfi31_qfi33->IsInside(qFi33[j], qFi31[i]))
                                            qtemp = 6.;
                                        else
                                            qtemp = 2.;
                                    }
                                    detector[countdet] = fi33;
                                    xdet[countdet] = xFi33[j];
                                    ydet[countdet] = yFi33[j];
                                    zdet[countdet] = 0.;
                                    qdet[countdet] = qtemp; // qFi33[j];
                                    tdet[countdet] = tTofd1rs[i1] - timeFi33[j];
                                    countdet++;
                                    fFi33[j] = true;
                                }
                            }
                        }
                    }
                    for (Int_t j = 0; j < mult23a; j++)
                    {
                        Double_t x31 = -6.625 * xFi23a[j] * 100. - 29.29; //-6.8672 * xFi23a[j] * 100. - 27.3507;
                        // if (abs(xFi31[i] * 100. - x31) < dx2 && xFi31[i] > -100 && xFi23a[j] > -100)
                        Double_t x23_ridge = -4.42618e+00 - 1.49472e-01 * xFi31[i] * 100;

                        fh_Fib31_vs_Fib23a->Fill(xFi23a[j] * 100., xFi31[i] * 100.);
                        fh_Fib31_vs_Fib23a_dx->Fill(xFi23a[j] * 100., xFi31[i] * 100. - xFi23a[j] * 100.);
                        fh_Fib31_vs_Fib23a_dt->Fill(xFi31[i] * 100. - xFi23a[j] * 100., timeFi31[i] - timeFi23a[j]);

                        tempFi23ar = false;
                        // if (fSimu)
                        // tempFi23ar = true;
                        if ((tTofd1rs[i1] - timeFi23a[j]) < 28 && (tTofd1rs[i1] - timeFi23a[j]) > -25.)
                            tempFi23ar = true;
                        // if (abs(xFi31[i] * 100. - x31) < dx2 && xFi31[i]*100. > -30 && xFi23a[j]*100. > -6 &&
                        if ((!fPairs) ||
                            //   (fPairs && (xFi23a[j] * 100.) < -0.26 && tempFi23ar)) // && abs(xFi31[i] * 100. - x31)
                            //   < dx2))
                            (fPairs && tempFi23ar && abs(xFi23a[j] * 100. - x23_ridge) < dx4))
                        {
                            // if (fGraphCuts && !cut_fi31_fi23a->IsInside(xFi23a[j] * 100., xFi31[i] * 100.))
                            //	continue;
                            fh_Fib31_vs_Fib23a_ac->Fill(xFi23a[j] * 100., xFi31[i] * 100.);
                            fh_Fib31_vs_Fib23a_dx_ac->Fill(xFi23a[j] * 100., xFi31[i] * 100. - xFi23a[j] * 100.);

                            fh_Fib_Time_ac[fi23a]->Fill(xFi23a[j] * 100., timeFi23a[j]);
                            fh_xy_Fib_ac[fi23a]->Fill(xFi23a[j] * 100., yFi23a[j] * 100.);
                            fh_ToT_Fib_ac[fi23a]->Fill(xFi23a[j] * 100., qFi23a[j]);
                            fh_Fibs_vs_Tofd_ac[fi23a]->Fill(xTofd1rs[i1] * 100., xFi23a[j] * 100.);
                            fh_Fib_ToF_ac[fi23a]->Fill(xFi23a[j] * 100., tTofd1rs[i1] - timeFi23a[j]);
                            fh_Fib23a_vs_Fib3X_tot->Fill(qFi23a[j], qFi31[i]);
                            fh_Fib_ToF_vs_Qtofd[fi23a]->Fill(qTofd1rs[i1], tTofd1rs[i1] - timeFi23a[j]);

                            if (debug_fibcut)
                                cout << "Fib23a after fi31 cuts" << j << " x: " << xFi23a[j] << " q: " << qFi23a[j]
                                     << endl;

                            cond2 = kTRUE;
                            if (debug_fibcut)
                                cout << "cond2" << endl;

                            if (!maxWerte)
                            {
                                if (!fFi23a[j])
                                {
                                    detector[countdet] = fi23a;
                                    xdet[countdet] = xFi23a[j];
                                    ydet[countdet] = yFi23a[j];
                                    zdet[countdet] = 0.;
                                    qdet[countdet] = qFi23a[j];
                                    tdet[countdet] = tTofd1rs[i1] - timeFi23a[j];
                                    countdet++;

                                    detector[countdet] = fi23ar;
                                    xdet[countdet] = xFi23a[j];
                                    ydet[countdet] = yFi23a[j];
                                    zdet[countdet] = 0.;
                                    qdet[countdet] = qFi23a[j];
                                    tdet[countdet] = tTofd1rs[i1] - timeFi23a[j];
                                    countdet++;

                                    fFi23a[j] = true;
                                }
                            }
                        }
                    }
                    for (Int_t j = 0; j < mult23b; j++)
                    {
                        // Double_t x31 = -6.8672 * xFi23b[j] * 100. - 27.3507;
                        // if (abs(xFi31[i] * 100. - x31) < dx2 && xFi31[i] > -100 && xFi23b[j] > -100)
                        fh_tofdright_vs_Fib23b->Fill(yFi23b[j] * 100., yTofd1rs[i1] * 100.);
                        fh_Fib31_vs_Fib23b_dx->Fill(yFi23b[j] * 100., yFi31[i] * 100. - yFi23b[j] * 100.);
                        fh_Fib31_vs_Fib23b_dt->Fill(yFi23b[j] * 100., timeFi31[i] - timeFi23b[j]);

                        tempFi23br = false;
                        //	if (fSimu)
                        //	tempFi23br = true;
                        if ((tTofd1rs[i1] - timeFi23b[j]) < 25. && (tTofd1rs[i1] - timeFi23b[j]) > -25.)
                            tempFi23br = true;

                        Double_t ycheck = 7.6358668 * yFi23b[j] * 100. - 0.0817257;
                        if (debug_fibcut)
                        {
                            cout << "Before Tofd vs fi23b: " << endl;
                            cout << "Tofd1r: " << yTofd1rs[i1] * 100 << "; " << qTofd1rs[i1] << ", " << endl;
                            cout << "Fib23b: " << yFi23b[j] * 100. << endl;
                            cout << "ycheck1:" << ycheck << ", " << abs(yTofd1rs[i1] * 100. - ycheck) << endl;
                        }

                        if (fPairs && std::abs(yTofd1rs[i1] * 100. - ycheck) > dytofdfib23b)
                            continue;

                        if ((!fPairs) ||
                            (fPairs && tempFi23br)) //(fPairs && abs(yFi23b[j] * 100.) > 0.26 && tempFi23br))
                        {
                            // if (fGraphCuts && !cut_fi30_fi23b->IsInside(yFi23b[j] * 100., xFi30[i] * 100.))
                            //	continue;
                            fh_tofdright_vs_Fib23b_ac->Fill(yFi23b[j] * 100., yTofd1rs[i1] * 100.);
                            fh_Fib31_vs_Fib23b_dx_ac->Fill(yFi23b[j] * 100., yFi31[i] * 100. - yFi23b[j] * 100.);

                            fh_Fib_Time_ac[fi23b]->Fill(yFi23b[j] * 100., timeFi23b[j]);
                            fh_xy_Fib_ac[fi23b]->Fill(xFi23b[j] * 100., yFi23b[j] * 100.);
                            fh_ToT_Fib_ac[fi23b]->Fill(yFi23b[j] * 100., qFi23b[j]);
                            fh_Fibs_vs_Tofd_ac[fi23b]->Fill(yTofd1rs[i1] * 100., yFi23b[j] * 100.);
                            fh_Fib_ToF_ac[fi23b]->Fill(yFi23b[j] * 100., tTofd1rs[i1] - timeFi23b[j]);
                            fh_Fib23b_vs_Fib3X_tot->Fill(qFi23b[j], qFi31[i]);
                            fh_Fib_ToF_vs_Qtofd[fi23b]->Fill(qTofd1rs[i1], tTofd1rs[i1] - timeFi23b[j]);

                            if (debug_fibcut)
                                cout << "Fib23b after fi31 cuts: " << j << " y: " << yFi23b[j] << " q: " << qFi23b[j]
                                     << endl;

                            cond3 = kTRUE;
                            if (debug_fibcut)
                                cout << "cond3" << endl;

                            if (!maxWerte)
                            {
                                if (!fFi23b[j])
                                {
                                    detector[countdet] = fi23b;
                                    xdet[countdet] = xFi23b[j];
                                    ydet[countdet] = yFi23b[j];
                                    zdet[countdet] = 0.;
                                    qdet[countdet] = qFi23b[j];
                                    tdet[countdet] = tTofd1rs[i1] - timeFi23b[j];
                                    countdet++;
                                    fFi23b[j] = true;
                                }
                            }
                        }
                    }
                }
            }

            if (debug_fibcut)
                cout << "**Loop over Tofd2r ****" << endl;

            for (Int_t i1 = 0; i1 < i2rc; i1++)
            {

                if (qTofd2rs[i1] > 0 && !cTofd2rs[i1])
                {
                    detector[countdet] = tofd2rc;
                    xdet[countdet] = xTofd2rs[i1];
                    ydet[countdet] = yTofd2rs[i1];
                    zdet[countdet] = 0.;
                    tdet[countdet] = tTofd2rs[i1];
                    qdet[countdet] = (int)(qTofd2rs[i1] + 0.5); // q for tracker must be integer

                    if (debug_fibcut)
                        cout << "Tofd2r saved in fibcuts: " << xdet[countdet] * 100 << "; " << ydet[countdet] * 100
                             << "; " << qdet[countdet] << endl;

                    // fh_xy_tofd_ac->Fill(x1[tofd1r] * 100., y1[tofd1r] * 100.);

                    countdet++;
                    cTofd2rs[i1] = true;
                }

                for (Int_t i = 0; i < mult31; i++)
                {
                    fh_tofd_vs_Fib31->Fill(xFi31[i] * 100., xTofd2rs[i1] * 100.);
                    fh_ytofd_vs_yFib31->Fill(yFi31[i] * 100., yTofd2rs[i1] * 100.);
                    fh_qtofd_vs_qFib[fi31]->Fill(qFi31[i], qTofd2rs[i1]);

                    Double_t xtemp31;
                    xtemp31 = -0.732 * xTofd2rs[i1] * 100. - 29.364; // erste ToFD Ebene

                    if (debug_fibcut)
                    {
                        cout << "Fib31 before cut: " << i << " x: " << xFi31[i] * 100. << " q: " << qFi31[i] << endl;
                        cout << "tofd2r before fi31 cut: " << xTofd2rs[i1] * 100. << ", " << yTofd2rs[i1] * 100. << ", "
                             << qTofd2rs[i1] << endl;
                        cout << "diff: " << std::abs(xtemp31 - xFi31[i] * 100.) << ", " << dx3 << endl;
                    }

                    if (fPairs && std::abs(xtemp31 - xFi31[i] * 100.) > dx3)
                        continue;
                    if ((tTofd2rs[i1] - timeFi31[i]) > 30. || tTofd2rs[i1] - timeFi31[i] < -20.)
                        continue;

                    if (debug_fibcut)
                    {
                        cout << "Fib31 after cut: " << i << " x: " << xFi31[i] * 100. << " q: " << qFi31[i] << endl;
                        cout << "tofd2r after fi31 cut: " << xTofd2rs[i1] * 100. << ", " << yTofd2rs[i1] * 100. << ", "
                             << qTofd2rs[i1] << endl;
                    }

                    fh_tofd_vs_Fib31_ac->Fill(xFi31[i] * 100., xTofd2rs[i1] * 100.);
                    fh_Fib_Time_ac[fi31]->Fill(xFi31[i] * 100., timeFi31[i]);
                    fh_xy_Fib_ac[fi31]->Fill(xFi31[i] * 100., yFi31[i] * 100.);
                    fh_ToT_Fib_ac[fi31]->Fill(xFi31[i] * 100., qFi31[i]);
                    fh_Fibs_vs_Tofd_ac[fi31]->Fill(xTofd2rs[i1] * 100. + randx, xFi31[i] * 100.);
                    fh_Fib_ToF_ac[fi31]->Fill(xFi31[i] * 100., tTofd2rs[i1] - timeFi31[i]);
                    fh_qtofd_vs_qFib_ac[fi31]->Fill(qFi31[i], qTofd2rs[i1]);
                    fh_Fib_ToF_vs_Qtofd[fi31]->Fill(qFi31[i], tTofd2rs[i1] - timeFi31[i]);

                    for (Int_t j = 0; j < mult33; j++)
                    {
                        Double_t x31 =
                            0.866 * xFi33[j] * 100. - 4.637 - offset_corr; // 1.1165 * xFi31[i] * 100. + 2.53856;

                        fh_Fib33_vs_Fib31->Fill(xFi31[i] * 100., xFi33[j] * 100.);
                        fh_Fib33_vs_Fib31_dx->Fill(xFi31[i] * 100., xFi33[j] * 100. - xFi31[i] * 100.);
                        fh_Fib33_vs_Fib31_dt->Fill(xFi33[j] * 100., timeFi31[i] - timeFi33[j]);
                        fh_qtofd_vs_qFib[fi33]->Fill(qFi33[j], qTofd2rs[i1]);

                        tempFi33 = false;
                        // if (fSimu)
                        //   tempFi33 = true;
                        if (abs(timeFi31[i] - timeFi33[j]) < dtft && (tTofd2rs[i1] - timeFi33[j]) < 30 &&
                            (tTofd2rs[i1] - timeFi33[j]) > -25.)
                            tempFi33 = true;
                        if (abs(xFi31[i] * 100. - x31) < dx1 && xFi31[i] * 100. > -30. && xFi33[j] * 100. > -30. &&
                            tempFi33)
                        {
                            fh_Fib33_vs_Fib31_ac->Fill(xFi31[i] * 100., xFi33[j] * 100.);
                            fh_Fib33_vs_Fib31_dx_ac->Fill(xFi31[i] * 100., xFi33[j] * 100. - xFi31[i] * 100.);

                            fh_Fib_Time_ac[fi33]->Fill(xFi33[j] * 100., timeFi33[j]);
                            fh_xy_Fib_ac[fi33]->Fill(xFi33[j] * 100., yFi33[j] * 100.);
                            fh_ToT_Fib_ac[fi33]->Fill(xFi33[j] * 100., qFi33[j]);
                            fh_Fibs_vs_Tofd_ac[fi33]->Fill(xTofd2rs[i1] * 100. + randx, xFi33[j] * 100.);
                            fh_Fib_ToF_ac[fi33]->Fill(xFi33[j] * 100., tTofd2rs[i1] - timeFi33[j]);
                            fh_qtofd_vs_qFib_ac[fi33]->Fill(qFi33[j], qTofd2rs[i1]);
                            fh_Fib33_vs_Fib31_tot->Fill(qFi33[j], qFi31[i]);
                            fh_Fib_ToF_vs_Qtofd[fi33]->Fill(qFi33[j], tTofd2rs[i1] - timeFi33[j]);

                            if (debug_fibcut)
                                cout << "Fib33 after fi31 cuts: " << j << " x: " << xFi33[j] << " q: " << qFi33[j]
                                     << endl;

                            cond1 = kTRUE;
                            if (debug_fibcut)
                                cout << "cond1" << endl;

                            if (!maxWerte)
                            {
                                if (!fFi31[i])
                                {
                                    Double_t qtemp = qFi31[i];
                                    if (fGraphCuts && !fSimu)
                                    {
                                        if (cut_qfi31_qfi33->IsInside(qFi33[j], qFi31[i]))
                                            qtemp = 6.;
                                        else
                                            qtemp = 2.;
                                    }
                                    detector[countdet] = fi31;
                                    xdet[countdet] = xFi31[i];
                                    ydet[countdet] = yFi31[i];
                                    zdet[countdet] = 0.;
                                    qdet[countdet] = qtemp; // qFi31[i];
                                    tdet[countdet] = tTofd2rs[i1] - timeFi31[i];
                                    countdet++;
                                    fFi31[i] = true;
                                }
                                if (!fFi33[j])
                                {
                                    Double_t qtemp = qFi33[j];
                                    if (fGraphCuts && !fSimu)
                                    {
                                        if (cut_qfi31_qfi33->IsInside(qFi33[j], qFi31[i]))
                                            qtemp = 6.;
                                        else
                                            qtemp = 2.;
                                    }
                                    detector[countdet] = fi33;
                                    xdet[countdet] = xFi33[j];
                                    ydet[countdet] = yFi33[j];
                                    zdet[countdet] = 0.;
                                    qdet[countdet] = qtemp; // qFi33[j];
                                    tdet[countdet] = tTofd2rs[i1] - timeFi33[j];
                                    countdet++;
                                    fFi33[j] = true;
                                }
                            }
                        }
                    }
                    for (Int_t j = 0; j < mult23a; j++)
                    {
                        Double_t x31 = -6.625 * xFi23a[j] * 100. - 29.29; //-6.8672 * xFi23a[j] * 100. - 27.3507;
                        // if (abs(xFi31[i] * 100. - x31) < dx2 && xFi31[i] > -100 && xFi23a[j] > -100)
                        Double_t x23_ridge = -4.42618e+00 - 1.49472e-01 * xFi31[i] * 100;

                        fh_Fib31_vs_Fib23a->Fill(xFi23a[j] * 100., xFi31[i] * 100.);
                        fh_Fib31_vs_Fib23a_dx->Fill(xFi23a[j] * 100., xFi31[i] * 100. - xFi23a[j] * 100.);
                        fh_Fib31_vs_Fib23a_dt->Fill(xFi31[i] * 100. - xFi23a[j] * 100., timeFi31[i] - timeFi23a[j]);

                        tempFi23ar = false;
                        //	if (fSimu)
                        //	tempFi23ar = true;
                        if ((tTofd2rs[i1] - timeFi23a[j]) < 28. && (tTofd2rs[i1] - timeFi23a[j]) > -25.)
                            tempFi23ar = true;
                        // if (abs(xFi31[i] * 100. - x31) < dx2 && xFi31[i]*100. > -30 && xFi23a[j]*100. > -6 &&
                        if ((!fPairs) ||
                            //   (fPairs && (xFi23a[j] * 100.) < -0.26 && tempFi23ar)) // && abs(xFi31[i] * 100. - x31)
                            //   < dx2))
                            (fPairs && tempFi23ar && abs(xFi23a[j] * 100. - x23_ridge) < dx4))
                        {

                            // if (fGraphCuts && !cut_fi31_fi23a->IsInside(xFi23a[j] * 100., xFi31[i] * 100.))
                            //	continue;
                            fh_Fib31_vs_Fib23a_ac->Fill(xFi23a[j] * 100., xFi31[i] * 100.);
                            fh_Fib31_vs_Fib23a_dx_ac->Fill(xFi23a[j] * 100., xFi31[i] * 100. - xFi23a[j] * 100.);

                            fh_Fib_Time_ac[fi23a]->Fill(xFi23a[j] * 100., timeFi23a[j]);
                            fh_xy_Fib_ac[fi23a]->Fill(xFi23a[j] * 100., yFi23a[j] * 100.);
                            fh_ToT_Fib_ac[fi23a]->Fill(xFi23a[j] * 100., qFi23a[j]);
                            fh_Fibs_vs_Tofd_ac[fi23a]->Fill(xTofd2rs[i1] * 100., xFi23a[j] * 100.);
                            fh_Fib_ToF_ac[fi23a]->Fill(xFi23a[j] * 100., tTofd2rs[i1] - timeFi23a[j]);
                            fh_Fib23a_vs_Fib3X_tot->Fill(qFi23a[j], qFi31[i]);
                            fh_Fib_ToF_vs_Qtofd[fi23a]->Fill(qTofd2rs[i1], tTofd2rs[i1] - timeFi23a[j]);

                            if (debug_fibcut)
                                cout << "Fib23a after fi31 cuts" << j << " x: " << xFi23a[j] << " q: " << qFi23a[j]
                                     << endl;

                            cond2 = kTRUE;
                            if (debug_fibcut)
                                cout << "cond2" << endl;

                            if (!maxWerte)
                            {
                                if (!fFi23a[j])
                                {
                                    detector[countdet] = fi23a;
                                    xdet[countdet] = xFi23a[j];
                                    ydet[countdet] = yFi23a[j];
                                    zdet[countdet] = 0.;
                                    qdet[countdet] = qFi23a[j];
                                    tdet[countdet] = tTofd2rs[i1] - timeFi23a[j];
                                    countdet++;

                                    detector[countdet] = fi23ar;
                                    xdet[countdet] = xFi23a[j];
                                    ydet[countdet] = yFi23a[j];
                                    zdet[countdet] = 0.;
                                    qdet[countdet] = qFi23a[j];
                                    tdet[countdet] = tTofd2rs[i1] - timeFi23a[j];
                                    countdet++;

                                    fFi23a[j] = true;
                                }
                            }
                        }
                    }
                    for (Int_t j = 0; j < mult23b; j++)
                    {
                        // Double_t x31 = -6.8672 * xFi23b[j] * 100. - 27.3507;
                        // if (abs(xFi31[i] * 100. - x31) < dx2 && xFi31[i] > -100 && xFi23b[j] > -100)
                        fh_tofdright_vs_Fib23b->Fill(yFi23b[j] * 100., yTofd2rs[i1] * 100.);
                        fh_Fib31_vs_Fib23b_dx->Fill(yFi23b[j] * 100., yFi31[i] * 100. - yFi23b[j] * 100.);
                        fh_Fib31_vs_Fib23b_dt->Fill(yFi23b[j] * 100., timeFi31[i] - timeFi23b[j]);

                        tempFi23br = false;
                        //	if (fSimu)
                        //		tempFi23br = true;
                        if ((tTofd2rs[i1] - timeFi23b[j]) < 25. && (tTofd2rs[i1] - timeFi23b[j]) > -25.)
                            tempFi23br = true;

                        Double_t ycheck = 7.6358668 * yFi23b[j] * 100. - 0.0817257;
                        if (debug_fibcut)
                        {
                            cout << "Before Tofd vs fi23b: " << endl;
                            cout << "Tofd2r: " << yTofd2rs[i1] * 100 << "; " << qTofd2rs[i1] << ", " << endl;
                            cout << "Fib23b: " << yFi23b[j] * 100. << endl;
                            cout << "ycheck1:" << ycheck << ", " << abs(yTofd2rs[i1] * 100. - ycheck) << endl;
                        }

                        if (fPairs && std::abs(yTofd2rs[i1] * 100. - ycheck) > dytofdfib23b)
                            continue;

                        if ((!fPairs) ||
                            (fPairs && tempFi23br)) //(fPairs && abs(yFi23b[j] * 100.) > 0.26 && tempFi23br))
                        {
                            // if (fGraphCuts && !cut_fi30_fi23b->IsInside(yFi23b[j] * 100., xFi30[i] * 100.))
                            //	continue;
                            fh_tofdright_vs_Fib23b_ac->Fill(yFi23b[j] * 100., yTofd2rs[i1] * 100.);
                            fh_Fib31_vs_Fib23b_dx_ac->Fill(yFi23b[j] * 100., yFi31[i] * 100. - yFi23b[j] * 100.);

                            fh_Fib_Time_ac[fi23b]->Fill(yFi23b[j] * 100., timeFi23b[j]);
                            fh_xy_Fib_ac[fi23b]->Fill(xFi23b[j] * 100., yFi23b[j] * 100.);
                            fh_ToT_Fib_ac[fi23b]->Fill(yFi23b[j] * 100., qFi23b[j]);
                            fh_Fibs_vs_Tofd_ac[fi23b]->Fill(yTofd2rs[i1] * 100., yFi23b[j] * 100.);
                            fh_Fib_ToF_ac[fi23b]->Fill(yFi23b[j] * 100., tTofd2rs[i1] - timeFi23b[j]);
                            fh_Fib23b_vs_Fib3X_tot->Fill(qFi23b[j], qFi31[i]);
                            fh_Fib_ToF_vs_Qtofd[fi23b]->Fill(qTofd2rs[i1], tTofd2rs[i1] - timeFi23b[j]);

                            if (debug_fibcut)
                                cout << "Fib23b after fi31 cuts: " << j << " y: " << yFi23b[j] << " q: " << qFi23b[j]
                                     << endl;

                            cond3 = kTRUE;
                            if (debug_fibcut)
                                cout << "cond3" << endl;

                            if (!maxWerte)
                            {
                                if (!fFi23b[j])
                                {
                                    detector[countdet] = fi23b;
                                    xdet[countdet] = xFi23b[j];
                                    ydet[countdet] = yFi23b[j];
                                    zdet[countdet] = 0.;
                                    qdet[countdet] = qFi23b[j];
                                    tdet[countdet] = tTofd2rs[i1] - timeFi23b[j];
                                    countdet++;
                                    fFi23b[j] = true;
                                }
                            }
                        }
                    }
                }
            }

            if (debug_fibcut)
                cout << "**Loop over Tofd3r (Tofi) ****" << endl;

            for (Int_t i1 = 0; i1 < i3rc; i1++)
            {

                if (qTofd3rs[i1] > 0 && !cTofd3rs[i1])
                {
                    detector[countdet] = tofd3rc;
                    xdet[countdet] = xTofd3rs[i1];
                    ydet[countdet] = yTofd3rs[i1];
                    zdet[countdet] = 0.;
                    tdet[countdet] = tTofd3rs[i1];
                    qdet[countdet] = (int)(qTofd3rs[i1] + 0.5); // q for tracker must be integer

                    if (debug_fibcut)
                        cout << "Tofd3r saved in fibcuts: " << xdet[countdet] * 100 << "; " << ydet[countdet] * 100
                             << "; " << qdet[countdet] << endl;

                    // fh_xy_tofd_ac->Fill(x1[tofd3r] * 100., y1[tofd3r] * 100.);

                    countdet++;
                    cTofd3rs[i1] = true;
                }
            }
            if (debug_fibcut)
                cout << "**** FibCuts ON THE LEFT SIDE *****" << endl;

            if (debug_fibcut)
                cout << "**Loop over Tofd1l ****" << endl;

            for (Int_t i1 = 0; i1 < i1lc; i1++)
            {
                if (qTofd1ls[i1] > 0 && !cTofd1ls[i1])
                {
                    detector[countdet] = tofd1lc;
                    xdet[countdet] = xTofd1ls[i1];
                    ydet[countdet] = yTofd1ls[i1];
                    zdet[countdet] = 0.;
                    tdet[countdet] = tTofd1ls[i1];
                    qdet[countdet] = (int)(qTofd1ls[i1] + 0.5); // q for tracker must be integer

                    if (debug_fibcut)
                        cout << "Tofd1l saved in fibcuts: " << xdet[countdet] * 100 << "; " << ydet[countdet] * 100
                             << "; " << qdet[countdet] << endl;

                    //  fh_xy_tofd_ac->Fill(x1[tofd1l] * 100., y1[tofd1l] * 100.);

                    countdet++;
                    cTofd1ls[i1] = true;
                }

                for (Int_t i = 0; i < mult30; i++)
                {
                    fh_tofd_vs_Fib30->Fill(xFi30[i] * 100., xTofd1ls[i1] * 100.);
                    fh_ytofd_vs_yFib30->Fill(yFi30[i] * 100., yTofd1ls[i1] * 100.);
                    fh_qtofd_vs_qFib[fi30]->Fill(qFi30[i], qTofd1ls[i1]);

                    Double_t xtemp30;
                    xtemp30 = 0.7796 * xTofd1ls[i1] * 100. - 28.669; // erste ToFD Ebene

                    if (debug_fibcut)
                    {
                        cout << "Fib30 before cut: " << i << " x: " << xFi30[i] * 100. << " q: " << qFi30[i] << endl;
                        cout << "tofd1l before fi30 cut: " << xTofd1ls[i1] * 100. << ", " << yTofd1ls[i1] * 100. << ", "
                             << qTofd1ls[i1] << endl;
                        cout << "Diff: " << std::abs(xtemp30 - xFi30[i] * 100.) << ", " << dx3 << endl;
                    }

                    if (fPairs && std::abs(xtemp30 - xFi30[i] * 100.) > dx3)
                        continue;
                    if ((tTofd1ls[i1] - timeFi30[i]) > 20. || (tTofd1ls[i1] - timeFi30[i]) < -25.)
                        continue;

                    if (debug_fibcut)
                    {
                        cout << "Fib30after fibcuts: " << i << " x: " << xFi30[i] << " q: " << qFi30[i] << endl;
                        cout << "Tofd1l after fi30 cuts: " << xTofd1ls[i1] * 100. << ", " << yTofd1ls[i1] * 100. << ", "
                             << qTofd1ls[i1] << endl;
                    }

                    fh_tofd_vs_Fib30_ac->Fill(xFi30[i] * 100., xTofd1ls[i1] * 100.);
                    fh_Fib_ToF_ac[fi30]->Fill(xFi30[i] * 100., tTofd1ls[i1] - timeFi30[i]);
                    fh_xy_Fib_ac[fi30]->Fill(xFi30[i] * 100., yFi30[i] * 100.);
                    fh_ToT_Fib_ac[fi30]->Fill(xFi30[i] * 100., qFi30[i]);
                    fh_Fibs_vs_Tofd_ac[fi30]->Fill(xTofd1ls[i1] * 100. + randx, xFi30[i] * 100.);
                    fh_Fib_Time_ac[fi30]->Fill(xFi30[i] * 100., timeFi30[i]);
                    fh_qtofd_vs_qFib_ac[fi30]->Fill(qFi30[i], qTofd1ls[i1]);
                    fh_Fib_ToF_vs_Qtofd[fi30]->Fill(qTofd1ls[i1], tTofd1ls[i1] - timeFi30[i]);

                    for (Int_t j = 0; j < mult32; j++)
                    {
                        fh_Fib32_vs_Fib30->Fill(xFi30[i] * 100., xFi32[j] * 100.);
                        fh_Fib32_vs_Fib30_dx->Fill(xFi30[i] * 100., xFi32[j] * 100. - xFi30[i] * 100.);
                        fh_Fib32_vs_Fib30_dt->Fill(xFi32[j] * 100., timeFi30[i] - timeFi32[j]);
                        fh_qtofd_vs_qFib[fi32]->Fill(qFi32[j], qTofd1ls[i1]);

                        Double_t x30 =
                            0.873 * xFi32[j] * 100. - 3.302 - offset_corr; // 1.10926 * xFi30[i] * 100. + 2.8943;
                        tempFi32 = false;
                        // if (fSimu)
                        // tempFi32 = true;
                        if (abs(timeFi30[i] - timeFi32[j]) < dtft && (tTofd1ls[i1] - timeFi32[j]) < 20. &&
                            (tTofd1ls[i1] - timeFi32[j]) > -25.)
                            tempFi32 = true;

                        if ((abs(xFi30[i] * 100. - x30) < dx1 && xFi32[j] * 100. > -30 && xFi30[i] * 100. > -30. &&
                             tempFi32))
                        {
                            fh_Fib32_vs_Fib30_ac->Fill(xFi30[i] * 100., xFi32[j] * 100.);
                            fh_Fib32_vs_Fib30_dx_ac->Fill(xFi30[i] * 100., xFi32[j] * 100. - xFi30[i] * 100.);

                            fh_Fib_Time_ac[fi32]->Fill(xFi32[j] * 100., timeFi32[j]);
                            fh_xy_Fib_ac[fi32]->Fill(xFi32[j] * 100., yFi32[j] * 100.);
                            fh_ToT_Fib_ac[fi32]->Fill(xFi32[j] * 100., qFi32[j]);
                            fh_Fibs_vs_Tofd_ac[fi32]->Fill(xTofd1ls[i1] * 100. + randx, xFi32[j] * 100.);
                            fh_Fib_ToF_ac[fi32]->Fill(xFi32[j] * 100., tTofd1ls[i1] - timeFi32[j]);
                            fh_qtofd_vs_qFib_ac[fi32]->Fill(qFi32[j], qTofd1ls[i1]);
                            fh_Fib32_vs_Fib30_tot->Fill(qFi32[j], qFi30[i]);
                            fh_Fib_ToF_vs_Qtofd[fi32]->Fill(qTofd1ls[i1], tTofd1ls[i1] - timeFi32[j]);

                            if (debug_fibcut)
                                cout << "Fib32 after fi30 cuts: " << j << " x: " << xFi32[j] << " q: " << qFi32[j]
                                     << endl;

                            cond4 = kTRUE;
                            if (debug_fibcut)
                                cout << "cond4" << endl;

                            if (!maxWerte)
                            {
                                if (!fFi30[i])
                                {
                                    Double_t qtemp = qFi30[i];
                                    if (fGraphCuts && !fSimu)
                                    {
                                        if (cut_qfi30_qfi32->IsInside(qFi32[j], qFi30[i]))
                                            qtemp = 6.;
                                        else
                                            qtemp = 2;
                                    }
                                    detector[countdet] = fi30;
                                    xdet[countdet] = xFi30[i];
                                    ydet[countdet] = yFi30[i];
                                    zdet[countdet] = 0.;
                                    qdet[countdet] = qtemp; // qFi30[i];
                                    tdet[countdet] = tTofd1ls[i1] - timeFi30[i];
                                    countdet++;
                                    fFi30[i] = true;
                                }
                                if (!fFi32[j])
                                {
                                    Double_t qtemp = qFi32[j];
                                    if (fGraphCuts && !fSimu)
                                    {
                                        if (cut_qfi30_qfi32->IsInside(qFi32[j], qFi30[i]))
                                            qtemp = 6.;
                                        else
                                            qtemp = 2;
                                    }
                                    detector[countdet] = fi32;
                                    xdet[countdet] = xFi32[j];
                                    ydet[countdet] = yFi32[j];
                                    zdet[countdet] = 0.;
                                    qdet[countdet] = qtemp; // qFi32[j];
                                    tdet[countdet] = tTofd1ls[i1] - timeFi32[j];
                                    countdet++;
                                    fFi32[j] = true;
                                }
                            }
                        }
                    }
                    for (Int_t j = 0; j < mult23a; j++)
                    {
                        fh_Fib30_vs_Fib23a->Fill(xFi23a[j] * 100., xFi30[i] * 100.);
                        fh_Fib30_vs_Fib23a_dx->Fill(xFi23a[j] * 100., xFi30[i] * 100. - xFi23a[j] * 100.);
                        fh_Fib30_vs_Fib23a_dt->Fill(xFi30[i] * 100. - xFi23a[j] * 100., timeFi30[i] - timeFi23a[j]);

                        Double_t x30 = 6.498 * xFi23a[j] * 100. - 28.93;
                        // if (abs(xFi30[i] * 100. - x30) < dx2 && xFi30[i] > -100 && xFi23a[j] > -100)
                        Double_t x23_ridge = 4.27564e+00 + 1.48326e-01 * xFi30[i] * 100;

                        tempFi23al = false;
                        // if (fSimu)
                        // tempFi23al = true;
                        if ((tTofd1ls[i1] - timeFi23a[j]) < 28. && (tTofd1ls[i1] - timeFi23a[j]) > -25.)
                            tempFi23al = true;

                        if ((!fPairs) ||
                            //  (fPairs && (xFi23a[j] * 100.) > 0.26 && tempFi23al)) // && abs(xFi30[i] * 100. - x30) <
                            //  dx2))
                            (fPairs && tempFi23al && abs(xFi23a[j] * 100. - x23_ridge) < dx4))
                        {
                            // if (fGraphCuts && !cut_fi30_fi23b->IsInside(yFi23b[j] * 100., xFi30[i] * 100.))
                            //	continue;

                            fh_Fib30_vs_Fib23a_ac->Fill(xFi23a[j] * 100., xFi30[i] * 100.);
                            fh_Fib30_vs_Fib23a_dx_ac->Fill(xFi23a[j] * 100., xFi30[i] * 100. - xFi23a[j] * 100.);

                            fh_Fib_Time_ac[fi23a]->Fill(xFi23a[j] * 100., timeFi23a[j]);
                            fh_xy_Fib_ac[fi23a]->Fill(xFi23a[j] * 100., yFi23a[j] * 100.);
                            fh_ToT_Fib_ac[fi23a]->Fill(xFi23a[j] * 100., qFi23a[j]);
                            fh_Fibs_vs_Tofd_ac[fi23a]->Fill(xTofd1ls[i1] * 100., xFi23a[j] * 100.);
                            fh_Fib_ToF_ac[fi23a]->Fill(xFi23a[j] * 100., tTofd1ls[i1] - timeFi23a[j]);
                            fh_Fib23a_vs_Fib3X_tot->Fill(qFi23a[j], qFi30[i]);
                            fh_Fib_ToF_vs_Qtofd[fi23a]->Fill(qTofd1ls[i1], tTofd1ls[i1] - timeFi23a[j]);

                            if (debug_fibcut)
                                cout << "Fib23a after fi30 cuts: " << j << " x: " << xFi23a[j] << " q: " << qFi23a[j]
                                     << endl;

                            cond5 = kTRUE;
                            if (debug_fibcut)
                                cout << "cond5" << endl;

                            if (!maxWerte)
                            {
                                if (!fFi23a[j])
                                {
                                    detector[countdet] = fi23a;
                                    xdet[countdet] = xFi23a[j];
                                    ydet[countdet] = yFi23a[j];
                                    zdet[countdet] = 0.;
                                    qdet[countdet] = qFi23a[j];
                                    tdet[countdet] = tTofd1ls[i1] - timeFi23a[j];
                                    countdet++;

                                    detector[countdet] = fi23al;
                                    xdet[countdet] = xFi23a[j];
                                    ydet[countdet] = yFi23a[j];
                                    zdet[countdet] = 0.;
                                    qdet[countdet] = qFi23a[j];
                                    tdet[countdet] = tTofd1ls[i1] - timeFi23a[j];
                                    countdet++;
                                    fFi23a[j] = true;
                                }
                            }
                        }
                    }
                    for (Int_t j = 0; j < mult23b; j++)
                    {
                        // Double_t x30 = 6.98386 * xFi23b[j] * 100. - 27.39897;
                        // if (abs(xFi30[i] * 100. - x30) < dx2 && xFi30[i] > -100 && xFi23b[j] > -100)

                        fh_tofdleft_vs_Fib23b->Fill(yFi23b[j] * 100., yTofd1ls[i1] * 100.);
                        fh_Fib30_vs_Fib23b_dx->Fill(yFi23b[j] * 100., yFi30[i] * 100. - yFi23b[j] * 100.);
                        fh_Fib30_vs_Fib23b_dt->Fill(yFi23b[j] * 100., timeFi30[i] - timeFi23b[j]);

                        tempFi23bl = false;
                        // if (fSimu)
                        //  tempFi23bl = true;
                        if ((tTofd1ls[i1] - timeFi23b[j]) < 25. && (tTofd1ls[i1] - timeFi23b[j]) > -25.)
                            tempFi23bl = true;

                        Double_t ycheck = 7.6358668 * yFi23b[j] * 100. - 0.0817257;

                        if (debug_fibcut)
                        {
                            cout << "Before Tofd vs fi23b: " << endl;
                            cout << "Tofd1l: " << yTofd1ls[i1] * 100 << "; " << qTofd1ls[i1] << ", " << endl;
                            cout << "Fib23b: " << yFi23b[j] * 100. << endl;
                            cout << "ycheck1:" << ycheck << ", " << abs(yTofd1ls[i1] * 100. - ycheck) << endl;
                        }

                        if (fPairs && std::abs(yTofd1ls[i1] * 100. - ycheck) > dytofdfib23b)
                            continue;

                        if ((!fPairs) ||
                            (fPairs && tempFi23bl)) //(fPairs && abs(yFi23b[j] * 100.) > 0.26 && tempFi23bl))
                        {

                            // if (fGraphCuts && !cut_fi30_fi23b->IsInside(yFi23b[j] * 100., xFi30[i] * 100.))
                            //	continue;
                            fh_tofdleft_vs_Fib23b_ac->Fill(yFi23b[j] * 100., yTofd1ls[i1] * 100.);
                            fh_Fib30_vs_Fib23b_dx_ac->Fill(yFi23b[j] * 100., yFi30[i] * 100. - yFi23b[j] * 100.);

                            fh_Fib_Time_ac[fi23b]->Fill(yFi23b[j] * 100., timeFi23b[j]);
                            fh_xy_Fib_ac[fi23b]->Fill(xFi23b[j] * 100., yFi23b[j] * 100.);
                            fh_ToT_Fib_ac[fi23b]->Fill(yFi23b[j] * 100., qFi23b[j]);
                            fh_Fibs_vs_Tofd_ac[fi23b]->Fill(yTofd1ls[i1] * 100., yFi23b[j] * 100.);
                            fh_Fib_ToF_ac[fi23b]->Fill(yFi23b[j] * 100., tTofd1ls[i1] - timeFi23b[j]);
                            fh_Fib23b_vs_Fib3X_tot->Fill(qFi23b[j], qFi30[i]);
                            fh_Fib_ToF_vs_Qtofd[fi23b]->Fill(qTofd1ls[i1], tTofd1ls[i1] - timeFi23b[j]);

                            if (debug_fibcut)
                                cout << "Fib23b after fi30 cuts: " << j << " y: " << yFi23b[j] << " q: " << qFi23b[j]
                                     << endl;

                            cond6 = kTRUE;
                            if (debug_fibcut)
                                cout << "cond6" << endl;

                            if (!maxWerte)
                            {
                                if (!fFi23b[j])
                                {
                                    detector[countdet] = fi23b;
                                    xdet[countdet] = xFi23b[j];
                                    ydet[countdet] = yFi23b[j];
                                    zdet[countdet] = 0.;
                                    qdet[countdet] = qFi23b[j];
                                    tdet[countdet] = tTofd1ls[i1] - timeFi23b[j];
                                    countdet++;
                                    fFi23b[j] = true;
                                }
                            }
                        }
                    }
                }
            }

            if (debug_fibcut)
                cout << "**Loop over Tofd2l ****" << endl;

            for (Int_t i1 = 0; i1 < i2lc; i1++)
            {
                if (qTofd2ls[i1] > 0 && !cTofd2ls[i1])
                {
                    detector[countdet] = tofd2lc;
                    xdet[countdet] = xTofd2ls[i1];
                    ydet[countdet] = yTofd2ls[i1];
                    zdet[countdet] = 0.;
                    tdet[countdet] = tTofd2ls[i1];
                    qdet[countdet] = (int)(qTofd2ls[i1] + 0.5); // q for tracker must be integer

                    if (debug_fibcut)
                        cout << "Tofd2l saved in fibcuts: " << xdet[countdet] * 100 << "; " << ydet[countdet] * 100
                             << "; " << qdet[countdet] << endl;

                    //  fh_xy_tofd_ac->Fill(x1[tofd1l] * 100., y1[tofd1l] * 100.);

                    countdet++;
                    cTofd2ls[i1] = true;
                }

                for (Int_t i = 0; i < mult30; i++)
                {

                    fh_tofd_vs_Fib30->Fill(xFi30[i] * 100., xTofd2ls[i1] * 100.);
                    fh_ytofd_vs_yFib30->Fill(yFi30[i] * 100., yTofd2ls[i1] * 100.);
                    fh_qtofd_vs_qFib[fi30]->Fill(qFi30[i], qTofd2ls[i1]);

                    Double_t xtemp30;
                    xtemp30 = 0.7796 * xTofd2ls[i1] * 100. - 28.669; // erste ToFD Ebene

                    if (debug_fibcut)
                    {
                        cout << "Fib30 before cut: " << i << " x: " << xFi30[i] * 100. << " q: " << qFi30[i] << endl;
                        cout << "tofd2l before fi30 cut: " << xTofd2ls[i1] * 100. << ", " << yTofd2ls[i1] * 100. << ", "
                             << qTofd2ls[i1] << endl;
                        cout << "Diff: " << std::abs(xtemp30 - xFi30[i] * 100.) << ", " << dx3 << endl;
                    }

                    if (fPairs && std::abs(xtemp30 - xFi30[i] * 100.) > dx3)
                        continue;
                    if ((tTofd2ls[i1] - timeFi30[i]) > 20. && (tTofd2ls[i1] - timeFi30[i]) < -25.)
                        continue;

                    if (debug_fibcut)
                    {
                        cout << "Fib30after fibcuts: " << i << " x: " << xFi30[i] << " q: " << qFi30[i] << endl;
                        cout << "Tofd2l after fi30 cuts: " << xTofd2ls[i1] * 100. << ", " << yTofd2ls[i1] * 100. << ", "
                             << qTofd2ls[i1] << endl;
                    }

                    fh_tofd_vs_Fib30_ac->Fill(xFi30[i] * 100., xTofd2ls[i1] * 100.);
                    fh_Fib_ToF_ac[fi30]->Fill(xFi30[i] * 100., tTofd2ls[i1] - timeFi30[i]);
                    fh_xy_Fib_ac[fi30]->Fill(xFi30[i] * 100., yFi30[i] * 100.);
                    fh_ToT_Fib_ac[fi30]->Fill(xFi30[i] * 100., qFi30[i]);
                    fh_Fibs_vs_Tofd_ac[fi30]->Fill(xTofd2ls[i1] * 100. + randx, xFi30[i] * 100.);
                    fh_Fib_Time_ac[fi30]->Fill(xFi30[i] * 100., timeFi30[i]);
                    fh_qtofd_vs_qFib_ac[fi30]->Fill(qFi30[i], qTofd2ls[i1]);
                    fh_Fib_ToF_vs_Qtofd[fi30]->Fill(qTofd2ls[i1], tTofd2ls[i1] - timeFi30[i]);

                    for (Int_t j = 0; j < mult32; j++)
                    {
                        fh_Fib32_vs_Fib30->Fill(xFi30[i] * 100., xFi32[j] * 100.);
                        fh_Fib32_vs_Fib30_dx->Fill(xFi30[i] * 100., xFi32[j] * 100. - xFi30[i] * 100.);
                        fh_Fib32_vs_Fib30_dt->Fill(xFi32[j] * 100., timeFi30[i] - timeFi32[j]);
                        fh_qtofd_vs_qFib[fi32]->Fill(qFi32[j], qTofd2ls[i1]);

                        Double_t x30 =
                            0.873 * xFi32[j] * 100. - 3.302 - offset_corr; // 1.10926 * xFi30[i] * 100. + 2.8943;
                        tempFi32 = false;
                        //	if (fSimu)
                        //	tempFi32 = true;
                        if (abs(timeFi30[i] - timeFi32[j]) < dtft && (tTofd2ls[i1] - timeFi32[j]) < 20. &&
                            (tTofd2ls[i1] - timeFi32[j]) > -25.)
                            tempFi32 = true;

                        if ((abs(xFi30[i] * 100. - x30) < dx1 && xFi32[j] * 100. > -30 && xFi30[i] * 100. > -30. &&
                             tempFi32))
                        {
                            fh_Fib32_vs_Fib30_ac->Fill(xFi30[i] * 100., xFi32[j] * 100.);
                            fh_Fib32_vs_Fib30_dx_ac->Fill(xFi30[i] * 100., xFi32[j] * 100. - xFi30[i] * 100.);

                            fh_Fib_Time_ac[fi32]->Fill(xFi32[j] * 100., timeFi32[j]);
                            fh_xy_Fib_ac[fi32]->Fill(xFi32[j] * 100., yFi32[j] * 100.);
                            fh_ToT_Fib_ac[fi32]->Fill(xFi32[j] * 100., qFi32[j]);
                            fh_Fibs_vs_Tofd_ac[fi32]->Fill(xTofd2ls[i1] * 100. + randx, xFi32[j] * 100.);
                            fh_Fib_ToF_ac[fi32]->Fill(xFi32[j] * 100., tTofd2ls[i1] - timeFi32[j]);
                            fh_qtofd_vs_qFib_ac[fi32]->Fill(qFi32[j], qTofd2ls[i1]);
                            fh_Fib32_vs_Fib30_tot->Fill(qFi32[j], qFi30[i]);
                            fh_Fib_ToF_vs_Qtofd[fi32]->Fill(qTofd2ls[i1], tTofd2ls[i1] - timeFi32[j]);

                            if (debug_fibcut)
                                cout << "Fib32 after fi30 cuts: " << j << " x: " << xFi32[j] << " q: " << qFi32[j]
                                     << endl;

                            cond4 = kTRUE;
                            if (debug_fibcut)
                                cout << "cond4" << endl;

                            if (!maxWerte)
                            {
                                if (!fFi30[i])
                                {
                                    Double_t qtemp = qFi30[i];
                                    if (fGraphCuts && !fSimu)
                                    {
                                        if (cut_qfi30_qfi32->IsInside(qFi32[j], qFi30[i]))
                                            qtemp = 6.;
                                        else
                                            qtemp = 2;
                                    }
                                    detector[countdet] = fi30;
                                    xdet[countdet] = xFi30[i];
                                    ydet[countdet] = yFi30[i];
                                    zdet[countdet] = 0.;
                                    qdet[countdet] = qtemp; // qFi30[i];
                                    tdet[countdet] = tTofd2ls[i1] - timeFi30[i];
                                    countdet++;
                                    fFi30[i] = true;
                                }
                                if (!fFi32[j])
                                {
                                    Double_t qtemp = qFi32[j];
                                    if (fGraphCuts && !fSimu)
                                    {
                                        if (cut_qfi30_qfi32->IsInside(qFi32[j], qFi30[i]))
                                            qtemp = 6.;
                                        else
                                            qtemp = 2;
                                    }
                                    detector[countdet] = fi32;
                                    xdet[countdet] = xFi32[j];
                                    ydet[countdet] = yFi32[j];
                                    zdet[countdet] = 0.;
                                    qdet[countdet] = qtemp; // qFi32[j];
                                    tdet[countdet] = tTofd2ls[i1] - timeFi32[j];
                                    countdet++;
                                    fFi32[j] = true;
                                }
                            }
                        }
                    }
                    for (Int_t j = 0; j < mult23a; j++)
                    {
                        fh_Fib30_vs_Fib23a->Fill(xFi23a[j] * 100., xFi30[i] * 100.);
                        fh_Fib30_vs_Fib23a_dx->Fill(xFi23a[j] * 100., xFi30[i] * 100. - xFi23a[j] * 100.);
                        fh_Fib30_vs_Fib23a_dt->Fill(xFi30[i] * 100. - xFi23a[j] * 100., timeFi30[i] - timeFi23a[j]);

                        Double_t x30 = 6.498 * xFi23a[j] * 100. - 28.93;
                        // if (abs(xFi30[i] * 100. - x30) < dx2 && xFi30[i] > -100 && xFi23a[j] > -100)
                        Double_t x23_ridge = 4.27564e+00 + 1.48326e-01 * xFi30[i] * 100;

                        tempFi23al = false;
                        // if (fSimu)
                        //  tempFi23al = true;
                        if ((tTofd2ls[i1] - timeFi23a[j]) < 28. && (tTofd2ls[i1] - timeFi23a[j]) > -25.)
                            tempFi23al = true;

                        if ((!fPairs) ||
                            //  (fPairs && (xFi23a[j] * 100.) > 0.26 && tempFi23al)) // && abs(xFi30[i] * 100. - x30) <
                            //  dx2))
                            (fPairs && tempFi23al && abs(xFi23a[j] * 100. - x23_ridge) < dx4))
                        {
                            // if (fGraphCuts && !cut_fi30_fi23b->IsInside(yFi23b[j] * 100., xFi30[i] * 100.))
                            //	continue;

                            fh_Fib30_vs_Fib23a_ac->Fill(xFi23a[j] * 100., xFi30[i] * 100.);
                            fh_Fib30_vs_Fib23a_dx_ac->Fill(xFi23a[j] * 100., xFi30[i] * 100. - xFi23a[j] * 100.);

                            fh_Fib_Time_ac[fi23a]->Fill(xFi23a[j] * 100., timeFi23a[j]);
                            fh_xy_Fib_ac[fi23a]->Fill(xFi23a[j] * 100., yFi23a[j] * 100.);
                            fh_ToT_Fib_ac[fi23a]->Fill(xFi23a[j] * 100., qFi23a[j]);
                            fh_Fibs_vs_Tofd_ac[fi23a]->Fill(xTofd2ls[i1] * 100., xFi23a[j] * 100.);
                            fh_Fib_ToF_ac[fi23a]->Fill(xFi23a[j] * 100., tTofd2ls[i1] - timeFi23a[j]);
                            fh_Fib23a_vs_Fib3X_tot->Fill(qFi23a[j], qFi30[i]);
                            fh_Fib_ToF_vs_Qtofd[fi23a]->Fill(qTofd2ls[i1], tTofd2ls[i1] - timeFi23a[j]);

                            if (debug_fibcut)
                                cout << "Fib23a after fi30 cuts: " << j << " x: " << xFi23a[j] << " q: " << qFi23a[j]
                                     << endl;

                            cond5 = kTRUE;
                            if (debug_fibcut)
                                cout << "cond5" << endl;

                            if (!maxWerte)
                            {
                                if (!fFi23a[j])
                                {
                                    detector[countdet] = fi23a;
                                    xdet[countdet] = xFi23a[j];
                                    ydet[countdet] = yFi23a[j];
                                    zdet[countdet] = 0.;
                                    qdet[countdet] = qFi23a[j];
                                    tdet[countdet] = tTofd2ls[i1] - timeFi23a[j];
                                    countdet++;

                                    detector[countdet] = fi23al;
                                    xdet[countdet] = xFi23a[j];
                                    ydet[countdet] = yFi23a[j];
                                    zdet[countdet] = 0.;
                                    qdet[countdet] = qFi23a[j];
                                    tdet[countdet] = tTofd2ls[i1] - timeFi23a[j];
                                    countdet++;
                                    fFi23a[j] = true;
                                }
                            }
                        }
                    }
                    for (Int_t j = 0; j < mult23b; j++)
                    {
                        // Double_t x30 = 6.98386 * xFi23b[j] * 100. - 27.39897;
                        // if (abs(xFi30[i] * 100. - x30) < dx2 && xFi30[i] > -100 && xFi23b[j] > -100)

                        fh_tofdleft_vs_Fib23b->Fill(yFi23b[j] * 100., yTofd2ls[i1] * 100.);
                        fh_Fib30_vs_Fib23b_dx->Fill(yFi23b[j] * 100., yFi30[i] * 100. - yFi23b[j] * 100.);
                        fh_Fib30_vs_Fib23b_dt->Fill(yFi23b[j] * 100., timeFi30[i] - timeFi23b[j]);

                        tempFi23bl = false;
                        // if (fSimu)
                        //  tempFi23bl = true;
                        if ((tTofd2ls[i1] - timeFi23b[j]) < 25. && (tTofd2ls[i1] - timeFi23b[j]) > -25.)
                            tempFi23bl = true;

                        Double_t ycheck = 7.6358668 * yFi23b[j] * 100. - 0.0817257;

                        if (debug_fibcut)
                        {
                            cout << "Before Tofd vs fi23b: " << endl;
                            cout << "Tofd1l: " << yTofd2ls[i1] * 100 << "; " << qTofd2ls[i1] << ", " << endl;
                            cout << "Fib23b: " << yFi23b[j] * 100. << endl;
                            cout << "ycheck1:" << ycheck << ", " << abs(yTofd2ls[i1] * 100. - ycheck) << endl;
                        }

                        if (fPairs && std::abs(yTofd2ls[i1] * 100. - ycheck) > dytofdfib23b)
                            continue;

                        if ((!fPairs) ||
                            (fPairs && tempFi23bl)) //(fPairs && abs(yFi23b[j] * 100.) > 0.26 && tempFi23bl))
                        {
                            // if (fGraphCuts && !cut_fi30_fi23b->IsInside(yFi23b[j] * 100., xFi30[i] * 100.))
                            //	continue;

                            fh_tofdleft_vs_Fib23b_ac->Fill(yFi23b[j] * 100., yTofd2ls[i1] * 100.);
                            fh_Fib30_vs_Fib23b_dx_ac->Fill(yFi23b[j] * 100., yFi30[i] * 100. - yFi23b[j] * 100.);

                            fh_Fib_Time_ac[fi23b]->Fill(yFi23b[j] * 100., timeFi23b[j]);
                            fh_xy_Fib_ac[fi23b]->Fill(xFi23b[j] * 100., yFi23b[j] * 100.);
                            fh_ToT_Fib_ac[fi23b]->Fill(yFi23b[j] * 100., qFi23b[j]);
                            fh_Fibs_vs_Tofd_ac[fi23b]->Fill(yTofd2ls[i1] * 100., yFi23b[j] * 100.);
                            fh_Fib_ToF_ac[fi23b]->Fill(yFi23b[j] * 100., tTofd2ls[i1] - timeFi23b[j]);
                            fh_Fib23b_vs_Fib3X_tot->Fill(qFi23b[j], qFi30[i]);
                            fh_Fib_ToF_vs_Qtofd[fi23b]->Fill(qTofd2ls[i1], tTofd2ls[i1] - timeFi23b[j]);

                            if (debug_fibcut)
                                cout << "Fib23b after fi30 cuts: " << j << " y: " << yFi23b[j] << " q: " << qFi23b[j]
                                     << endl;

                            cond6 = kTRUE;
                            if (debug_fibcut)
                                cout << "cond6" << endl;

                            if (!maxWerte)
                            {
                                if (!fFi23b[j])
                                {
                                    detector[countdet] = fi23b;
                                    xdet[countdet] = xFi23b[j];
                                    ydet[countdet] = yFi23b[j];
                                    zdet[countdet] = 0.;
                                    qdet[countdet] = qFi23b[j];
                                    tdet[countdet] = tTofd2ls[i1] - timeFi23b[j];
                                    countdet++;
                                    fFi23b[j] = true;
                                }
                            }
                        }
                    }
                }
            }

            if (debug_fibcut)
                cout << "**Loop over Tofd3l (Tofi) ****" << endl;

            for (Int_t i1 = 0; i1 < i3lc; i1++)
            {
                if (qTofd3ls[i1] > 0 && !cTofd3ls[i1])
                {
                    detector[countdet] = tofd3lc;
                    xdet[countdet] = xTofd3ls[i1];
                    ydet[countdet] = yTofd3ls[i1];
                    zdet[countdet] = 0.;
                    tdet[countdet] = tTofd3ls[i1];
                    qdet[countdet] = (int)(qTofd3ls[i1] + 0.5); // q for tracker must be integer

                    if (debug_fibcut)
                        cout << "Tofd3l saved in fibcuts: " << xdet[countdet] * 100 << "; " << ydet[countdet] * 100
                             << "; " << qdet[countdet] << endl;

                    //  fh_xy_tofd_ac->Fill(x1[tofd3l] * 100., y1[tofd3l] * 100.);

                    countdet++;
                    cTofd3ls[i1] = true;
                }
            }
        }

        if (countdet > 50)
        {
            countdet50 += 1;
            if (debug_write)
                cout << "Too many hits!!!" << endl;
            continue;
        }

        Int_t ijk = 0;
        for (Int_t i = 0; i < countdet; i++)
        {
            if (detector[i] == fi23ar || detector[i] == fi23al)
                ijk++;
        }
        Nhits_after_fibcor = Nhits_after_fibcor + countdet - ijk;

        Bool_t temp_cond = false;

        if (fPairs && !fibCuts && (mult30 + mult32) != 1 && (mult31 + mult33) != 1 && mult23a > 1 && mult23b > 1)
            temp_cond = true;

        //         if (fPairs &&  fibCuts && (mult30 + mult32) != 1 && (mult31 + mult33) != 1 && mult23a > 1 && mult23b
        //         > 1 &&
        //             ((cond1 && cond2 && cond3) || (cond4 && cond5 && cond6))) temp_cond = true;

        if (fPairs && fibCuts && (mult30 + mult32) != 1 && (mult31 + mult33) != 1)
            temp_cond = true;

        if (!fPairs && !fibCuts && maxWerte && ((maxWertFi30 && maxWertFi32) || (maxWertFi31 && maxWertFi33)))
            temp_cond = true;

        if (!fPairs && fibCuts && maxWerte && ((maxWertFi30 && maxWertFi32) || (maxWertFi31 && maxWertFi33)) &&
            ((cond1 && cond2 && cond3) || (cond4 && cond5 && cond6)))
            temp_cond = true;

        if (temp_cond)
            Neventafterfibcor++;

        if (fB != -1710 && !fPairs && temp_cond)
        {
            Int_t ncount[10] = { 0 };

            Bool_t fcut = false;
            Bool_t fcut1 = false;
            Bool_t fcut2 = false;
            for (Int_t i = 0; i < countdet; i++)
            {
                ncount[detector[i]] += 1;

                if ((abs(fB) < 1710. && (detector[i] == 2 && (xdet[i] * 100. > fX1min && xdet[i] * 100. < fX1max &&
                                                              ydet[i] * 100. < 4. && ydet[i] * 100. > -4.))) ||
                    (abs(fB) > 1710. && (detector[i] == 3 && (xdet[i] * 100. > fX1min && xdet[i] * 100. < fX1max &&
                                                              ydet[i] * 100. < 4. && ydet[i] * 100. > -4.))))
                {
                    fcut1 = true;
                    if (debug_write)
                        cout << "Fib30/31 && fcut: " << detector[i] << ", " << xdet[i] * 100. << ", " << ydet[i] * 100.
                             << "; " << fcut << ", " << abs(fB) << ", " << fX1min << ";" << fX1max << endl;
                }
                if ((abs(fB) < 1710. && (detector[i] == 4 && (ydet[i] * 100. < 4. && ydet[i] * 100. > -4. &&
                                                              xdet[i] * 100. > fX2min && xdet[i] * 100. < fX2max))) ||
                    (abs(fB) > 1710. && (detector[i] == 5 && (ydet[i] * 100. < 4. && ydet[i] * 100. > -4. &&
                                                              xdet[i] * 100. > fX2min && xdet[i] * 100. < fX2max))))
                {
                    fcut2 = true;
                    if (debug_write)
                        cout << "Fib32/33 && fcut: " << detector[i] << ", " << xdet[i] * 100. << ", " << ydet[i] * 100.
                             << "; " << fcut << ", " << abs(fB) << ", " << fX2min << ";" << fX2max << endl;
                }
            }
            if (fcut1 && fcut2)
                fcut = true;

            if (fSimu)
                fcut = true;

            countdet_s = 0;

            Double_t qdet_s7 = 0., qdet_s6 = 0., qdet_s8 = 0., qdet_s9 = 0.;

            if (fSimu ||
                (!fSimu && ((ncount[9] > 0 && ncount[7] > 0 && ncount[2] > 0 && ncount[4] > 0 && ncount[0] > 0 &&
                             ncount[1] > 0 && ncount[3] < 1 && ncount[5] < 1 && ncount[9] == ncount[7] &&
                             ncount[7] == ncount[2] && ncount[2] == ncount[4] && ncount[6] < 1 && ncount[8] < 1 &&
                             ncount[3] < 1 && ncount[5] < 1 && abs(fB) < 1710. && fcut) ||
                            (ncount[6] > 0 && ncount[8] > 0 && ncount[3] > 0 && ncount[5] > 0 && ncount[0] > 0 &&
                             ncount[1] > 0 && ncount[2] < 1 && ncount[4] < 1 && ncount[8] == ncount[6] &&
                             ncount[6] == ncount[3] && ncount[3] == ncount[5] && ncount[7] < 1 && ncount[9] < 1 &&
                             ncount[2] < 1 && ncount[4] < 1 && abs(fB) > 1710. && fcut))))
            // if(1 == 1)
            {
                Bool_t goodQ = true;

                for (Int_t i = 0; i < countdet; i++)
                {

                    if (debug_write)
                        cout << "After if - Count: " << i << ", det: " << detector[i] << ", qdet: " << qdet[i]
                             << ", xdet: " << xdet[i] << endl;

                    Bool_t iffib33 = true;
                    if (detector[i] > 5 && detector[i] < 10 && goodQ)
                    {
                        xdet_s[countdet_s] = xdet[i];
                        ydet_s[countdet_s] = ydet[i];
                        zdet_s[countdet_s] = zdet[i];
                        qdet_s[countdet_s] = qdet[i];
                        tdet_s[countdet_s] = tdet[i];
                        detector_s[countdet_s] = detector[i];

                        nsum[detector_s[countdet_s]] += 1;
                        xdet_sum[detector_s[countdet_s]] += xdet[i];
                        ydet_sum[detector_s[countdet_s]] += ydet[i];
                        zdet_sum[detector_s[countdet_s]] += zdet[i];
                        qdet_sum[detector_s[countdet_s]] += qdet[i];
                        tdet_sum[detector_s[countdet_s]] += tdet[i];

                        if (debug_write)
                        {
                            cout << "tofd det:   " << detector[i] << ", count " << i << ", x: " << xdet[i]
                                 << ", q: " << qdet[i] << endl;
                            cout << "tofd det_s: " << detector_s[countdet_s] << ", count " << countdet_s
                                 << ", x: " << xdet_s[countdet_s] << ", q: " << qdet_s[countdet_s] << endl;
                        }
                        countdet_s++;
                    }

                    else if (detector[i] > 1 && detector[i] < 6 && goodQ)
                    {
                        xdet_s[countdet_s] = xdet[i];
                        ydet_s[countdet_s] = ydet[i];
                        zdet_s[countdet_s] = zdet[i];
                        qdet_s[countdet_s] = 8.; // qdet[i];
                        tdet_s[countdet_s] = tdet[i];
                        detector_s[countdet_s] = detector[i];

                        nsum[detector_s[countdet_s]] += 1;
                        xdet_sum[detector_s[countdet_s]] += xdet[i];
                        ydet_sum[detector_s[countdet_s]] += ydet[i];
                        zdet_sum[detector_s[countdet_s]] += zdet[i];
                        qdet_sum[detector_s[countdet_s]] += 8.; // qdet[i];
                        tdet_sum[detector_s[countdet_s]] += tdet[i];

                        if (debug_write)
                        {
                            cout << "fib  det:   " << detector[i] << ", count " << i << ", x: " << xdet[i]
                                 << ", q: " << qdet[i] << endl;
                            cout << "fib  det_s: " << detector_s[countdet_s] << ", count " << countdet_s
                                 << ", x: " << xdet_s[countdet_s] << ", q: " << qdet_s[countdet_s] << endl;
                        }
                        countdet_s++;
                    }
                    else if (detector[i] < 2 && goodQ)
                    {
                        xdet_s[countdet_s] = xdet[i];
                        ydet_s[countdet_s] = ydet[i];
                        zdet_s[countdet_s] = zdet[i];
                        qdet_s[countdet_s] = 8; // qdet[i];
                        tdet_s[countdet_s] = tdet[i];
                        detector_s[countdet_s] = detector[i];

                        nsum[detector_s[countdet_s]] += 1;
                        xdet_sum[detector_s[countdet_s]] += xdet[i];
                        ydet_sum[detector_s[countdet_s]] += ydet[i];
                        zdet_sum[detector_s[countdet_s]] += zdet[i];
                        qdet_sum[detector_s[countdet_s]] += qdet[i];
                        tdet_sum[detector_s[countdet_s]] += tdet[i];

                        if (debug_write)
                        {
                            cout << "fib  det:   " << detector[i] << ", count " << i << ", x: " << xdet[i]
                                 << ", q: " << qdet[i] << endl;
                            cout << "fib  det_s: " << detector_s[countdet_s] << ", count " << countdet_s
                                 << ", x: " << xdet_s[countdet_s] << ", q: " << qdet_s[countdet_s] << endl;
                        }

                        countdet_s++;
                    }

                    fNeventselect += 1;
                }

                // Here write hit data of all detectors
                Double_t ax = 0.;
                Double_t ay = 0.;
                Double_t aq = 0.;
                Double_t at = 0.;
                Double_t cx = 0.;
                Double_t cy = 0.;
                Double_t cq = 0.;
                Double_t ct = 0.;

                counter1++;

                //  cout<<"************* counter1: "<<counter1<<endl;
                // if (counter1 < 801 && !fAverage)
                if (!fAverage)
                {
                    for (Int_t i = 0; i < countdet_s; i++)
                    {
                        if (debug_write)
                            cout << "event: " << fNEvents << ", counter1 :" << counter1 << " #" << i
                                 << " Det: " << detector_s[i] << " x: " << xdet_s[i] * 100.
                                 << " y: " << ydet_s[i] * 100. << " q: " << qdet_s[i] << endl;
                        if (detector_s[i] == 2)
                            counter2++;
                        if (detector_s[i] == 7)
                            counter3++;

                        if (detector_s[i] == 0)
                        {
                            countdet_written[0]++;
                            // write fiber detector hits
                            new ((*fFi23aHitItems)[fNofFi23aHitItems++]) R3BFiberMAPMTHitData(
                                0, xdet_s[i] * 100., ydet_s[i] * 100., qdet_s[i], tdet_s[i], 0, 0, 0., 0, 0.);
                            fh_check_QvsX[detector_s[i]]->Fill(xdet_s[i] * 100, qdet_s[i]);
                            fh_check_TvsX[detector_s[i]]->Fill(xdet_s[i] * 100, tdet_s[i]);
                            fh_check_XvsY[detector_s[i]]->Fill(xdet_s[i] * 100, ydet_s[i] * 100);
                        }
                        if (detector_s[i] == 1)
                        {
                            countdet_written[1]++;
                            // write fiber detector hits
                            new ((*fFi23bHitItems)[fNofFi23bHitItems++]) R3BFiberMAPMTHitData(
                                0, xdet_s[i] * 100., ydet_s[i] * 100., qdet_s[i], tdet_s[i], 0, 0, 0., 0, 0.);
                            fh_check_QvsX[detector_s[i]]->Fill(ydet_s[i] * 100, qdet_s[i]);
                            fh_check_TvsX[detector_s[i]]->Fill(ydet_s[i] * 100, tdet_s[i]);
                            fh_check_XvsY[detector_s[i]]->Fill(xdet_s[i] * 100, ydet_s[i] * 100);
                        }
                        if (detector_s[i] == 2 && abs(fB) < 1710.)
                        {
                            countdet_written[2]++;
                            // write fiber detector hits
                            new ((*fFi30HitItems)[fNofFi30HitItems++]) R3BFiberMAPMTHitData(
                                0, xdet_s[i] * 100., ydet_s[i] * 100., qdet_s[i], tdet_s[i], 0, 0, 0., 0, 0.);
                            fh_check_QvsX[detector_s[i]]->Fill(xdet_s[i] * 100, qdet_s[i]);
                            fh_check_TvsX[detector_s[i]]->Fill(xdet_s[i] * 100, tdet_s[i]);
                            fh_check_XvsY[detector_s[i]]->Fill(xdet_s[i] * 100, ydet_s[i] * 100);
                        }
                        if (detector_s[i] == 3 && abs(fB) > 1710.)
                        {
                            countdet_written[3]++;
                            // write fiber detector hits
                            new ((*fFi31HitItems)[fNofFi31HitItems++]) R3BFiberMAPMTHitData(
                                0, xdet_s[i] * 100., ydet_s[i] * 100., qdet_s[i], tdet_s[i], 0, 0, 0., 0, 0.);
                            fh_check_QvsX[detector_s[i]]->Fill(xdet_s[i] * 100, qdet_s[i]);
                            fh_check_TvsX[detector_s[i]]->Fill(xdet_s[i] * 100, tdet_s[i]);
                            fh_check_XvsY[detector_s[i]]->Fill(xdet_s[i] * 100, ydet_s[i] * 100);
                        }
                        if (detector_s[i] == 4 && abs(fB) < 1710.)
                        {
                            countdet_written[4]++;
                            // write fiber detector hits
                            new ((*fFi32HitItems)[fNofFi32HitItems++]) R3BFiberMAPMTHitData(
                                0, xdet_s[i] * 100., ydet_s[i] * 100., qdet_s[i], tdet_s[i], 0, 0, 0., 0, 0.);
                            fh_check_QvsX[detector_s[i]]->Fill(xdet_s[i] * 100, qdet_s[i]);
                            fh_check_TvsX[detector_s[i]]->Fill(xdet_s[i] * 100, tdet_s[i]);
                            fh_check_XvsY[detector_s[i]]->Fill(xdet_s[i] * 100, ydet_s[i] * 100);
                        }
                        if (detector_s[i] == 5 && abs(fB) > 1710.)
                        {
                            countdet_written[5]++;
                            // write fiber detector hits
                            //  if(isumdet5 == 0){
                            new ((*fFi33HitItems)[fNofFi33HitItems++]) R3BFiberMAPMTHitData(
                                0, xdet_s[i] * 100., ydet_s[i] * 100., qdet_s[i], tdet_s[i], 0, 0, 0., 0, 0.);
                            fh_check_QvsX[detector_s[i]]->Fill(xdet_s[i] * 100, qdet_s[i]);
                            fh_check_TvsX[detector_s[i]]->Fill(xdet_s[i] * 100, tdet_s[i]);
                            fh_check_XvsY[detector_s[i]]->Fill(xdet_s[i] * 100, ydet_s[i] * 100);
                            //	}
                            //      isumdet5 += 1;
                        }
                        if (detector_s[i] == 6 || detector_s[i] == 8 || detector_s[i] == 7 || detector_s[i] == 9)
                        {

                            if (detector_s[i] == 6)
                                countdet_written[6]++;
                            if (detector_s[i] == 7)
                                countdet_written[7]++;
                            if (detector_s[i] == 8)
                                countdet_written[8]++;
                            if (detector_s[i] == 9)
                                countdet_written[9]++;

                            if (ftrackerType == 0) // TofD data written out for each plane separately
                            {
                                Int_t ipl;
                                if (detector_s[i] == 6 || detector_s[i] == 7)
                                    ipl = 1;
                                else
                                    ipl = 2;

                                if (((detector_s[i] == 6 || detector_s[i] == 8) && abs(fB) > 1710.) ||
                                    ((detector_s[i] == 7 || detector_s[i] == 9) && abs(fB) < 1710.))
                                {
                                    new ((*fTofdHitItems)[fNofTofdHitItems++]) R3BTofdHitData(tdet_s[i],
                                                                                              xdet_s[i] * 100.,
                                                                                              ydet_s[i] * 100.,
                                                                                              qdet_s[i],
                                                                                              -5,
                                                                                              qdet_s[i],
                                                                                              ipl,
                                                                                              1,
                                                                                              0);

                                    fh_check_QvsX[6]->Fill(xdet_s[i] * 100, qdet_s[i]);
                                    fh_check_TvsX[6]->Fill(xdet_s[i] * 100, tdet_s[i]);
                                    fh_check_XvsY[6]->Fill(xdet_s[i] * 100, ydet_s[i] * 100);
                                }
                            }
                            else // first and second plane will be written as one detector; this makes Dima's tracker
                                 // faster
                            {
                                if (((detector_s[i] == 6 || detector_s[i] == 8) && abs(fB) > 1710.) ||
                                    ((detector_s[i] == 7 || detector_s[i] == 9) && abs(fB) < 1710.))
                                {
                                    ax += xdet_s[i];
                                    ay += ydet_s[i];
                                    at += tdet_s[i];
                                    aq += qdet_s[i];
                                }
                            }
                        }
                    }
                    if (counter2 != counter1)
                    {
                        // cout << "Fib 30: event: "<<fNEvents<<", counters1-3 :" << counter1 <<", "<<counter2<<",
                        // "<<counter3<<endl;
                        counter2 = counter1;
                    }
                    if (counter3 != counter1)
                    {
                        // cout << "TofD1: event: "<<fNEvents<<", counters1-3 :" << counter1 <<", "<<counter2<<",
                        // "<<counter3<<endl;
                        counter3 = counter1;
                    }
                    counter1mem = counter1;

                    if (ftrackerType == 1)
                    {
                        ax = ax / 2.;
                        ay = ay / 2.;
                        at = at / 2.;
                        aq = aq / 2.;

                        new ((*fTofdHitItems)[fNofTofdHitItems++])
                            R3BTofdHitData(at, ax * 100., ay * 100., aq, -5., aq, 10, 1, 0);

                        // Fill check spectra:
                        fh_check_QvsX[6]->Fill(ax * 100, aq);
                        fh_check_TvsX[6]->Fill(ax * 100, at);
                        fh_check_XvsY[6]->Fill(ax * 100, ay * 100);
                    }
                }
            }
        }

        if (fPairs && temp_cond)
        {
            Int_t ncount[n_det] = { 0 }, ncounttemp6 = 0, ncounttemp7 = 0, ncounttemp8 = 0, ncounttemp9 = 0,
                  ncounttemp10 = 0, ncounttemp11 = 0;
            Double_t x23ltemp[40], x23rtemp[40];

            Double_t qsum6 = 0., qsum7 = 0, qsum8 = 0, qsum9 = 0, qsum10 = 0, qsum11 = 0;
            for (Int_t i = 0; i < 40; i++)
            {
                xtofdtemp6[i] = -1000;
                xtofdtemp7[i] = -1000;
                xtofdtemp8[i] = -1000;
                xtofdtemp9[i] = -1000;
                xtofdtemp10[i] = -1000;
                xtofdtemp11[i] = -1000;
                ytofdtemp6[i] = -1000;
                ytofdtemp7[i] = -1000;
                ytofdtemp8[i] = -1000;
                ytofdtemp9[i] = -1000;
                ytofdtemp10[i] = -1000;
                ytofdtemp11[i] = -1000;
                qtofdtemp6[i] = -1000;
                qtofdtemp7[i] = -1000;
                qtofdtemp8[i] = -1000;
                qtofdtemp9[i] = -1000;
                qtofdtemp10[i] = -1000;
                qtofdtemp11[i] = -1000;
                ttofdtemp6[i] = -1000;
                ttofdtemp7[i] = -1000;
                ttofdtemp8[i] = -1000;
                ttofdtemp9[i] = -1000;
                ttofdtemp10[i] = -1000;
                ttofdtemp11[i] = -1000;
            }

            fNeventpair++;

            Bool_t bar24plane2 = false;
            Bool_t bar21plane2 = false;
            Bool_t edgetofd = false;
            for (Int_t i = 0; i < countdet; i++)
            {
                if (detector[i] > 5 && detector[i] < 10)
                    fh_xy_tofd_ac->Fill(xdet[i] * 100., ydet[i] * 100.);

                if (detector[i] == 14)
                    x23rtemp[ncount[detector[i]]] = xdet[i];
                if (detector[i] == 15)
                    x23ltemp[ncount[detector[i]]] = xdet[i];

                if (detector[i] == 6)
                {
                    xtofdtemp6[ncounttemp6] = xdet[i];
                    ytofdtemp6[ncounttemp6] = ydet[i];
                    qtofdtemp6[ncounttemp6] = qdet[i];
                    ttofdtemp6[ncounttemp6] = tdet[i];
                    if (xtofdtemp6[ncounttemp6] * 100. < -58.91)
                        edgetofd = true;
                    if (debug_write)
                        cout << "tofdtemp6: " << xtofdtemp6[ncounttemp6] * 100 << ", " << ytofdtemp6[ncounttemp6] * 100
                             << ", " << qtofdtemp6[ncounttemp6] << endl;
                    ncounttemp6 += 1;
                }
                if (detector[i] == 7)
                {
                    xtofdtemp7[ncounttemp7] = xdet[i];
                    ytofdtemp7[ncounttemp7] = ydet[i];
                    qtofdtemp7[ncounttemp7] = qdet[i];
                    ttofdtemp7[ncounttemp7] = tdet[i];
                    ncounttemp7 += 1;
                }
                if (detector[i] == 8)
                {
                    xtofdtemp8[ncounttemp8] = xdet[i];
                    ytofdtemp8[ncounttemp8] = ydet[i];
                    qtofdtemp8[ncounttemp8] = qdet[i];
                    ttofdtemp8[ncounttemp8] = tdet[i];
                    if (xtofdtemp8[ncounttemp8] * 100. < -4.0 && xtofdtemp8[ncounttemp8] * 100. > -4.2)
                        bar21plane2 = true;
                    if (xtofdtemp8[ncounttemp8] * 100. > -5.48)
                        edgetofd = true;
                    if (debug_write)
                        cout << "tofdtemp8: " << xtofdtemp8[ncounttemp8] * 100 << ", " << ytofdtemp8[ncounttemp8] * 100
                             << ", " << qtofdtemp8[ncounttemp8] << endl;
                    ncounttemp8 += 1;
                }
                if (detector[i] == 9)
                {
                    xtofdtemp9[ncounttemp9] = xdet[i];
                    ytofdtemp9[ncounttemp9] = ydet[i];
                    qtofdtemp9[ncounttemp9] = qdet[i];
                    ttofdtemp9[ncounttemp9] = tdet[i];
                    if (xtofdtemp9[ncounttemp9] * 100. > 4.0 && xtofdtemp9[ncounttemp9] * 100. < 4.2)
                        bar24plane2 = true;
                    if (xtofdtemp9[ncounttemp9] * 100. > 57.54 || xtofdtemp9[ncounttemp9] * 100. < 5.48)
                        edgetofd = true;
                    ncounttemp9 += 1;
                }
                if (detector[i] == 10)
                {
                    xtofdtemp10[ncounttemp10] = xdet[i];
                    ytofdtemp10[ncounttemp10] = ydet[i];
                    qtofdtemp10[ncounttemp10] = qdet[i];
                    ttofdtemp10[ncounttemp10] = tdet[i];
                    if (debug_write)
                        cout << "tofdtemp10: " << xtofdtemp10[ncounttemp10] * 100 << ", "
                             << ytofdtemp10[ncounttemp10] * 100 << ", " << qtofdtemp10[ncounttemp10] << endl;
                    ncounttemp10 += 1;
                }
                if (detector[i] == 11)
                {
                    xtofdtemp11[ncounttemp11] = xdet[i];
                    ytofdtemp11[ncounttemp11] = ydet[i];
                    qtofdtemp11[ncounttemp11] = qdet[i];
                    ttofdtemp11[ncounttemp11] = tdet[i];
                    ncounttemp11 += 1;
                }
                ncount[detector[i]] += 1;
            }
            // Combine plane 2 (tofd) and plane 3 (tofi)
            Double_t qdiff = 100, xdiff = 100.;
            Bool_t checktof = false;

            Bool_t tofd11used[ncounttemp11];
            for (Int_t iw8 = 0; iw8 < ncounttemp11; iw8++)
            {
                tofd11used[iw8] = false;
            }

            if (debug_write)
            {
                if (ncounttemp11 > 0 || ncounttemp10 > 0)
                    checktof = true;
                if (ncounttemp11 < 1 && bar24plane2)
                    checktof = true;
            }

            if (checktof) //&& debug_write)
            {
                cout << "******************************************" << endl;

                for (Int_t iw8 = 0; iw8 < ncounttemp7; iw8++)
                {
                    cout << "Before change det7: " << xtofdtemp7[iw8] * 100. << "; " << ytofdtemp7[iw8] * 100. << " "
                         << qtofdtemp7[iw8] << endl;
                }
                for (Int_t iw8 = 0; iw8 < ncounttemp9; iw8++)
                {
                    cout << "Before change det9: " << xtofdtemp9[iw8] * 100. << "; " << ytofdtemp9[iw8] * 100. << " "
                         << qtofdtemp9[iw8] << endl;
                }
                for (Int_t iw8 = 0; iw8 < ncounttemp11; iw8++)
                {
                    cout << "Before change det11: " << xtofdtemp11[iw8] * 100. << "; " << ytofdtemp11[iw8] * 100. << " "
                         << qtofdtemp11[iw8] << endl;
                }
            }

            if (ncounttemp9 < 1 &&
                ncounttemp11 > 0) // trying to recover hits "lost" in bar 24 due to high counting rate
            {
                ncounttemp9 = 0;
                for (Int_t iw11 = 0; iw11 < ncounttemp11; iw11++)
                {
                    if (checktof)
                        cout << "First select in: " << xtofdtemp9[ncounttemp9] * 100 << ", " << xtofdtemp11[iw11] * 100
                             << endl;
                    xtofdtemp9[ncounttemp9] = xtofdtemp11[iw11];
                    ytofdtemp9[ncounttemp9] = ytofdtemp11[iw11];
                    qtofdtemp9[ncounttemp9] = qtofdtemp11[iw11];
                    ttofdtemp9[ncounttemp9] = ttofdtemp11[iw11];
                    ncounttemp9++;
                    if (checktof)
                        cout << "First select: " << ncounttemp11 << ", " << ncounttemp9 << ", "
                             << xtofdtemp9[ncounttemp9 - 1] * 100 << ", " << qtofdtemp9[ncounttemp9 - 1] << endl;
                }
            }
            else if (ncounttemp9 > 0)
            {
                for (Int_t iw9 = 0; iw9 < ncounttemp9; iw9++)
                {
                    if (xtofdtemp9[iw9] * 100. > 2.74 && xtofdtemp9[iw9] * 100. < 5.48) // plane 2, bar 24
                    {
                        if (ncounttemp11 == 0) // no corresponding hit in plane 3, ignore bar 24
                        {
                            if (checktof)
                                cout << "Second select in: " << xtofdtemp9[iw9] * 100 << ", event num: " << fNEvents
                                     << endl;
                            qtofdtemp9[iw9] = 0;
                            xtofdtemp9[iw9] = -1000;
                            ytofdtemp9[iw9] = -1000;

                            if (checktof)
                                cout << "Second select: " << ncounttemp11 << ", " << ncounttemp9 << ", "
                                     << xtofdtemp9[ncounttemp9] * 100 << ", " << qtofdtemp9[ncounttemp9] << endl;
                        }
                        for (Int_t iw11 = 0; iw11 < ncounttemp11; iw11++)
                        {
                            qdiff = std::abs(qtofdtemp11[iw11] - qtofdtemp9[iw9]);
                            xdiff = std::abs(xtofdtemp11[iw11] - xtofdtemp9[iw9]) * 100.;
                            if (checktof)
                                cout << "3rd sel, xdiff, qdiff: " << ncounttemp9 << "; " << ncounttemp11 << "; "
                                     << qdiff << ", " << xdiff << ", " << tofd11used[iw11] << endl;
                            if (checktof)
                                cout << "Third select first: " << iw9 << "; " << iw11 << ", " << xtofdtemp9[iw9] * 100
                                     << ", " << xtofdtemp11[iw11] * 100 << " " << qtofdtemp9[iw9] << ", "
                                     << qtofdtemp11[iw11] << endl;
                            if (!tofd11used[iw11])
                            {
                                if (xdiff < 1.38)
                                {
                                    if (checktof)
                                        cout << "Third select in: " << xtofdtemp9[iw9] * 100 << ", "
                                             << xtofdtemp11[iw11] * 100 << " " << qtofdtemp9[iw9] << ", "
                                             << qtofdtemp11[iw11] << endl;
                                    qtofdtemp9[iw9] = qtofdtemp11[iw11];
                                    xtofdtemp9[iw9] = xtofdtemp11[iw11];
                                    ytofdtemp9[iw9] = ytofdtemp11[iw11];
                                    tofd11used[iw11] = true;
                                    if (checktof)
                                        cout << "Third select: " << ncounttemp11 << ", " << ncounttemp9 << ", "
                                             << xtofdtemp9[iw9] * 100 << ", " << qtofdtemp9[iw9] << endl;
                                }
                                else
                                {
                                    if (checktof)
                                        cout << "Fourth select in: " << xtofdtemp9[iw9] * 100 << ", "
                                             << xtofdtemp11[iw11] * 100 << " " << qtofdtemp9[iw9] << ", "
                                             << qtofdtemp11[iw11] << endl;
                                    qtofdtemp9[iw9] = 0;
                                    xtofdtemp9[iw9] = -1000;
                                    ytofdtemp9[iw9] = -1000;
                                    tofd11used[iw11] = true;
                                    if (checktof)
                                        cout << "Fourth select: " << ncounttemp11 << ", " << ncounttemp9 << ", "
                                             << xtofdtemp9[ncounttemp9] * 100 << ", " << qtofdtemp9[ncounttemp9]
                                             << endl;
                                }
                            }
                        }
                    }

                    if (xtofdtemp9[iw9] * 100. > 5.48 && !bar24plane2) // plane 2, bar 25
                    {
                        for (Int_t iw11 = 0; iw11 < ncounttemp11; iw11++)
                        {
                            qdiff = std::abs(qtofdtemp11[iw11] - qtofdtemp9[iw9]);
                            xdiff = std::abs(xtofdtemp11[iw11] - xtofdtemp9[iw9]) * 100.;
                            if (checktof && debug_write)
                                cout << "Fifth select in: " << xtofdtemp9[iw9] * 100 << ", " << xtofdtemp11[iw11] * 100
                                     << " " << qtofdtemp9[iw9] << ", " << qtofdtemp11[iw11] << endl;
                            if (qdiff < 1.5 && xdiff < 1.38 && !tofd11used[iw11])
                            {
                                xtofdtemp9[iw9] = (xtofdtemp11[iw11] + xtofdtemp9[iw9]) / 2.;
                                tofd11used[iw11] = true;

                                if (checktof && debug_write)
                                    cout << "Fifth select: " << ncounttemp11 << ", " << ncounttemp9 << ", "
                                         << xtofdtemp9[iw9] * 100 << ", " << qtofdtemp9[iw9] << endl;
                            }
                            if (xtofdtemp11[iw11] * 100 < 5.54 && !tofd11used[iw11])
                            {
                                xtofdtemp9[ncounttemp9] = xtofdtemp11[iw11];
                                ytofdtemp9[ncounttemp9] = ytofdtemp11[iw11];
                                qtofdtemp9[ncounttemp9] = qtofdtemp11[iw11];
                                ttofdtemp9[ncounttemp9] = ttofdtemp11[iw11];
                                tofd11used[iw11] = true;

                                if (checktof && debug_write)
                                    cout << "Fifth select: " << ncounttemp11 << ", " << ncounttemp9 << ", "
                                         << xtofdtemp9[ncounttemp9] * 100 << ", " << qtofdtemp9[ncounttemp9] << endl;
                                ncounttemp9++;
                            }
                        }
                    }
                }
            }

            for (Int_t iw9 = 0; iw9 < ncounttemp9; iw9++)
            {
                if (xtofdtemp9[iw9] * 100. == 4.11)
                {
                    xtofdtemp9[iw9] = -1000;
                    ytofdtemp9[iw9] = -1000;
                    qtofdtemp9[iw9] = 0;
                }
            }

            Int_t imem = 0;
            Bool_t wo9[ncounttemp9];
            for (Int_t iw9 = 0; iw9 < ncounttemp9; iw9++)
            {
                wo9[iw9] = true;
                if (qtofdtemp9[iw9] < 1)
                {
                    wo9[iw9] = false;
                }
            }
            for (Int_t iw9 = 0; iw9 < ncounttemp9; iw9++)
            {
                if (checktof)
                    cout << "reordering 1: " << iw9 << "; " << wo9[iw9] << ", " << xtofdtemp9[iw9] << endl;
                if (wo9[iw9])
                {
                    xtofdtemp9[imem] = xtofdtemp9[iw9];
                    ytofdtemp9[imem] = ytofdtemp9[iw9];
                    qtofdtemp9[imem] = qtofdtemp9[iw9];
                    ttofdtemp9[imem] = ttofdtemp9[iw9];
                    if (checktof)
                        cout << "reordering 2: " << imem << ", " << iw9 << "; " << wo9[iw9] << ", "
                             << ", " << xtofdtemp9[imem] << ", " << xtofdtemp9[iw9] << endl;
                    imem++;
                }
            }

            ncounttemp9 = imem;

            if (checktof)
            {
                for (Int_t iw8 = 0; iw8 < ncounttemp7; iw8++)
                {
                    cout << "After change det7: " << xtofdtemp7[iw8] * 100. << "; " << ytofdtemp7[iw8] * 100. << " "
                         << qtofdtemp7[iw8] << endl;
                }
                for (Int_t iw8 = 0; iw8 < ncounttemp9; iw8++)
                {
                    cout << "After change det9: " << xtofdtemp9[iw8] * 100. << "; " << ytofdtemp9[iw8] * 100. << " "
                         << qtofdtemp9[iw8] << endl;
                }
                for (Int_t iw8 = 0; iw8 < ncounttemp11; iw8++)
                {
                    cout << "After change det11: " << xtofdtemp11[iw8] * 100. << "; " << ytofdtemp11[iw8] * 100. << " "
                         << qtofdtemp11[iw8] << endl;
                }
            }

            if (checktof && debug_write)
            {
                for (Int_t iw8 = 0; iw8 < ncounttemp6; iw8++)
                {
                    cout << "Before change det6: " << xtofdtemp6[iw8] * 100. << "; " << ytofdtemp6[iw8] * 100. << " "
                         << qtofdtemp6[iw8] << endl;
                }
                for (Int_t iw8 = 0; iw8 < ncounttemp8; iw8++)
                {
                    cout << "Before change det8: " << xtofdtemp8[iw8] * 100. << "; " << ytofdtemp8[iw8] * 100. << " "
                         << qtofdtemp8[iw8] << endl;
                }
                for (Int_t iw8 = 0; iw8 < ncounttemp10; iw8++)
                {
                    cout << "Before change det10: " << xtofdtemp10[iw8] * 100. << "; " << ytofdtemp10[iw8] * 100. << " "
                         << qtofdtemp10[iw8] << endl;
                }
            }
            Bool_t tof10used[ncounttemp10];
            for (Int_t i = 0; i < ncounttemp10; i++)
            {
                tof10used[i] = false;
            }
            if (ncounttemp10 > 0 && ncounttemp8 < 1)
            {
                for (Int_t iw10 = 0; iw10 < ncounttemp10; iw10++)
                {
                    xtofdtemp8[ncounttemp8] = xtofdtemp10[iw10];
                    ytofdtemp8[ncounttemp8] = ytofdtemp10[iw10];
                    qtofdtemp8[ncounttemp8] = qtofdtemp10[iw10];
                    ttofdtemp8[ncounttemp8] = ttofdtemp10[iw10];
                    ncounttemp8++;
                }
            }
            else if (ncounttemp8 > 0)
            {
                for (Int_t iw8 = 0; iw8 < ncounttemp8; iw8++)
                {
                    if (xtofdtemp8[iw8] * 100. < -2.74 && xtofdtemp8[iw8] * 100. > -5.48) // bar 21 plane 2
                    {
                        if (ncounttemp10 == ncounttemp8 - 1 && ncounttemp6 == ncounttemp8 - 1)
                            xtofdtemp8[iw8] = -3.71 / 100.;

                        for (Int_t iw10 = 0; iw10 < ncounttemp10; iw10++)
                        {
                            qdiff = std::abs(qtofdtemp10[iw10] - qtofdtemp8[iw8]);
                            xdiff = std::abs(xtofdtemp10[iw10] - xtofdtemp8[iw8]) * 100;

                            if (qdiff < 1 && xdiff < 1.8 && !tof10used[iw10])
                            {
                                xtofdtemp8[iw8] = (xtofdtemp8[iw8] + xtofdtemp10[iw10]) / 2.;
                                tof10used[iw10] = true;
                            }
                        }
                    }
                    if (xtofdtemp8[iw8] * 100. < -5.49 && xtofdtemp8[iw8] * 100. > -8.21) // bar 20 plane 2
                    {
                        for (Int_t iw10 = 0; iw10 < ncounttemp10; iw10++)
                        {
                            qdiff = std::abs(qtofdtemp10[iw10] - qtofdtemp8[iw8]);
                            xdiff = std::abs(xtofdtemp10[iw10] - xtofdtemp8[iw8]) * 100;

                            if (qdiff < 1 && xdiff < 1.38 && !tof10used[iw10])
                            {
                                xtofdtemp8[iw8] = (xtofdtemp8[iw8] + xtofdtemp10[iw10]) / 2.;
                                tof10used[iw10] = true;
                            }
                            if (!bar21plane2 && !tof10used[iw10])
                            {
                                xtofdtemp8[ncounttemp8] = xtofdtemp10[iw10];
                                ytofdtemp8[ncounttemp8] = ytofdtemp10[iw10];
                                qtofdtemp8[ncounttemp8] = qtofdtemp10[iw10];
                                ttofdtemp8[ncounttemp8] = ttofdtemp10[iw10];
                                tof10used[iw10] = true;
                                ncounttemp8++;
                            }
                        }
                    }
                }
            }

            if (checktof && debug_write)
            {
                for (Int_t iw8 = 0; iw8 < ncounttemp6; iw8++)
                {
                    cout << "After change det6: " << xtofdtemp6[iw8] * 100. << "; " << ytofdtemp6[iw8] * 100. << " "
                         << qtofdtemp6[iw8] << endl;
                }
                for (Int_t iw8 = 0; iw8 < ncounttemp8; iw8++)
                {
                    cout << "After change det8: " << xtofdtemp8[iw8] * 100. << "; " << ytofdtemp8[iw8] * 100. << " "
                         << qtofdtemp8[iw8] << endl;
                }
                for (Int_t iw8 = 0; iw8 < ncounttemp10; iw8++)
                {
                    cout << "After change det10: " << xtofdtemp10[iw8] * 100. << "; " << ytofdtemp10[iw8] * 100. << " "
                         << qtofdtemp10[iw8] << endl;
                }
            }

            for (Int_t iw8 = 0; iw8 < ncounttemp8; iw8++)
            {
                qsum8 += qtofdtemp8[iw8];
            }

            for (Int_t iw8 = 0; iw8 < ncounttemp6; iw8++)
            {
                qsum6 += qtofdtemp6[iw8];
            }

            for (Int_t iw8 = 0; iw8 < ncounttemp9; iw8++)
            {
                qsum9 += qtofdtemp9[iw8];
            }

            for (Int_t iw8 = 0; iw8 < ncounttemp7; iw8++)
            {
                qsum7 += qtofdtemp7[iw8];
            }

            Bool_t goodQ = false;
            Bool_t goodQtofd = false;
            Int_t isumtofd = ncounttemp6 + ncounttemp7 + ncounttemp8 + ncounttemp9;
            Double_t Qsum = qsum6 + qsum7 + qsum8 + qsum9;
            ncount[6] = ncounttemp6;
            ncount[7] = ncounttemp7;
            ncount[8] = ncounttemp8;
            ncount[9] = ncounttemp9;
            if (debug_write)
                cout << "Sum of Tofd charges after fibcuts: qsum6: " << qsum6 << ", qsum7: " << qsum7
                     << ", qsum8: " << qsum8 << ", qsum9: " << qsum9 << ", total: " << Qsum << endl;
            if (debug_write)
                cout << "Tofd counts after fibcuts: " << ncount[6] << "; " << ncount[7] << "; " << ncount[8] << ", "
                     << ncount[9] << ", total: " << isumtofd << endl;

            if (isumtofd == 2 && Qsum == 8.)
                goodQtofd = true;
            if (isumtofd == 3 && (Qsum == 10 || Qsum == 14) && edgetofd)
                goodQtofd = true;
            if (isumtofd == 4 && Qsum == 16)
                goodQtofd = true;
            if (isumtofd == 5 && (Qsum == 18 || Qsum == 22) && edgetofd)
                goodQtofd = true;
            if (isumtofd == 6 && (Qsum == 20 || Qsum == 28))
                goodQtofd = true;
            if (isumtofd == 8 && Qsum == 32)
                goodQtofd = true;

            Double_t Qsum30 = 0, Qsum31 = 0, Qsum32 = 0, Qsum33 = 0;
            for (Int_t i = 0; i < countdet; i++)
            {
                if (detector[i] == 2)
                    Qsum30 += qdet[i];
                if (detector[i] == 3)
                    Qsum31 += qdet[i];
                if (detector[i] == 4)
                    Qsum32 += qdet[i];
                if (detector[i] == 5)
                    Qsum33 += qdet[i];
            }

            if (goodQtofd && Qsum30 == Qsum32 && Qsum31 == Qsum33 && (Qsum30 + Qsum32 + Qsum31 + Qsum33) > 15)
                goodQ = true;

            if (debug_write && !goodQ)
                cout << "!goodQ: "
                     << "Qsum30: " << Qsum30 << ", QSum32: " << Qsum32 << ", QSum31: " << Qsum31
                     << ", QSum33: " << Qsum33 << ", QSumTod: " << Qsum << endl;

            Bool_t goodFi3032counts = false;
            Bool_t goodFi3133counts = false;
            Bool_t goodFi3X = false;
            // sum of hits in fi30+fi32 or fi31+fi33  must be even, and the total sum must be higher than zero:
            if (int((ncount[2] + ncount[4]) % 2) == 0)
                goodFi3032counts = true;
            if (int((ncount[3] + ncount[5]) % 2) == 0)
                goodFi3133counts = true;
            if (goodFi3032counts && goodFi3133counts && (ncount[2] + ncount[3] + ncount[4] + ncount[5]) > 0)
                goodFi3X = true;

            if (goodQ)
                Nevent_goodQ++;
            if (goodQtofd)
                Nevent_goodQtofd++;

            if (goodQtofd && ncount[0] < 2)
                Nevent_zerofi23a_corel++;
            if (goodQtofd && ncount[0] > 1)
                Nevent_notzerofi23a_corel++;

            if (goodQtofd && ncount[1] < 2)
                Nevent_zerofi23b_corel++;
            if (goodQtofd && ncount[1] > 1)
                Nevent_notzerofi23b_corel++;

            if (goodQtofd && !goodFi3032counts)
                Nevent_zerofi3032_corel++;
            if (goodQtofd && goodFi3032counts)
                Nevent_notzerofi3032_corel++;

            if (goodQtofd && !goodFi3133counts)
                Nevent_zerofi3133_corel++;
            if (goodQtofd && goodFi3133counts)
                Nevent_notzerofi3133_corel++;

            if (goodQtofd)
                fh_qsum_mult_fi23a->Fill(ncount[0], Qsum);
            if (goodQtofd)
                fh_qsum_mult_fi23b->Fill(ncount[1], Qsum);

            countdet_s = 0;

            Bool_t goodEvent = false;

            // Ideal case (He on one side, C on another)
            if (fidealData &&
                (ncount[9] == ncount[7] && ncount[6] == ncount[8] && ncount[8] == ncount[9] && ncount[2] > 0 &&
                 ncount[3] > 0 && ncount[4] > 0 && ncount[5] > 0 && ncount[0] > 1 && ncount[1] > 1 &&
                 ncount[0] > ncount[2] && ncount[0] > ncount[3] && ncount[1] > ncount[2] && ncount[1] > ncount[3] &&
                 ncount[2] == ncount[4] && ncount[3] == ncount[5]))
                goodEvent = true;

            // Normal case
            if (!fidealData && goodQ && goodFi3X && ncount[0] > 1 && ncount[1] > 1 &&
                ((cond1 && cond2 && cond3) || (cond4 && cond5 && cond6)))
                goodEvent = true;

            if (goodEvent)
            {

                fNeventselect += 1;

                if (1 == 1)
                {
                    for (Int_t i = 0; i < ncounttemp6; i++)
                    {
                        for (Int_t j = 0; j < ncounttemp8; j++)
                        {
                            if (qtofdtemp6[i] == qtofdtemp8[j])
                                fh_yTofd_vs_yTofd->Fill(ytofdtemp6[i] * 100., ytofdtemp8[j] * 100);
                        }
                    }
                    for (Int_t i = 0; i < ncounttemp7; i++)
                    {
                        for (Int_t j = 0; j < ncounttemp9; j++)
                        {
                            if (qtofdtemp7[i] == qtofdtemp9[j])
                                fh_yTofd_vs_yTofd->Fill(ytofdtemp7[i] * 100., ytofdtemp9[j] * 100);
                        }
                    }
                    for (Int_t j = 0; j < ncount[14]; j++)
                    {
                        for (Int_t i = 0; i < ncounttemp6; i++)
                        {
                            fh_Fi23aX_TofdX_bw1->Fill(xtofdtemp6[i] * 100., x23rtemp[j] * 100.);
                        }
                        for (Int_t i = 0; i < ncounttemp8; i++)
                        {
                            fh_Fi23aX_TofdX_bw2->Fill(xtofdtemp8[i] * 100., x23rtemp[j] * 100.);
                        }
                    }
                    for (Int_t j = 0; j < ncount[15]; j++)
                    {
                        for (Int_t i = 0; i < ncounttemp7; i++)
                        {
                            fh_Fi23aX_TofdX_bw1->Fill(xtofdtemp7[i] * 100., x23ltemp[j] * 100.);
                        }
                        for (Int_t i = 0; i < ncounttemp9; i++)
                        {
                            fh_Fi23aX_TofdX_bw2->Fill(xtofdtemp9[i] * 100., x23ltemp[j] * 100.);
                        }
                    }
                }

                for (Int_t i = 0; i < countdet; i++)
                {
                    if (detector[i] > 5 && detector[i] < 10)
                    {
                        xdet_s[countdet_s] = xdet[i];
                        ydet_s[countdet_s] = ydet[i];
                        zdet_s[countdet_s] = zdet[i];
                        qdet_s[countdet_s] = qdet[i];
                        tdet_s[countdet_s] = tdet[i];
                        detector_s[countdet_s] = detector[i];

                        if (debug_write)
                        {
                            //  cout << "tofd det:   " << detector[i] << ", count " << i << ", x: " << xdet[i]
                            //       << ", q: " << qdet[i] << endl;
                            cout << "tofd det_s: " << detector_s[countdet_s] << ", count " << countdet_s
                                 << ", x: " << xdet_s[countdet_s] * 100. << ", y: " << ydet_s[countdet_s] * 100.
                                 << ", q: " << qdet_s[countdet_s] << ", fNEvent: " << fNEvents << ", nHits: " << nHits
                                 << endl;
                        }
                        countdet_s++;
                    }
                    else if (detector[i] > 1 && detector[i] < 6)
                    {
                        // if (fSimu ||
                        //     (!fSimu && (((detector[i] == 2 || detector[i] == 4) && abs(qdet_s7 - qdet[i]) < 2.1) ||
                        //                ((detector[i] == 3 || detector[i] == 5) && abs(qdet_s6 - qdet[i]) < 2.1))))
                        if (1 == 1)
                        {
                            xdet_s[countdet_s] = xdet[i];
                            ydet_s[countdet_s] = ydet[i];
                            zdet_s[countdet_s] = zdet[i];
                            qdet_s[countdet_s] = qdet[i];
                            tdet_s[countdet_s] = tdet[i];
                            detector_s[countdet_s] = detector[i];

                            if (debug_write && 1 == 2)
                            {
                                //  cout << "fib3x  det:   " << detector[i] << ", count " << i << ", x: " << xdet[i]
                                //     << ", q: " << qdet[i] << endl;
                                cout << "fib3x  det_s: " << detector_s[countdet_s] << ", count " << countdet_s
                                     << ", x: " << xdet_s[countdet_s] * 100. << ", y: " << ydet_s[countdet_s] * 100.
                                     << ", q: " << qdet_s[countdet_s] << endl;
                            }
                            countdet_s++;
                        }
                    }
                    else if (detector[i] < 2)
                    {
                        xdet_s[countdet_s] = xdet[i];
                        ydet_s[countdet_s] = ydet[i];
                        zdet_s[countdet_s] = zdet[i];
                        qdet_s[countdet_s] = qdet[i];
                        tdet_s[countdet_s] = tdet[i];
                        detector_s[countdet_s] = detector[i];

                        if (debug_write && 1 == 2)
                        {
                            //  cout << "fib2x  det:   " << detector[i] << ", count " << i << ", x: " << xdet[i] * 100.
                            //     << ", y: " << ydet[i] * 100. << ", q: " << qdet[i] << endl;
                            cout << "fib2x  det_s: " << detector_s[countdet_s] << ", count " << countdet_s
                                 << ", x: " << xdet_s[countdet_s] * 100. << ", y: " << ydet_s[countdet_s] * 100.
                                 << ", q: " << qdet_s[countdet_s] << endl;
                        }
                        countdet_s++;
                    }
                }

                // Here write hit data of all detectors

                counter1++;
                // cout << "* counter1: " << counter1 << ", from Nevents: " << fNEvents << ", " << fNeventselect <<
                // endl;

                Int_t nfin[7] = { 0 }, nC = 0, nHe = 0;
                Double_t yfi23b[15], xfi23a[15], xfi30[10], yfi30[10], xfi31[10], yfi31[10], xfi32[10], yfi32[10],
                    xfi33[10], yfi33[10], ytofdC[15], xtofdC[15], qtofdC[15], ytofdHe[15], xtofdHe[15], qtofdHe[15],
                    ttofdHe[15], ttofdC[15], qfi30[10], qfi31[10], qfi32[10], qfi33[10], qfi23b[15], qfi23a[15],
                    tfi23b[15], tfi23a[15];

                Int_t count = 0;
                Double_t xbinavr;
                Double_t ffact;

                if (1 == 1)
                {
                    counterwo++;
                    //   cout << "* counter written out: " << counterwo << ", from Nevents: " << fNEvents << ", " <<
                    //   fNeventselect << endl;
                    SelectedPairs = true;

                    for (Int_t i = 0; i < countdet_s; i++)
                    {
                        //  if (debug_write)
                        //    cout << "#" << i << " Det: " << detector_s[i] << " x: " << xdet_s[i] * 100.
                        //     << " y: " << ydet_s[i] * 100. << " q: " << qdet_s[i] << " t: " << tdet_s[i] << endl;

                        if (detector_s[i] == 0)
                        {
                            xfi23a[nfin[0]] = xdet_s[i] * 100.;
                            qfi23a[nfin[0]] = qdet_s[i];
                            tfi23a[nfin[0]] = tdet_s[i];
                            nfin[0] += 1;
                            countdet_written[0]++;

                            // write fiber detector hits
                            new ((*fFi23aHitItems)[fNofFi23aHitItems++]) R3BFiberMAPMTHitData(
                                0, xdet_s[i] * 100., ydet_s[i] * 100., qdet_s[i], tdet_s[i], 0, 0, 0., 0, 0.);

                            if (debug_write && 1 == 2)
                                cout << "Fib23a write out: " << fNEvents << ", " << xdet_s[i] * 100 << "; "
                                     << ydet_s[i] * 100. << "; " << qdet_s[i] << ", " << tdet_s[i] << ", " << mult23a
                                     << endl;

                            fh_check_QvsX[detector_s[i]]->Fill(xdet_s[i] * 100, qdet_s[i]);
                            fh_check_TvsX[detector_s[i]]->Fill(xdet_s[i] * 100, tdet_s[i]);
                            fh_check_XvsY[detector_s[i]]->Fill(xdet_s[i] * 100, ydet_s[i] * 100);
                        }
                        if (detector_s[i] == 1)
                        {
                            countdet_written[1]++;

                            yfi23b[nfin[1]] = ydet_s[i] * 100.;
                            qfi23b[nfin[0]] = qdet_s[i];
                            tfi23b[nfin[0]] = tdet_s[i];
                            nfin[1] += 1;
                            // write fiber detector hits
                            new ((*fFi23bHitItems)[fNofFi23bHitItems++]) R3BFiberMAPMTHitData(
                                0, xdet_s[i] * 100., ydet_s[i] * 100., qdet_s[i], tdet_s[i], 0, 0, 0., 0, 0.);
                            if (debug_write && 1 == 2)
                                cout << "Fib23b write out: " << fNEvents << ", " << xdet_s[i] * 100 << "; "
                                     << ydet_s[i] * 100. << "; " << qdet_s[i] << ", " << tdet_s[i] << ", " << mult23b
                                     << endl;

                            fh_check_QvsX[detector_s[i]]->Fill(ydet_s[i] * 100, qdet_s[i]);
                            fh_check_TvsX[detector_s[i]]->Fill(ydet_s[i] * 100, tdet_s[i]);
                            fh_check_XvsY[detector_s[i]]->Fill(xdet_s[i] * 100, ydet_s[i] * 100);
                        }
                        if (detector_s[i] == 2)
                        {
                            countdet_written[2]++;

                            xfi30[nfin[2]] = xdet_s[i] * 100.;
                            yfi30[nfin[2]] = ydet_s[i] * 100.;
                            qfi30[nfin[2]] = qdet_s[i];
                            nfin[2] += 1;
                            // write fiber detector hits
                            new ((*fFi30HitItems)[fNofFi30HitItems++]) R3BFiberMAPMTHitData(
                                0, xdet_s[i] * 100., ydet_s[i] * 100., qdet_s[i], tdet_s[i], 0, 0, 0., 0, 0.);
                            if (debug_write && 1 == 2)
                                cout << "Fib30 write out: " << fNEvents << ", " << xdet_s[i] * 100 << "; "
                                     << ydet_s[i] * 100. << "; " << qdet_s[i] << ", " << mult30 << endl;
                            fh_check_QvsX[detector_s[i]]->Fill(xdet_s[i] * 100, qdet_s[i]);
                            fh_check_TvsX[detector_s[i]]->Fill(xdet_s[i] * 100, tdet_s[i]);
                            fh_check_XvsY[detector_s[i]]->Fill(xdet_s[i] * 100, ydet_s[i] * 100);
                        }
                        if (detector_s[i] == 3)
                        {
                            countdet_written[3]++;

                            xfi31[nfin[3]] = xdet_s[i] * 100.;
                            yfi31[nfin[3]] = ydet_s[i] * 100.;
                            qfi31[nfin[2]] = qdet_s[i];
                            nfin[3] += 1;
                            // write fiber detector hits
                            new ((*fFi31HitItems)[fNofFi31HitItems++]) R3BFiberMAPMTHitData(
                                0, xdet_s[i] * 100., ydet_s[i] * 100., qdet_s[i], tdet_s[i], 0, 0, 0., 0, 0.);
                            if (debug_write && 1 == 2)
                                cout << "Fib31 write out: " << fNEvents << ", " << xdet_s[i] * 100 << "; "
                                     << ydet_s[i] * 100. << "; " << qdet_s[i] << ", " << mult31 << endl;
                            fh_check_QvsX[detector_s[i]]->Fill(xdet_s[i] * 100, qdet_s[i]);
                            fh_check_TvsX[detector_s[i]]->Fill(xdet_s[i] * 100, tdet_s[i]);
                            fh_check_XvsY[detector_s[i]]->Fill(xdet_s[i] * 100, ydet_s[i] * 100);
                        }
                        if (detector_s[i] == 4)
                        {
                            countdet_written[4]++;

                            xfi32[nfin[4]] = xdet_s[i] * 100.;
                            yfi32[nfin[4]] = ydet_s[i] * 100.;
                            qfi32[nfin[2]] = qdet_s[i];
                            nfin[4] += 1;
                            // write fiber detector hits
                            new ((*fFi32HitItems)[fNofFi32HitItems++]) R3BFiberMAPMTHitData(
                                0, xdet_s[i] * 100., ydet_s[i] * 100., qdet_s[i], tdet_s[i], 0, 0, 0., 0, 0.);
                            if (debug_write && 1 == 2)
                                cout << "Fib32 write out: " << fNEvents << ", " << xdet_s[i] * 100 << "; "
                                     << ydet_s[i] * 100. << "; " << qdet_s[i] << ", " << mult32 << endl;
                            fh_check_QvsX[detector_s[i]]->Fill(xdet_s[i] * 100, qdet_s[i]);
                            fh_check_TvsX[detector_s[i]]->Fill(xdet_s[i] * 100, tdet_s[i]);
                            fh_check_XvsY[detector_s[i]]->Fill(xdet_s[i] * 100, ydet_s[i] * 100);
                        }
                        if (detector_s[i] == 5)
                        {
                            countdet_written[5]++;

                            // write fiber detector hits
                            xfi33[nfin[5]] = xdet_s[i] * 100.;
                            yfi33[nfin[5]] = ydet_s[i] * 100.;
                            qfi33[nfin[5]] = qdet_s[i];
                            nfin[5] += 1;
                            new ((*fFi33HitItems)[fNofFi33HitItems++]) R3BFiberMAPMTHitData(
                                0, xdet_s[i] * 100., ydet_s[i] * 100., qdet_s[i], tdet_s[i], 0, 0, 0., 0, 0.);
                            if (debug_write && 1 == 2)
                                cout << "Fib33 write out: " << fNEvents << ", " << xdet_s[i] * 100 << "; "
                                     << ydet_s[i] * 100. << "; " << qdet_s[i] << ", " << mult33 << endl;
                            fh_check_QvsX[detector_s[i]]->Fill(xdet_s[i] * 100, qdet_s[i]);
                            fh_check_TvsX[detector_s[i]]->Fill(xdet_s[i] * 100, tdet_s[i]);
                            fh_check_XvsY[detector_s[i]]->Fill(xdet_s[i] * 100, ydet_s[i] * 100);
                        }
                        if ((detector_s[i] == 6 || detector_s[i] == 8 || detector_s[i] == 7 || detector_s[i] == 9) &&
                            ftrackerType == 0) // TofD data written out for each plane separately
                        {
                            countdet_written[detector_s[i]]++;

                            Int_t ipl;
                            if (detector_s[i] == 6 || detector_s[i] == 7)
                                ipl = 1;
                            else
                                ipl = 2;

                            new ((*fTofdHitItems)[fNofTofdHitItems++]) R3BTofdHitData(
                                tdet_s[i], xdet_s[i] * 100., ydet_s[i] * 100., qdet_s[i], -5, qdet_s[i], ipl, 1, 0);
                            if (debug_write)
                                cout << "ToFD write out: " << detector_s[i] << ", " << ipl << "; " << xdet_s[i] * 100
                                     << "; " << ydet_s[i] * 100. << "; " << qdet_s[i] << endl;

                            fh_check_QvsX[6]->Fill(xdet_s[i] * 100, qdet_s[i]);
                            fh_check_TvsX[6]->Fill(xdet_s[i] * 100, tdet_s[i]);
                            fh_check_XvsY[6]->Fill(xdet_s[i] * 100, ydet_s[i] * 100);
                        }
                    }
                }

                fh_mult_Fib_ac[fi23a]->Fill(nfin[0]);
                fh_mult_Fib_ac[fi23b]->Fill(nfin[1]);
                fh_mult_Fib_ac[fi30]->Fill(nfin[2]);
                fh_mult_Fib_ac[fi31]->Fill(nfin[3]);
                fh_mult_Fib_ac[fi32]->Fill(nfin[4]);
                fh_mult_Fib_ac[fi33]->Fill(nfin[5]);

                Int_t isumHe = 0, isumC = 0;
                if (ftrackerType == 1)
                {
                    nfin[6] = 0;
                    Double_t qleft[20] = { 0 }, qright[20] = { 0 }, xleft[20], yleft[20], tleft[20], xright[20],
                             yright[20], tright[20];
                    qdiff = 100.;
                    xdiff = 100.;
                    Int_t iright = 0, ileft = 0;
                    Bool_t imemleft[20], imemright[20];
                    for (Int_t i = 0; i < 20; i++)
                    {
                        imemleft[i] = false;
                        imemright[i] = false;
                        xleft[i] = -1000;
                        yleft[i] = -1000;
                        tleft[i] = -1000;
                        xright[i] = -1000;
                        yright[i] = -1000;
                        tright[i] = -1000;
                    }

                    countdet_written[6] += ncounttemp6;
                    countdet_written[7] += ncounttemp7;
                    countdet_written[8] += ncounttemp8;
                    countdet_written[9] += ncounttemp9;

                    if (ncounttemp8 == ncounttemp6)
                    {
                        for (Int_t iw8 = 0; iw8 < ncounttemp8; iw8++)
                        {
                            for (Int_t iw6 = 0; iw6 < ncounttemp6; iw6++)
                            {
                                qdiff = std::abs(qtofdtemp8[iw8] - qtofdtemp6[iw6]);
                                xdiff = std::abs(xtofdtemp8[iw8] - xtofdtemp6[iw6]) * 100;
                                if (qdiff < 0.5 && xdiff < 1.38)
                                {
                                    qright[iright] = (qtofdtemp8[iw8] + qtofdtemp6[iw6]) / 2.;
                                    // if(xtofdtemp8[iw8]*100. < -5.62)
                                    xright[iright] = (xtofdtemp8[iw8] + xtofdtemp6[iw6]) / 2.;
                                    //	else
                                    //	xright[iright] = xtofdtemp8[iw8];
                                    yright[iright] = (ytofdtemp8[iw8] + ytofdtemp6[iw6]) / 2.;
                                    tright[iright] = (ttofdtemp8[iw8] + ttofdtemp6[iw6]) / 2.;
                                    iright++;
                                }
                            }
                        }
                    }
                    if (ncounttemp8 > ncounttemp6)
                    {
                        for (Int_t iw8 = 0; iw8 < ncounttemp8; iw8++)
                        {
                            for (Int_t iw6 = 0; iw6 < ncounttemp6; iw6++)
                            {
                                qdiff = std::abs(qtofdtemp8[iw8] - qtofdtemp6[iw6]);
                                xdiff = std::abs(xtofdtemp8[iw8] - xtofdtemp6[iw6]) * 100;

                                if (qdiff < 0.5 && xdiff < 1.38)
                                {
                                    qright[iright] = (qtofdtemp8[iw8] + qtofdtemp6[iw6]) / 2.;
                                    //  if(xtofdtemp8[iw8]*100. < -5.62)
                                    xright[iright] = (xtofdtemp8[iw8] + xtofdtemp6[iw6]) / 2.;
                                    // else
                                    //	xright[iright] = xtofdtemp8[iw8];
                                    yright[iright] = (ytofdtemp8[iw8] + ytofdtemp6[iw6]) / 2.;
                                    tright[iright] = (ttofdtemp8[iw8] + ttofdtemp6[iw6]) / 2.;
                                    imemright[iw8] = true;
                                    iright++;
                                }
                            }
                        }
                        for (Int_t iw8 = 0; iw8 < ncounttemp8; iw8++)
                        {
                            if (imemright[iw8] == false)
                            {
                                qright[iright] = qtofdtemp8[iw8];
                                xright[iright] = xtofdtemp8[iw8];
                                yright[iright] = ytofdtemp8[iw8];
                                tright[iright] = ttofdtemp8[iw8];
                                iright++;
                            }
                        }
                    }
                    if (ncounttemp8 < ncounttemp6)
                    {
                        for (Int_t iw6 = 0; iw6 < ncounttemp6; iw6++)
                        {
                            for (Int_t iw8 = 0; iw8 < ncounttemp8; iw8++)
                            {
                                qdiff = std::abs(qtofdtemp8[iw8] - qtofdtemp6[iw6]);
                                xdiff = std::abs(xtofdtemp8[iw8] - xtofdtemp6[iw6]) * 100;
                                if (qdiff < 0.5 && xdiff < 1.38)
                                {
                                    qright[iright] = (qtofdtemp8[iw8] + qtofdtemp6[iw6]) / 2.;
                                    // if(xtofdtemp8[iw8]*100. < -5.62)
                                    xright[iright] = (xtofdtemp8[iw8] + xtofdtemp6[iw6]) / 2.;
                                    //	else
                                    //	xright[iright] = xtofdtemp8[iw8];
                                    yright[iright] = (ytofdtemp8[iw8] + ytofdtemp6[iw6]) / 2.;
                                    tright[iright] = (ttofdtemp8[iw8] + ttofdtemp6[iw6]) / 2.;
                                    imemright[iw6] = true;
                                    iright++;
                                }
                            }
                        }
                        for (Int_t iw6 = 0; iw6 < ncounttemp6; iw6++)
                        {
                            if (imemright[iw6] == false)
                            {
                                qright[iright] = qtofdtemp6[iw6];
                                xright[iright] = xtofdtemp6[iw6];
                                yright[iright] = ytofdtemp6[iw6];
                                tright[iright] = ttofdtemp6[iw6];
                                iright++;
                            }
                        }
                    }

                    if (debug_write)
                    {
                        cout << "**** FINAL VALUES ****" << endl;
                        for (Int_t i = 0; i < iright; i++)
                        {
                            cout << "Right side values: " << xright[i] * 100 << ", " << yright[i] * 100 << ", "
                                 << qright[i] << endl;
                        }
                    }

                    if (ncounttemp9 == ncounttemp7)
                    {
                        for (Int_t iw9 = 0; iw9 < ncounttemp9; iw9++)
                        {
                            for (Int_t iw7 = 0; iw7 < ncounttemp7; iw7++)
                            {
                                qdiff = std::abs(qtofdtemp9[iw9] - qtofdtemp7[iw7]);
                                xdiff = std::abs(xtofdtemp9[iw9] - xtofdtemp7[iw7]) * 100;
                                if (qdiff < 0.5 && xdiff < 1.38)
                                {
                                    qleft[ileft] = (qtofdtemp9[iw9] + qtofdtemp7[iw7]) / 2.;
                                    if (xtofdtemp9[iw9] * 100 > 5.83)
                                        xleft[ileft] = (xtofdtemp9[iw9] + xtofdtemp7[iw7]) / 2.;
                                    else
                                        xleft[ileft] = xtofdtemp9[iw9];
                                    yleft[ileft] = (ytofdtemp9[iw9] + ytofdtemp7[iw7]) / 2.;
                                    tleft[ileft] = (ttofdtemp9[iw9] + ttofdtemp7[iw7]) / 2.;
                                    ileft++;
                                }
                            }
                        }
                    }

                    if (ncounttemp9 > ncounttemp7)
                    {
                        for (Int_t iw9 = 0; iw9 < ncounttemp9; iw9++)
                        {
                            for (Int_t iw7 = 0; iw7 < ncounttemp7; iw7++)
                            {
                                qdiff = std::abs(qtofdtemp9[iw9] - qtofdtemp7[iw7]);
                                xdiff = std::abs(xtofdtemp9[iw9] - xtofdtemp7[iw7]) * 100;
                                if (qdiff < 0.5 && xdiff < 1.38)
                                {
                                    qleft[ileft] = (qtofdtemp9[iw9] + qtofdtemp7[iw7]) / 2.;
                                    if (xtofdtemp9[iw9] * 100 > 5.83)
                                        xleft[ileft] = (xtofdtemp9[iw9] + xtofdtemp7[iw7]) / 2.;
                                    else
                                        xleft[ileft] = xtofdtemp9[iw9];
                                    yleft[ileft] = (ytofdtemp9[iw9] + ytofdtemp7[iw7]) / 2.;
                                    tleft[ileft] = (ttofdtemp9[iw9] + ttofdtemp7[iw7]) / 2.;
                                    imemleft[iw9] = true;
                                    ileft++;
                                }
                            }
                        }
                        for (Int_t iw9 = 0; iw9 < ncounttemp9; iw9++)
                        {
                            if (imemleft[iw9] == false)
                            {
                                qleft[ileft] = qtofdtemp9[iw9];
                                xleft[ileft] = xtofdtemp9[iw9];
                                yleft[ileft] = ytofdtemp9[iw9];
                                tleft[ileft] = ttofdtemp9[iw9];
                                ileft++;
                            }
                        }
                    }
                    if (ncounttemp9 < ncounttemp7)
                    {
                        for (Int_t iw7 = 0; iw7 < ncounttemp7; iw7++)
                        {
                            for (Int_t iw9 = 0; iw9 < ncounttemp9; iw9++)
                            {
                                qdiff = std::abs(qtofdtemp9[iw9] - qtofdtemp7[iw7]);
                                xdiff = std::abs(xtofdtemp9[iw9] - xtofdtemp7[iw7]) * 100;
                                if (qdiff < 0.5 && xdiff < 1.38)
                                {
                                    qleft[ileft] = (qtofdtemp9[iw9] + qtofdtemp7[iw7]) / 2.;
                                    if (xtofdtemp9[iw9] * 100 > 5.83)
                                        xleft[ileft] = (xtofdtemp9[iw9] + xtofdtemp7[iw7]) / 2.;
                                    else
                                        xleft[ileft] = xtofdtemp9[iw9];
                                    yleft[ileft] = (ytofdtemp9[iw9] + ytofdtemp7[iw7]) / 2.;
                                    tleft[ileft] = (ttofdtemp9[iw9] + ttofdtemp7[iw7]) / 2.;
                                    imemleft[iw7] = true;
                                    ileft++;
                                }
                            }
                        }
                        for (Int_t iw7 = 0; iw7 < ncounttemp7; iw7++)
                        {
                            if (imemleft[iw7] == false)
                            {
                                qleft[ileft] = qtofdtemp7[iw7];
                                xleft[ileft] = xtofdtemp7[iw7];
                                yleft[ileft] = ytofdtemp7[iw7];
                                tleft[ileft] = ttofdtemp7[iw7];
                                ileft++;
                            }
                        }
                    }

                    if (debug_write)
                    {
                        cout << "***** FINAL VALUES ********" << endl;
                        for (Int_t i = 0; i < ileft; i++)
                        {
                            cout << "Left side values: " << xleft[i] * 100 << ", " << yleft[i] * 100 << ", " << qleft[i]
                                 << endl;
                        }
                    }

                    for (Int_t l = 0; l < ileft; l++)
                    {

                        if (qleft[l] == 2)
                        {
                            xtofdHe[isumHe] = xleft[l] * 100.;
                            ytofdHe[isumHe] = yleft[l] * 100.;
                            qtofdHe[isumHe] = 2;
                            ttofdHe[isumHe] = tleft[l];
                            isumHe++;
                        }
                        if (qleft[l] == 6)
                        {
                            xtofdC[isumC] = xleft[l] * 100.;
                            ytofdC[isumC] = yleft[l] * 100.;
                            qtofdC[isumC] = 6;
                            ttofdC[isumC] = tleft[l];
                            isumC++;
                        }

                        new ((*fTofdHitItems)[fNofTofdHitItems++]) R3BTofdHitData(
                            tleft[l], xleft[l] * 100., yleft[l] * 100., qleft[l], -5., qleft[l], 10, 1, 0);
                        if (debug_write)
                            cout << "Tofd write out: " << fNEvents << ", " << xleft[l] * 100 << "; " << yleft[l] * 100.
                                 << "; " << qleft[l] << endl;
                        fh_check_QvsX[6]->Fill(xleft[l] * 100, qleft[l]);
                        fh_check_TvsX[6]->Fill(xleft[l] * 100, tleft[l]);
                        fh_check_XvsY[6]->Fill(xleft[l] * 100, yleft[l] * 100);
                    }

                    for (Int_t l = 0; l < iright; l++)
                    {

                        if (qright[l] == 2)
                        {
                            xtofdHe[isumHe] = xright[l] * 100.;
                            ytofdHe[isumHe] = yright[l] * 100.;
                            qtofdHe[isumHe] = 2;
                            ttofdHe[isumHe] = tright[l];
                            isumHe++;
                        }
                        if (qright[l] == 6)
                        {
                            xtofdC[isumC] = xright[l] * 100.;
                            ytofdC[isumC] = yright[l] * 100.;
                            qtofdC[isumC] = 6;
                            ttofdC[isumC] = tright[l];
                            isumC++;
                        }

                        new ((*fTofdHitItems)[fNofTofdHitItems++]) R3BTofdHitData(
                            tright[l], xright[l] * 100., yright[l] * 100., qright[l], -5., qright[l], 10, 1, 0);
                        if (debug_write)
                            cout << "Tofd write out: " << fNEvents << ", " << xright[l] * 100 << "; "
                                 << yright[l] * 100. << "; " << qright[l] << endl;
                        fh_check_QvsX[6]->Fill(xright[l] * 100, qright[l]);
                        fh_check_TvsX[6]->Fill(xright[l] * 100, tright[l]);
                        fh_check_XvsY[6]->Fill(xright[l] * 100, yright[l] * 100);
                    }

                    for (Int_t i = 0; i < isumHe; i++)
                    {
                        for (Int_t k = 0; k < isumC; k++)
                        {
                            fh_xTofdHe_vs_xTofdC->Fill(xtofdHe[i], xtofdC[k]);
                            fh_yTofdHe_vs_yTofdC->Fill(ytofdHe[i], ytofdC[k]);
                            if ((xtofdHe[i] > 0 && xtofdC[k] > 0) || (xtofdHe[i] < 0 && xtofdC[k] < 0))
                                ncorellx += 1;
                            if ((xtofdHe[i] < 0 && xtofdC[k] > 0) || (xtofdHe[i] > 0 && xtofdC[k] < 0))
                                nanticorellx += 1;
                            if ((ytofdHe[i] > 0 && ytofdC[k] > 0) || (ytofdHe[i] < 0 && ytofdC[k] < 0))
                                ncorelly += 1;
                            if ((ytofdHe[i] < 0 && ytofdC[k] > 0) || (ytofdHe[i] > 0 && ytofdC[k] < 0))
                                nanticorelly += 1;
                        }
                    }
                }

                for (int i = 0; i < nfin[0]; i++)
                {
                    for (int j = 0; j < nfin[1]; j++)
                    {
                        fh_xy_target->Fill(xfi23a[i], yfi23b[j]);
                    }
                }

                for (int j = 0; j < isumHe; j++)
                {
                    for (int i = 0; i < nfin[0]; i++)
                    {
                        fh_Fi23aX_TofdX_He->Fill(xtofdHe[j], xfi23a[i]);
                        fh_Fi23aToF_Q[nfin[0]]->Fill(qfi23a[i], ttofdHe[j] - tfi23a[i]);
                    }
                    for (int i = 0; i < nfin[1]; i++)
                    {
                        fh_Fi23bY_TofdY_He->Fill(ytofdHe[j], yfi23b[i]);
                        fh_Fi23bToF_Q[nfin[1]]->Fill(qfi23b[i], ttofdHe[j] - tfi23b[i]);
                    }
                    for (int i = 0; i < nfin[2]; i++)
                    {
                        if (qfi30[i] == 2)
                            fh_Fi301X_TofdX_He->Fill(xtofdHe[j], xfi30[i]);

                        if (qfi30[i] == 2)
                            fh_Fi301Y_TofdY_He->Fill(ytofdHe[j], yfi30[i]);
                    }
                    for (int i = 0; i < nfin[3]; i++)
                    {
                        if (qfi31[i] == 2)
                            fh_Fi301X_TofdX_He->Fill(xtofdHe[j], xfi31[i]);

                        if (qfi31[i] == 2)
                            fh_Fi301Y_TofdY_He->Fill(ytofdHe[j], yfi31[i]);
                    }
                    for (int i = 0; i < nfin[4]; i++)
                    {
                        if (qfi32[i] == 2)
                            fh_Fi323X_TofdX_He->Fill(xtofdHe[j], xfi32[i]);

                        if (qfi32[i] == 2)
                            fh_Fi323Y_TofdY_He->Fill(ytofdHe[j], yfi32[i]);
                    }
                    for (int i = 0; i < nfin[5]; i++)
                    {
                        if (qfi33[i] == 2)
                            fh_Fi323X_TofdX_He->Fill(xtofdHe[j], xfi33[i]);

                        if (qfi33[i] == 2)
                            fh_Fi323Y_TofdY_He->Fill(ytofdHe[j], yfi33[i]);
                    }
                }
                for (int j = 0; j < isumC; j++)
                {
                    for (int i = 0; i < nfin[0]; i++)
                    {
                        fh_Fi23aX_TofdX_C->Fill(xtofdC[j], xfi23a[i]);
                        fh_Fi23aToF_Q[nfin[0]]->Fill(qfi23a[i], ttofdC[j] - tfi23a[i]);
                    }
                    for (int i = 0; i < nfin[1]; i++)
                    {
                        fh_Fi23bY_TofdY_C->Fill(ytofdC[j], yfi23b[i]);
                        fh_Fi23bToF_Q[nfin[1]]->Fill(qfi23b[i], ttofdC[j] - tfi23b[i]);
                    }
                    for (int i = 0; i < nfin[2]; i++)
                    {
                        if (qfi30[i] == 6)
                            fh_Fi301X_TofdX_C->Fill(xtofdC[j], xfi30[i]);

                        if (qfi30[i] == 6)
                            fh_Fi301Y_TofdY_C->Fill(ytofdC[j], yfi30[i]);
                    }
                    for (int i = 0; i < nfin[3]; i++)
                    {
                        if (qfi31[i] == 6)
                            fh_Fi301X_TofdX_C->Fill(xtofdC[j], xfi31[i]);

                        if (qfi31[i] == 6)
                            fh_Fi301Y_TofdY_C->Fill(ytofdC[j], yfi31[i]);
                    }
                    for (int i = 0; i < nfin[4]; i++)
                    {
                        if (qfi32[i] == 6)
                            fh_Fi323X_TofdX_C->Fill(xtofdC[j], xfi32[i]);

                        if (qfi32[i] == 6)
                            fh_Fi323Y_TofdY_C->Fill(ytofdC[j], yfi32[i]);
                    }
                    for (int i = 0; i < nfin[5]; i++)
                    {
                        if (qfi33[i] == 6)
                            fh_Fi323X_TofdX_C->Fill(xtofdC[j], xfi33[i]);

                        if (qfi33[i] == 6)
                            fh_Fi323Y_TofdY_C->Fill(ytofdC[j], yfi33[i]);
                    }
                }
            }
        }

        for (int i = 0; i < n_det; i++)
        {
            x[i] = -1000.;
            y[i] = -1000.;
            z[i] = -1000.;
            q[i] = -1000.;
            t[i] = -1000.;

            x1[i] = -1000.;
            y1[i] = -1000.;
            z1[i] = -1000.;
            q1[i] = -1000.;
            t1[i] = -1000.;

            x2[i] = -1000.;
            y2[i] = -1000.;
            z2[i] = -1000.;
            q2[i] = -1000.;
            t2[i] = -1000.;

            xMax[i] = -1000.;
            yMax[i] = -1000.;
            zMax[i] = -1000.;
            qMax[i] = -1000.;
            tMax[i] = -1000.;
        }

        if (fHitItems.at(DET_CALIFA) && SelectedPairs && !fSimu)
        {
            // CALIFA
            auto detCalifa = fHitItems.at(DET_CALIFA);
            nHitsCalifa = detCalifa->GetEntriesFast();
            for (Int_t ihitC = 0; ihitC < nHitsCalifa; ihitC++)
            {
                R3BCalifaClusterData* hitCalifa = (R3BCalifaClusterData*)detCalifa->At(ihitC);
                if (!hitCalifa)
                    continue;
                Double_t theta = hitCalifa->GetTheta();
                Double_t Energy = hitCalifa->GetEnergy();
                Double_t Energy_dc = hitCalifa->GetEnergy() * GAMMA * (1 - BETA * TMath::Cos(theta));

                if (IS_NAN(Energy))
                    fh_califa_overflow->Fill(hitCalifa->GetClusterType(), hitCalifa->GetMotherCrystal());
                else
                    fh_califa_good->Fill(hitCalifa->GetClusterType(), hitCalifa->GetMotherCrystal());

                if (hitCalifa->GetClusterType() == 1 && !(IS_NAN(Energy)))
                {
                    ULong64_t timeCalifa = hitCalifa->GetTime();
                    Double_t timerelCalifa = (double)(timeCalifa - time);
                    if (hitCalifa->GetMotherCrystal() > 927 && hitCalifa->GetMotherCrystal() < 1953) // gammas in barrel
                    {
                        fh_califa_tofd_barrel->Fill(timerelCalifa, tStart);
                        fh_califa_energy_select_barrel->Fill(hitCalifa->GetCrystalList().size(), Energy_dc);
                    }
                    if (hitCalifa->GetMotherCrystal() > 1952 &&
                        hitCalifa->GetMotherCrystal() < 2433) // gammas in barrel
                    {
                        fh_califa_tofd_iphos->Fill(timerelCalifa, tStart);
                        fh_califa_energy_select_iphos->Fill(hitCalifa->GetCrystalList().size(), Energy_dc);
                    }

                    new ((*fCalifaHitItems)[fNofCalifaHitItems++]) R3BCalifaClusterData(hitCalifa->GetCrystalList(),
                                                                                        hitCalifa->GetEnergyList(),
                                                                                        hitCalifa->GetTimeList(),
                                                                                        hitCalifa->GetEnergy(),
                                                                                        hitCalifa->GetNf(),
                                                                                        hitCalifa->GetNs(),
                                                                                        hitCalifa->GetTheta(),
                                                                                        hitCalifa->GetPhi(),
                                                                                        hitCalifa->GetTime(),
                                                                                        hitCalifa->GetClusterType());
                }
            }
        }

    } // end ToFD loop

    if (multTofd > 0)
        fh_tofd_mult_ac->Fill(multTofd);
}
void R3BPreTrackS494::FinishEvent()
{
    fNofTofdHitItems = 0;
    fTofdHitItems->Clear();
    fNofFi23aHitItems = 0;
    fFi23aHitItems->Clear();
    fNofFi23bHitItems = 0;
    fFi23bHitItems->Clear();
    fNofFi30HitItems = 0;
    fFi30HitItems->Clear();
    fNofFi31HitItems = 0;
    fFi31HitItems->Clear();
    fNofFi32HitItems = 0;
    fFi32HitItems->Clear();
    fNofFi33HitItems = 0;
    fFi33HitItems->Clear();
    if (fHitItems.at(DET_CALIFA))
    {
        fNofCalifaHitItems = 0;
        fCalifaHitItems->Clear();
    }
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

void R3BPreTrackS494::FinishTask()
{

    if (!fPairs && fAverage)
    {
        Double_t ax = 0.;
        Double_t ay = 0.;
        Double_t aq = 0.;
        Double_t at = 0.;

        for (Int_t i = 0; i < 10; i++)
        {
            // if (debug_write)
            cout << "output :"
                 << " #" << i << " Det: " << i << " x: " << xdet_sum[i] * 100. / nsum[i]
                 << " y: " << ydet_sum[i] * 100. / nsum[i] << " q: " << qdet_sum[i] / nsum[i] << endl;

            xdet_sum[i] = xdet_sum[i] / nsum[i];
            ydet_sum[i] = ydet_sum[i] / nsum[i];
            qdet_sum[i] = qdet_sum[i] / nsum[i];
            tdet_sum[i] = tdet_sum[i] / nsum[i];

            if (i == 0)
            {
                // write fiber detector hits
                new ((*fFi23aHitItems)[fNofFi23aHitItems++]) R3BFiberMAPMTHitData(
                    0, xdet_sum[i] * 100., ydet_sum[i] * 100., qdet_sum[i], tdet_sum[i], 0, 0, 0., 0, 0.);
                fh_check_QvsX[i]->Fill(xdet_sum[i] * 100, qdet_sum[i]);
                fh_check_TvsX[i]->Fill(xdet_sum[i] * 100, tdet_sum[i]);
                fh_check_XvsY[i]->Fill(xdet_sum[i] * 100, ydet_sum[i] * 100);
            }
            if (i == 1)
            {
                // write fiber detector hits
                new ((*fFi23bHitItems)[fNofFi23bHitItems++]) R3BFiberMAPMTHitData(
                    0, xdet_sum[i] * 100., ydet_sum[i] * 100., qdet_sum[i], tdet_sum[i], 0, 0, 0., 0, 0.);
                fh_check_QvsX[i]->Fill(xdet_sum[i] * 100, qdet_sum[i]);
                fh_check_TvsX[i]->Fill(xdet_sum[i] * 100, tdet_sum[i]);
                fh_check_XvsY[i]->Fill(xdet_sum[i] * 100, ydet_sum[i] * 100);
            }
            if (i == 2 && abs(fB) < 1710.)
            {
                // write fiber detector hits
                new ((*fFi30HitItems)[fNofFi30HitItems++]) R3BFiberMAPMTHitData(
                    0, xdet_sum[i] * 100., ydet_sum[i] * 100., qdet_sum[i], tdet_sum[i], 0, 0, 0., 0, 0.);
                fh_check_QvsX[i]->Fill(xdet_sum[i] * 100, qdet_sum[i]);
                fh_check_TvsX[i]->Fill(xdet_sum[i] * 100, tdet_sum[i]);
                fh_check_XvsY[i]->Fill(xdet_sum[i] * 100, ydet_sum[i] * 100);
            }
            if (i == 3 && abs(fB) > 1710.)
            {
                // write fiber detector hits
                new ((*fFi31HitItems)[fNofFi31HitItems++]) R3BFiberMAPMTHitData(
                    0, xdet_sum[i] * 100., ydet_sum[i] * 100., qdet_sum[i], tdet_sum[i], 0, 0, 0., 0, 0.);
                fh_check_QvsX[i]->Fill(xdet_sum[i] * 100, qdet_sum[i]);
                fh_check_TvsX[i]->Fill(xdet_sum[i] * 100, tdet_sum[i]);
                fh_check_XvsY[i]->Fill(xdet_sum[i] * 100, ydet_sum[i] * 100);
            }
            if (i == 4 && abs(fB) < 1710.)
            {
                // write fiber detector hits
                new ((*fFi32HitItems)[fNofFi32HitItems++]) R3BFiberMAPMTHitData(
                    0, xdet_sum[i] * 100., ydet_sum[i] * 100., qdet_sum[i], tdet_sum[i], 0, 0, 0., 0, 0.);
                fh_check_QvsX[i]->Fill(xdet_sum[i] * 100, qdet_sum[i]);
                fh_check_TvsX[i]->Fill(xdet_sum[i] * 100, tdet_sum[i]);
                fh_check_XvsY[i]->Fill(xdet_sum[i] * 100, ydet_sum[i] * 100);
            }
            if (i == 5 && abs(fB) > 1710.)
            {
                // write fiber detector hits
                //  if(isumdet5 == 0){
                new ((*fFi33HitItems)[fNofFi33HitItems++]) R3BFiberMAPMTHitData(
                    0, xdet_sum[i] * 100., ydet_sum[i] * 100., qdet_sum[i], tdet_sum[i], 0, 0, 0., 0, 0.);
                fh_check_QvsX[i]->Fill(xdet_sum[i] * 100, qdet_sum[i]);
                fh_check_TvsX[i]->Fill(xdet_sum[i] * 100, tdet_sum[i]);
                fh_check_XvsY[i]->Fill(xdet_sum[i] * 100, ydet_sum[i] * 100);
                //	}
                //      isumdet5 += 1;
            }
            if (i == 6 || i == 8 || i == 7 || i == 9)
            {
                if (ftrackerType == 0) // TofD data written out for each plane separately
                {
                    Int_t ipl;
                    if (i == 6 || i == 7)
                        ipl = 1;
                    else
                        ipl = 2;

                    if (((i == 6 || i == 8) && abs(fB) > 1710.) || ((i == 7 || i == 9) && abs(fB) < 1710.))
                    {
                        new ((*fTofdHitItems)[fNofTofdHitItems++]) R3BTofdHitData(tdet_sum[i],
                                                                                  xdet_sum[i] * 100.,
                                                                                  ydet_sum[i] * 100.,
                                                                                  qdet_sum[i],
                                                                                  -5,
                                                                                  qdet_sum[i],
                                                                                  ipl,
                                                                                  1,
                                                                                  0);

                        fh_check_QvsX[6]->Fill(xdet_sum[i] * 100, qdet_sum[i]);
                        fh_check_TvsX[6]->Fill(xdet_sum[i] * 100, tdet_sum[i]);
                        fh_check_XvsY[6]->Fill(xdet_sum[i] * 100, ydet_sum[i] * 100);
                    }
                }
                else // first and second plane will be written as one detector; this makes Dima's tracker
                     // faster
                {
                    if (((i == 6 || i == 8) && abs(fB) > 1710.) || ((i == 7 || i == 9) && abs(fB) < 1710.))
                    {
                        ax += xdet_sum[i];
                        ay += ydet_sum[i];
                        at += tdet_sum[i];
                        aq += qdet_sum[i];
                    }
                }
            }
        }

        if (ftrackerType == 1)
        {
            ax = ax / 2.;
            ay = ay / 2.;
            at = at / 2.;
            aq = aq / 2.;

            new ((*fTofdHitItems)[fNofTofdHitItems++]) R3BTofdHitData(at, ax * 100., ay * 100., aq, -5., aq, 10, 1, 0);

            // Fill check spectra:
            fh_check_QvsX[6]->Fill(ax * 100, aq);
            fh_check_TvsX[6]->Fill(ax * 100, at);
            fh_check_XvsY[6]->Fill(ax * 100, ay * 100);
        }
    }

    //  finish_from_cpp_();
    cout << "Statistics:" << endl;
    cout << "Total number of Events:         " << fNEvents << endl;
    cout << "Wrong Trigger:                  " << counterWrongTrigger << endl;
    cout << "Wrong Tpat:                     " << counterWrongTpat << endl;
    cout << "ROLU veto:                      " << counterRolu << endl;
    cout << "Califa true:                    " << counterCalifa << endl;
    cout << "Evets with nHits Tofd > 0:      " << counterTofd << endl;
    cout << "TofD average multi:             " << (float)counterTofd / (float)counterTofdMulti << endl;
    cout << "Number of single hits in Tofd:  " << Nevent_singles << endl;
    cout << "Number of hits with wrong Q:    " << Nevent_wrongcharge << endl;
    cout << "Events before fibs analysis:    " << hits1 << endl;
    cout << "Hits before fibcorel:           " << Nhits_before_fibcor << endl;
    cout << "Hits after fibcorel:            " << Nhits_after_fibcor << endl;
    cout << "Events after fibcorel:          " << Neventafterfibcor << endl;
    //    cout << "Events before goodQ:            " << fNeventpair <<endl;
    cout << "Events with goodQtofd:          " << Nevent_goodQtofd << endl;
    cout << "Events with goodQ:              " << Nevent_goodQ << endl;
    cout << "ncounts[0] < 2 for goodQ:       " << Nevent_zerofi23a_corel << endl;
    cout << "ncounts[0] > 1 for goodQ:       " << Nevent_notzerofi23a_corel << endl;
    cout << "ncounts[1] < 2 for goodQ:       " << Nevent_zerofi23b_corel << endl;
    cout << "ncounts[1] > 1 for goodQ:       " << Nevent_notzerofi23b_corel << endl;
    cout << "fi30/32 wrong counts for goodQ: " << Nevent_zerofi3032_corel << endl;
    cout << "fi30/32 good counts for goodQ:  " << Nevent_notzerofi3032_corel << endl;
    cout << "fi31/33 wrong counts for goodQ: " << Nevent_zerofi3133_corel << endl;
    cout << "fi31/33 good counts for goodQ:  " << Nevent_notzerofi3133_corel << endl;
    cout << "Final number of events:         " << fNeventselect << endl;
    cout << "Correlated pairs in x:          " << ncorellx << endl;
    cout << "Anticorrelated pairs in x:      " << nanticorellx << endl;
    cout << "Correlated pairs in y:          " << ncorelly << endl;
    cout << "Anticorrelated pairs in y:      " << nanticorelly << endl;
    cout << "Hits with countddet>50:         " << countdet50 << endl;

    /*
        cout << "Eff. Fi30 min: " << hits30 << "  " << hits30 / hits1 << endl;
        cout << "Eff. Fi31 min: " << hits31 << "  " << hits31 / hits1 << endl;
        cout << "Eff. Fi32 min: " << hits32 << "  " << hits32 / hits1 << endl;
        cout << "Eff. Fi33 min: " << hits33 << "  " << hits33 / hits1 << endl;

        cout << "Eff. Fi30 max: " << hits30bc << "  " << hits30bc / hits1 << endl;
        cout << "Eff. Fi31 max: " << hits31bc << "  " << hits31bc / hits1 << endl;
        cout << "Eff. Fi32 max: " << hits32bc << "  " << hits32bc / hits1 << endl;
        cout << "Eff. Fi33 max: " << hits33bc << "  " << hits33bc / hits1 << endl;
    */
    cout << "Number of written hits for selected nEvents: " << counter1 << ", counterwo: " << counterwo << endl;
    cout << "Fib23a: " << countdet_written[0] << endl;
    cout << "Fib23b: " << countdet_written[1] << endl;
    cout << "Fib30:  " << countdet_written[2] << endl;
    cout << "Fib31:  " << countdet_written[3] << endl;
    cout << "Fib32:  " << countdet_written[4] << endl;
    cout << "Fib33:  " << countdet_written[5] << endl;
    cout << "Tofd6:  " << countdet_written[6] << endl;
    cout << "Tofd7:  " << countdet_written[7] << endl;
    cout << "Tofd8:  " << countdet_written[8] << endl;
    cout << "Tofd9:  " << countdet_written[9] << endl;

    if (fHisto)
    {
        fh_Tpat->Write();
        fh_Trigger->Write();
        if (fMappedItems.at(DET_BMON))
        {
            fh_TOFDOR->Write();
            fh_SEE->Write();
            fh_IC->Write();
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
            if (fHitItems.at(ifibcount + DET_FI_FIRST))
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
                fh_qtofd_vs_qFib[ifibcount]->Write();
                fh_qtofd_vs_qFib_ac[ifibcount]->Write();
                fh_Fib_ToF_vs_Qtofd[ifibcount]->Write();
            }
        }

        fh_TimePreviousEventFi23a->Write();
        fh_TimePreviousEventFi23b->Write();
        fh_tofd_vs_Fib31->Write();
        fh_ytofd_vs_yFib31->Write();
        fh_tofd_vs_Fib31_ac->Write();
        fh_tofd_vs_Fib30->Write();
        fh_ytofd_vs_yFib30->Write();
        fh_tofd_vs_Fib30_ac->Write();
        fh_Fib33_vs_Fib31->Write();
        fh_Fib33_vs_Fib31_dx->Write();
        fh_Fib33_vs_Fib31_dt->Write();
        fh_Fib31_vs_Fib23a->Write();
        fh_Fib31_vs_Fib23a_dx->Write();
        fh_Fib31_vs_Fib23a_dt->Write();
        fh_Fib32_vs_Fib30->Write();
        fh_Fib32_vs_Fib30_dx->Write();
        fh_Fib32_vs_Fib30_dt->Write();
        fh_tofdleft_vs_Fib23b->Write();
        fh_Fib30_vs_Fib23b_dx->Write();
        fh_Fib30_vs_Fib23b_dt->Write();
        fh_Fib30_vs_Fib23a->Write();
        fh_Fib30_vs_Fib23a_dx->Write();
        fh_Fib30_vs_Fib23a_dt->Write();
        fh_tofdright_vs_Fib23b->Write();
        fh_Fib31_vs_Fib23b_dx->Write();
        fh_Fib31_vs_Fib23b_dt->Write();
        fh_Fib33_vs_Fib31_ac->Write();
        fh_Fib33_vs_Fib31_dx_ac->Write();
        fh_Fib31_vs_Fib23a_ac->Write();
        fh_Fib31_vs_Fib23a_dx_ac->Write();
        fh_Fib32_vs_Fib30_ac->Write();
        fh_Fib32_vs_Fib30_dx_ac->Write();
        fh_tofdleft_vs_Fib23b_ac->Write();
        fh_Fib30_vs_Fib23b_dx_ac->Write();
        fh_Fib30_vs_Fib23a_ac->Write();
        fh_Fib30_vs_Fib23a_dx_ac->Write();
        fh_tofdright_vs_Fib23b_ac->Write();
        fh_Fib31_vs_Fib23b_dx_ac->Write();
        fh_Fib32_vs_Fib30_tot->Write();
        fh_Fib33_vs_Fib31_tot->Write();
        fh_Fib23a_vs_Fib3X_tot->Write();
        fh_Fib23b_vs_Fib3X_tot->Write();
        fh_Fi23bY_TofdY_He->Write();
        fh_Fi23bY_TofdY_C->Write();
        fh_Fi23aX_TofdX_He->Write();
        fh_Fi23aX_TofdX_C->Write();
        fh_Fi301Y_TofdY_He->Write();
        fh_Fi301Y_TofdY_C->Write();
        fh_Fi301X_TofdX_He->Write();
        fh_Fi301X_TofdX_C->Write();
        fh_Fi323Y_TofdY_He->Write();
        fh_Fi323Y_TofdY_C->Write();
        fh_Fi323X_TofdX_He->Write();
        fh_Fi323X_TofdX_C->Write();
        fh_Fi23aX_TofdX_bw1->Write();
        fh_Fi23aX_TofdX_bw2->Write();
        fh_yTofd_vs_yTofd->Write();
        fh_yTofdHe_vs_yTofdC->Write();
        fh_xTofdHe_vs_xTofdC->Write();
        fh_xy_target->Write();
        fh_qsum_mult_fi23a->Write();
        fh_qsum_mult_fi23b->Write();
        fh_qtofi_vs_qtofd->Write();
        /*
           for(Int_t i=0; i<21;i++){
               fh_Fi23aToF_Q[i]->Write();
               fh_Fi23bToF_Q[i]->Write();
           }
           for (int i = 0; i < 6; i++)
           {
               fh_tofd_x_vs_y_z[i]->Write();
           }
       */
        for (int i = 0; i < 7; i++)
        {
            fh_check_QvsX[i]->Write();
            fh_check_TvsX[i]->Write();
            fh_check_XvsY[i]->Write();
        }
        if (fHitItems.at(DET_CALIFA))
        {
            fh_califa_energy->Write();
            fh_califa_energy_dc->Write();
            fh_califa_energy_select_barrel->Write();
            fh_califa_energy_select_iphos->Write();
            fh_califa_time->Write();
            fh_califa_tofd_barrel->Write();
            fh_califa_tofd_iphos->Write();
            fh_califa_cluster->Write();
            fh_califa_overflow_rolu->Write();
            fh_califa_overflow_antirolu->Write();
            fh_califa_cluster_rolu->Write();
            fh_califa_overflow->Write();
            fh_califa_good->Write();
            fh_califa_good_rolu->Write();
            fh_califa_good_antirolu->Write();
            fh_califa_motherid_clustertype->Write();
        }
    }
}

ClassImp(R3BPreTrackS494)
