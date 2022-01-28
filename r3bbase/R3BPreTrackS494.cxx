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
    , fNEvents(0)
    , fTofdHitItems(new TClonesArray("R3BTofdHitData"))
    , fFi23aHitItems(new TClonesArray("R3BFiberMAPMTHitData"))
    , fFi23bHitItems(new TClonesArray("R3BFiberMAPMTHitData"))
    , fFi30HitItems(new TClonesArray("R3BFiberMAPMTHitData"))
    , fFi31HitItems(new TClonesArray("R3BFiberMAPMTHitData"))
    , fFi32HitItems(new TClonesArray("R3BFiberMAPMTHitData"))
    , fFi33HitItems(new TClonesArray("R3BFiberMAPMTHitData"))
    , fNofTofdHitItems(0)
    , fNofFi23aHitItems(0)
    , fNofFi23bHitItems(0)
    , fNofFi30HitItems(0)
    , fNofFi31HitItems(0)
    , fNofFi32HitItems(0)
    , fNofFi33HitItems(0)
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

    LOG(INFO) << "R3BPreTrackS494::Init ";

    // try to get a handle on the EventHeader. EventHeader may not be
    // present though and hence may be null. Take care when using.

    FairRootManager* mgr = FairRootManager::Instance();
    if (NULL == mgr)
        LOG(fatal) << "FairRootManager not found";

    header = (R3BEventHeader*)mgr->GetObject("R3BEventHeader");
    FairRunOnline* run = FairRunOnline::Instance();

    // Get objects for detectors on all levels
    fMCTrack = (TClonesArray*)mgr->GetObject("MCTrack");
    if (fMCTrack)
        mgr->Register("MCTrack", "Monte Carlo Tracks", fMCTrack, kTRUE);

    assert(DET_MAX + 1 == sizeof(fDetectorNames) / sizeof(fDetectorNames[0]));
    printf("Have %d fiber detectors.\n", NOF_FIB_DET);
    fMappedItems.push_back((TClonesArray*)mgr->GetObject(Form("%sMappedData", fDetectorNames[0])));
    if (NULL == fMappedItems.at(0))
    {
        printf("Could not find mapped data for '%s'.\n", fDetectorNames[0]);
    }
    fCalItems.push_back((TClonesArray*)mgr->GetObject(Form("%sCrystalCalData", fDetectorNames[0])));
    if (NULL == fCalItems.at(0))
    {
        printf("Could not find Cal data for '%s'.\n", fDetectorNames[0]);
    }
    fHitItems.push_back((TClonesArray*)mgr->GetObject(Form("%sHit", fDetectorNames[0])));
    if (NULL == fHitItems.at(0))
    {
        printf("Could not find hit data for '%s'.\n", fDetectorNames[0]);
    }
    for (int det = 1; det < DET_MAX; det++)
    {
        fMappedItems.push_back((TClonesArray*)mgr->GetObject(Form("%sMapped", fDetectorNames[det])));
        if (NULL == fMappedItems.at(det))
        {
            printf("Could not find mapped data for '%s'.\n", fDetectorNames[det]);
        }
        if (det == 9)
            maxevent = mgr->CheckMaxEventNo();
        fCalItems.push_back((TClonesArray*)mgr->GetObject(Form("%sCal", fDetectorNames[det])));
        if (NULL == fCalItems.at(det))
        {
            printf("Could not find Cal data for '%s'.\n", fDetectorNames[det]);
        }
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

    //------------------------------------------------------------------------
    // graphical cuts
    //------------------------------------------------------------------------
    if (fGraphCuts)
    {
        cut_fi31_fi23a = NULL;
        cut_fi30_fi23b = NULL;
        cut_Fi33vsFi31 = NULL;
        cut_Fi30vsFi32 = NULL;

        TFile* f = TFile::Open("fiberCuts.root", "read");
        cut_Fi33vsFi31 = dynamic_cast<TCutG*>(f->Get("cut_dx_Fi33_Fi31"));
        cut_Fi30vsFi32 = dynamic_cast<TCutG*>(f->Get("cut_dx_Fi30_Fi32"));

        TFile* f23a = TFile::Open("myCutsFi23a.root", "read");
        cut_fi31_fi23a = dynamic_cast<TCutG*>(f23a->Get("cut_fi31_fi23a"));

        TFile* f23b = TFile::Open("myCutsFi23b.root", "read");
        cut_fi30_fi23b = dynamic_cast<TCutG*>(f23b->Get("cut_fi32_fi23b"));
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
    LOG(DEBUG) << fsens_SEE << ", " << fexp << ", " << fpow << ", " << calib_SEE << endl;

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

            LOG(DEBUG) << "I am creating canvas " << detName << endl;

            // xy:
            fh_xy_Fib[ifibcount] =
                new TH2F(Form("%s_xy", detName), Form("%s xy", detName), 1000, -50., 50., 1000, -50., 50);
            fh_xy_Fib[ifibcount]->GetXaxis()->SetTitle("x / cm ");
            fh_xy_Fib[ifibcount]->GetYaxis()->SetTitle("y / cm");

            fh_xy_Fib_ac[ifibcount] =
                new TH2F(Form("%s_xy_ac", detName), Form("%s xy after cuts", detName), 1000, -50., 50, 1000, -50., 50);
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
                                                1000,
                                                -50.,
                                                50,
                                                400,
                                                0.,
                                                400.);
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

            // hit fiber vs. fiber position:

        } // end if(Mapped)

    } // end for(ifibcount)
    fh_Fib33_vs_Fib31 = new TH2F("fib33_vs_fib31", "Fiber 33 vs. Fiber 31", 1000, -50, 50, 1000, -50., 50.);
    fh_Fib33_vs_Fib31->GetYaxis()->SetTitle("Fiber33");
    fh_Fib33_vs_Fib31->GetXaxis()->SetTitle("Fiber31");

    fh_Fib31_vs_Fib23a = new TH2F("fib31_vs_fib23a", "Fiber 31 vs. Fiber 23a", 1000, -50, 50, 1000, -50., 50.);
    fh_Fib31_vs_Fib23a->GetYaxis()->SetTitle("Fiber31");
    fh_Fib31_vs_Fib23a->GetXaxis()->SetTitle("Fiber23a");

    fh_Fib32_vs_Fib30 = new TH2F("fib32_vs_fib30", "Fiber 32 vs. Fiber 30", 1000, -50, 50, 1000, -50., 50.);
    fh_Fib32_vs_Fib30->GetYaxis()->SetTitle("Fiber32");
    fh_Fib32_vs_Fib30->GetXaxis()->SetTitle("Fiber30");

    fh_Fib30_vs_Fib23b = new TH2F("fib30_vs_fib23b", "Fiber 30 vs. Fiber 23b", 1000, -50, 50, 1000, -50., 50.);
    fh_Fib30_vs_Fib23b->GetYaxis()->SetTitle("Fiber30");
    fh_Fib30_vs_Fib23b->GetXaxis()->SetTitle("Fiber23b");

    fh_Fib30_vs_Fib23a = new TH2F("fib30_vs_fib23a", "Fiber 30 vs. Fiber 23a", 1000, -50, 50, 1000, -50., 50.);
    fh_Fib30_vs_Fib23a->GetYaxis()->SetTitle("Fiber30");
    fh_Fib30_vs_Fib23a->GetXaxis()->SetTitle("Fiber23a");

    fh_Fib31_vs_Fib23b = new TH2F("fib31_vs_fib23b", "Fiber 31 vs. Fiber 23b", 1000, -50, 50, 1000, -50., 50.);
    fh_Fib31_vs_Fib23b->GetYaxis()->SetTitle("Fiber31");
    fh_Fib31_vs_Fib23b->GetXaxis()->SetTitle("Fiber23b");

    // dx between fibers vs x
    fh_Fib33_vs_Fib31_dx = new TH2F("fib33_fib31_dx", "dx of Fiber 33 and Fiber 31", 1000, -50, 50, 1000, -50., 50.);
    fh_Fib33_vs_Fib31_dx->GetYaxis()->SetTitle("xFi33 - xFi31 / cm");
    fh_Fib33_vs_Fib31_dx->GetXaxis()->SetTitle("x Fi31 / cm");

    fh_Fib31_vs_Fib23a_dx = new TH2F("fib31_fib23a_dx", "dx of Fiber 31 and Fiber 23a", 1000, -50, 50, 1000, -50., 50.);
    fh_Fib31_vs_Fib23a_dx->GetYaxis()->SetTitle("xFi31 - xFi23a / cm");
    fh_Fib31_vs_Fib23a_dx->GetXaxis()->SetTitle("x Fi23a / cm");

    fh_Fib32_vs_Fib30_dx = new TH2F("fib32_fib30_dx", "dx of Fiber 32 and Fiber 30", 1000, -50, 50, 1000, -50., 50.);
    fh_Fib32_vs_Fib30_dx->GetYaxis()->SetTitle("xFi32 - xFi30 / cm");
    fh_Fib32_vs_Fib30_dx->GetXaxis()->SetTitle("x Fi30 / cm");

    fh_Fib30_vs_Fib23b_dx = new TH2F("fib32_fib23b_dx", "dx of Fiber 32 and Fiber 23b", 1000, -50, 50, 1000, -50., 50.);
    fh_Fib30_vs_Fib23b_dx->GetYaxis()->SetTitle("xFi32 - xFi23b / cm");
    fh_Fib30_vs_Fib23b_dx->GetXaxis()->SetTitle("x Fi23b / cm");

    fh_Fib31_vs_Fib23b_dx = new TH2F("fib31_fib23b_dx", "dx of Fiber 31 and Fiber 23b", 1000, -50, 50, 1000, -50., 50.);
    fh_Fib31_vs_Fib23b_dx->GetYaxis()->SetTitle("xFi31 - xFi23b / cm");
    fh_Fib31_vs_Fib23b_dx->GetXaxis()->SetTitle("x Fi23b / cm");

    fh_Fib30_vs_Fib23a_dx = new TH2F("fib30_fib23a_dx", "dx of Fiber 30 and Fiber 23a", 1000, -50, 50, 1000, -50., 50.);
    fh_Fib30_vs_Fib23a_dx->GetYaxis()->SetTitle("xFi30 - xFi23a / cm");
    fh_Fib30_vs_Fib23a_dx->GetXaxis()->SetTitle("x Fi23a / cm");

    //---------------------------------------------------------------------------------------------------
    // TofD detector

    if (fHitItems.at(DET_TOFD))
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
    for (int i = 0; i < 6; i++)
    {
        fh_check_QvsX[i] = new TH2F(Form("fhCheckQvsX%d", i), Form("Check Q vs X det%d", i), 1000, -50, 50, 20, 0, 10);
        fh_check_TvsX[i] =
            new TH2F(Form("fhCheckTvsX%d", i), Form("Check T vs X det%d", i), 1000, -50, 50, 500, -100, 100);
        fh_check_XvsY[i] =
            new TH2F(Form("fhCheckXvsY%d", i), Form("Check X vs Y det%d", i), 600, -30, 30, 600, -30, 30);
    }
    fh_check_QvsX[6] = new TH2F(Form("fhCheckQvsX%d", 6), Form("Check Q vs X det%d", 6), 200, -100., 100., 20, 0, 10);
    fh_check_TvsX[6] =
        new TH2F(Form("fhCheckTvsX%d", 6), Form("Check T vs X det%d", 6), 200, -100., 100., 500, -100, 100);
    fh_check_XvsY[6] =
        new TH2F(Form("fhCheckXvsY%d", 6), Form("Check X vs Y det%d", 6), 100, -135., 135., 200, -100., 100.);

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
    if (fNEvents / 10000. == (int)fNEvents / 10000)
        std::cout << "\rEvents: " << fNEvents << " / " << maxevent << " (" << (int)(fNEvents * 100. / maxevent)
                  << " %) " << std::flush;
    // cout << "Event: " << fNEvents << endl;
    fNEvents += 1;

    FairRootManager* mgr = FairRootManager::Instance();
    if (NULL == mgr)
        LOG(fatal) << "FairRootManager not found";

    if (header && !fSimu)
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

    if (fMappedItems.at(DET_BMON))
    {
        unsigned long IC;
        unsigned long SEETRAM_raw;
        Double_t SEETRAM;
        unsigned long TOFDOR;

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
        // cout << "ROLU cut applied !!!" << endl;
        counterRolu++;
        return;
    }

    Bool_t CalifaHit = false;
    if (fMappedItems.at(DET_CALIFA))
    {
        // CALIFA
        auto detCalifa = fMappedItems.at(DET_CALIFA);
        Int_t nHitsCalifa = detCalifa->GetEntriesFast();
        // cout<<"Califa hits: "<<nHitsCalifa<<endl;

        for (Int_t ihit = 0; ihit < nHitsCalifa; ihit++)
        {
            R3BCalifaMappedData* hitCalifa = (R3BCalifaMappedData*)detCalifa->At(ihit);
            if (!hitCalifa)
                continue;

            Int_t Crystal = hitCalifa->GetCrystalId();
            Int_t Energy = hitCalifa->GetEnergy();
            // cout << "Califa: " << Crystal << " Energy: " << Energy << endl;
            if (Energy > 0)
            {
                fh_califa_energy->Fill(Crystal, Energy);
                CalifaHit = true;
            }
        }
    }
    if (CalifaHit)
    {
        counterCalifa++;
        //		return;
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
            LOG(DEBUG) << "PID " << PID << endl;
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

                    LOG(DEBUG) << "******************************************" << endl;
                    LOG(DEBUG) << "Track In 4He"
                               << "x " << XHes << " y " << YHes << " z " << ZHes << endl;
                    LOG(DEBUG) << "px " << pHexs << " py " << pHeys << " z " << pHezs << endl;
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

                    LOG(DEBUG) << "******************************************" << endl;
                    LOG(DEBUG) << "Track In 12C"
                               << "x " << XCs << " y " << YCs << " z " << ZCs << endl;
                    LOG(DEBUG) << "px " << pCxs << " py " << pCys << " z " << pCzs << endl;
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

                    LOG(DEBUG) << "******************************************" << endl;
                    LOG(DEBUG) << "Track In 16O"
                               << "x " << Xf << " y " << Yf << " z " << Zf << endl;
                    LOG(DEBUG) << "px " << Pxf << " py " << Pyf << " z " << Pzf << endl;
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
    Double_t zdet_s[max];
    Double_t tdet[max];
    Int_t qdet[max];
    Double_t qdet_s[max] = { 0 };
    Double_t xdet_s[max] = { 0 };
    Double_t ydet_s[max] = { 0 };
    Double_t tdet_s[max] = { 0 };
    Int_t detector_s[max] = { 0 };
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
    for (int i = 0; i < max; i++)
    {
        fFi33[i] = false;
        fFi32[i] = false;
        fFi31[i] = false;
        fFi30[i] = false;
        fFi23a[i] = false;
        fFi23b[i] = false;
    }

    countdet = 0;

    Int_t n_det = 10;

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
    Bool_t debug = false;
    Bool_t debug2 = false;
    Bool_t debug_in = false;
    Bool_t debug3 = false;
    Bool_t twice = false; // true: one part after another tracked; false: both part traacked simult.
    Bool_t fibCuts = true;
    // only consider fiber with maximum energy deposit, only for sweep runs with beam
    Bool_t maxWerte = false;
    // if (!fPairs && fB != -1710 && !fSimu)
    //   maxWerte = true;
    // if fibCuts true, dx1 (Fib3xvs3x), dx2(Fib3xvs2x), dx3(Fib3xvsTofd) used for cuts in xfib correlations
    Double_t dx1 = 2;   // cm
    Double_t dx2 = 5.;  // cm
    Double_t dx3 = 4.;  // cm
    Double_t dtft = 8.; // ns
                        // END CHOOSING OPTIONS **********

    for (int i = 0; i < n_det; i++)
    {
        x[i] = -1000.;
        y[i] = -1000.;
        z[i] = 0.;
        q[i] = 0.;
        t[i] = -1000.;

        x1[i] = -1000.;
        y1[i] = -1000.;
        z1[i] = 0.;
        q1[i] = 0.;
        t1[i] = -1000.;

        x2[i] = -1000.;
        y2[i] = -1000.;
        z2[i] = 0.;
        q2[i] = 0.;
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
    Int_t tofd1r = 6;
    Int_t tofd1l = 7;
    Int_t tofd2r = 8;
    Int_t tofd2l = 9;

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

    if (!fPairs)
    {
        auto detHit23a = fHitItems.at(DET_FI23A);
        Int_t nHits23a = detHit23a->GetEntriesFast();
        if (nHits23a > 0)
            return;
        auto detHit23b = fHitItems.at(DET_FI23B);
        Int_t nHits23b = detHit23b->GetEntriesFast();
        if (nHits23b > 0)
            return;
    }

    //   cout<<"*** Entering analysis ***"<<endl;
    for (Int_t i = 0; i < 10; i++)
    {
        tPrev[i] = -1000.;
        detPrev[i] = -1;
    }

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

    if (nHits > 100)
        return;

    Int_t multTofd = 0;

    // loop over ToFD
    for (Int_t ihit = 0; ihit < nHits; ihit++)
    {

        R3BTofdHitData* hitTofd = (R3BTofdHitData*)detTofd->At(ihit);
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

        if (fCuts && (ttt < -100. || ttt > 100.) && !fSimu) // change time cuts
        {                                                   // trigger window -1500, 1500
            if (debug_in)
                cout << "No trigger particle!" << endl;
            continue;
        }

        Double_t qqq = hitTofd->GetEloss(); // / 1.132;
        // if (fB == -1710)
        //{
        //    qqq = qqq * 1.11; // / 1.132;
        //}

        if (fSimu)
            qqq = qqq * 26.76 * 3.2706460;

        fh_tofd_charge->Fill(qqq);

        Double_t xxx = hitTofd->GetX();
        Double_t yyy = hitTofd->GetY();

        if (!fSimu && !fPairs)
        {
            yyy = yyy * (-1.); // -1 until we solve problem with y direction
            if (qqq > 1.4 && qqq < 2.6)
                yyy = yyy * 0.5;
        }

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

            if (!fPairs) // fPairs given in macro, true 2 particles, false 16O tracking
            {
                if (qqq < 7.2 || qqq > 8.8)
                {
                    if (debug_in)
                        cout << "Not the right charge! Charge <> 8" << endl;

                    continue;
                }
                else
                {
                    qqq = 8.; // this istemp, as hitpar from He run are not quite good any more.
                }
            }
            if (fPairs && !(qqq > 1.5 && qqq < 2.5) && !(qqq > 5.5 && qqq < 6.5))
            {
                if (debug_in)
                    cout << "Not the right charge! Charge = " << qqq << " Charge <> 2 and Charge <> 6" << endl;
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
        x2[det2] = xxx / 100.; // for tracker everything in meters
        y2[det2] = yyy / 100. + y_corr;
        if (y2[det2] < -0.8 || y2[det2] > 0.8) // in meters!
        {
            if (debug_in)
                cout << "Wrong ytoFD-position! " << y2[det2] << " for charge " << q2[det2] << endl;
            continue;
        }

        z2[det2] = 0.;
        q2[det2] = qqq;
        t2[det2] = hitTofd->GetTime();

        fh_xy_tofd_ac->Fill(x2[det2] * 100. + 0. * randx, y2[det2] * 100.);
        fh_tofd_charge_ac->Fill(q2[det2]);
        fh_tofd_time_ac->Fill(t2[det2]);

        // register hits for tracker as long a time is in the coincidence window
        if ((abs(t2[det2] - t1[det1]) < 50.) || first) // change back to 5.
        {
            if (debug_in)
                cout << "2 particle within 5 ns   " << first << endl;
            // register point for tracker
            detector[countdet] = det2;
            xdet[countdet] = x2[det2];
            ydet[countdet] = y2[det2];
            zdet[countdet] = z2[det2];
            tdet[countdet] = t2[det2];
            qdet[countdet] = (int)(q2[det2] + 0.5); // q for tracker must be integer

            if (debug_in)
            {
                cout << "registered"
                     << " x: " << x2[det2] << " y: " << y2[det2] << " q: " << q2[det2] << " t: " << t2[det2] << " ID "
                     << id2 << endl;
            }
            if (abs(qdet[countdet] - 2.) < 0.5)
            {
                alpha = true;
                x_4He = xdet[countdet];
                y_4He = ydet[countdet];
                z_4He = zdet[countdet];
            }
            if (abs(qdet[countdet] - 6.) < 0.5)
            {
                carbon = true;
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
        {
            cout << "found good event, x_12C " << x_12C * 100. << " y_12C " << y_12C * 100. << endl;
            cout << "found good event, x_4He " << x_4He * 100. << " y_4He " << y_4He * 100. << endl;
            cout << "Delta: " << delta << endl;
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

        multTofd++;

        counterTofdMulti++;

        // cut in ToT for Fibers
        Double_t cutQ = 0.;
        if (!fPairs || fB != -1710)
            // check if cut can be applied
            cutQ = -10.;

        if (debug_in)
            cout << "start fiber analysis" << endl;

        // loop over fiber 33
        auto detHit33 = fHitItems.at(DET_FI33);
        Int_t nHits33 = detHit33->GetEntriesFast();
        LOG(DEBUG) << "Fi33 hits: " << nHits33 << endl;

        Int_t mult33 = 0;
        for (Int_t ihit33 = 0; ihit33 < nHits33; ihit33++)
        {
            det = fi33;
            R3BBunchedFiberHitData* hit33 = (R3BBunchedFiberHitData*)detHit33->At(ihit33);
            x1[det] = hit33->GetX() / 100.;
            y1[det] = hit33->GetY() / 100.;
            z1[det] = 0.;
            q1[det] = hit33->GetEloss();

            /*     if (fSimu)
                 {
                     if (q1[det] > 7.) // these cuts are for simulations
                         q1[det] = 8.;
                     else if (q1[det] > 5.)
                         q1[det] = 6.;
                     else if (q1[det] > 0.)
                         q1[det] = 2.;
                 }*/

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

            if (debug3)
                cout << "Fi33 bc: " << ihit33 << " x1: " << x1[det] << " y1: " << y1[det] << " q1: " << q1[det]
                     << " t1: " << t1[det] << endl;

            hits33bc++;

            // Cuts on Fi33
            if (fCuts && (t1[det] < -30 || t1[det] > 30) && !fSimu && fB == -1710)
                continue;
            if (fCuts && (t1[det] < -40 || t1[det] > 40) && !fSimu && fB != -1710)
                continue;
            if (fCuts && (x1[det] < -0.3 || x1[det] > 0.3))
                continue;
            if (fCuts && (y1[det] < -0.3 || y1[det] > 0.3))
                continue;
            if (fCuts && !fPairs && q1[det] < cutQ)
                continue;
            if (fCuts && (tof < -50 || tof > 50) && !fSimu && fB == -1710)
                continue;
            if (fCuts && (tof < 20 || tof > 90) && !fSimu && fB != -1710)
                continue;
            //            if (fGraphCuts && !cut_Fi33vsTofd->IsInside(x1[tofd1r] * 100., x1[det] * 100.))
            //                continue;

            hits33++;

            xFi33[mult33] = x1[det];
            yFi33[mult33] = y1[det];
            if (q1[det] > 7. && (q1[tofd1r] > 5.5 || q1[tofd2r] > 5.5))
                qFi33[mult33] = 6;
            else
                qFi33[mult33] = 0;
            tFi33[mult33] = tof;
            timeFi33[mult33] = t1[det];

            mult33++;
            if (mult33 > 100)
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
                countdet++;
            }

            if (debug3)
                cout << "Fi33: " << ihit33 << " x1: " << x1[det] << " y1: " << y1[det] << " q1: " << q1[det]
                     << " t1: " << t1[det] << endl;
        }

        if (mult33 > 0 && maxWerte)
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

            countdet++;
        }

        // loop over fiber 31
        auto detHit31 = fHitItems.at(DET_FI31);
        Int_t nHits31 = detHit31->GetEntriesFast();
        LOG(DEBUG) << "Fi31 hits: " << nHits31 << endl;
        Int_t mult31 = 0;
        for (Int_t ihit31 = 0; ihit31 < nHits31; ihit31++)
        {
            det = fi31;
            R3BBunchedFiberHitData* hit31 = (R3BBunchedFiberHitData*)detHit31->At(ihit31);
            x1[det] = hit31->GetX() / 100.;
            y1[det] = hit31->GetY() / 100.;
            z1[det] = 0.;
            q1[det] = hit31->GetEloss();

            /*  if (fSimu)
              {
                  if (q1[det] > 7.)
                      q1[det] = 8.;
                  else if (q1[det] > 5.)
                      q1[det] = 6.;
                  else if (q1[det] > 0.)
                      q1[det] = 2.;
              }*/

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

            if (debug3)
                cout << "Fi31 bc: " << ihit31 << " x1: " << x1[det] << " y1: " << y1[det] << " q1: " << q1[det]
                     << " t1: " << t1[det] << endl;
            hits31bc++;

            // Cuts on Fi31
            //            if (fCuts && x1[det] * 100. < -24.4)
            // if (fCuts && x1[det] * 100. < -25.75)
            // continue;
            if (fCuts && (t1[det] < -30 || t1[det] > 30) && !fSimu && fB == -1710)
                continue;
            if (fCuts && (t1[det] < -40 || t1[det] > 40) && !fSimu && fB != -1710)
                continue;
            if (fCuts && (x1[det] < -0.3 || x1[det] > 0.3))
                continue;
            if (fCuts && (y1[det] < -0.3 || y1[det] > 0.3))
                continue;
            if (fCuts && !fPairs && q1[det] < cutQ)
                continue;
            if (fCuts && (tof < -50 || tof > 50) && !fSimu && fB == -1710)
                continue;
            if (fCuts && (tof < 20 || tof > 90) && !fSimu && fB != -1710)
                continue;
            //            if (fGraphCuts && !cut_Fi33vsTofd->IsInside(x1[tofd1r] * 100., x1[det] * 100.))
            //                continue;

            hits31++;

            xFi31[mult31] = x1[det];
            yFi31[mult31] = y1[det];
            if (q1[det] > 6.7 && (q1[tofd1r] > 5.5 || q1[tofd2r] > 5.5))
                qFi31[mult31] = 6.;
            else
                qFi31[mult31] = 0;
            tFi31[mult31] = tof;
            timeFi31[mult31] = t1[det];
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
                countdet++;
            }

            if (debug3)
                cout << "Fi31: " << ihit31 << " x1: " << x1[det] << " y1: " << y1[det] << " q1: " << q1[det]
                     << " t1: " << t1[det] << endl;
        }

        if (mult31 > 0 && maxWerte)
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

            countdet++;
        }

        // loop over fiber 32
        auto detHit32 = fHitItems.at(DET_FI32);
        Int_t nHits32 = detHit32->GetEntriesFast();
        LOG(DEBUG) << "Fi32 hits: " << nHits32 << endl;
        Int_t mult32 = 0;
        for (Int_t ihit32 = 0; ihit32 < nHits32; ihit32++)
        {
            det = fi32;
            R3BBunchedFiberHitData* hit32 = (R3BBunchedFiberHitData*)detHit32->At(ihit32);
            x1[det] = hit32->GetX() / 100.;
            y1[det] = hit32->GetY() / 100.;
            z1[det] = 0.;
            q1[det] = hit32->GetEloss();

            /*  if (fSimu)
              {
                  if (q1[det] > 7.)
                      q1[det] = 8.;
                  else if (q1[det] > 5.)
                      q1[det] = 6.;
                  else if (q1[det] > 0.)
                      q1[det] = 2.;
              }*/

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

            // Cuts on Fi32
            if (fCuts && (t1[det] < -30 || t1[det] > 30) && !fSimu && fB == -1710)
                continue;
            if (fCuts && (t1[det] < -40 || t1[det] > 40) && !fSimu && fB != -1710)
                continue;
            if (fCuts && (x1[det] < -0.3 || x1[det] > 0.3))
                continue;
            if (fCuts && !fPairs && q1[det] < cutQ)
                continue;
            if (fCuts && (tof < -50 || tof > 50) && !fSimu && fB == -1710)
                continue;
            if (fCuts && (tof < 20 || tof > 90) && !fSimu && fB != -1710)
                continue;

            hits32++;

            xFi32[mult32] = x1[det];
            yFi32[mult32] = y1[det];
            if (q1[det] > 7 && (q1[tofd1l] > 5.5 || q1[tofd2l] > 5.5))
                qFi32[mult32] = 6.;
            else
                qFi32[mult32] = 0;
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
                countdet++;
            }

            if (debug3)
                cout << "Fi32: " << ihit32 << " x1: " << x1[det] << " y1: " << y1[det] << " q1: " << q1[det]
                     << " t1: " << t1[det] << endl;
        }

        if (mult32 > 0 && maxWerte)
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

            countdet++;
        }

        // loop over fiber 30
        auto detHit30 = fHitItems.at(DET_FI30);
        Int_t nHits30 = detHit30->GetEntriesFast();
        LOG(DEBUG) << "Fi30 hits: " << nHits30 << endl;
        Int_t mult30 = 0;
        for (Int_t ihit30 = 0; ihit30 < nHits30; ihit30++)
        {
            det = fi30;
            R3BBunchedFiberHitData* hit30 = (R3BBunchedFiberHitData*)detHit30->At(ihit30);
            x1[det] = hit30->GetX() / 100.;
            y1[det] = hit30->GetY() / 100.;
            z1[det] = 0.;
            q1[det] = hit30->GetEloss();

            /*   if (fSimu)
               {
                   if (q1[det] > 7.)
                       q1[det] = 8.;
                   else if (q1[det] > 5.)
                       q1[det] = 6.;
                   else if (q1[det] > 0.)
                       q1[det] = 2.;
               }*/

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

            // Cuts on Fi30
            if (fCuts && (t1[det] < -30 || t1[det] > 30) && !fSimu && fB == -1710)
                continue;
            if (fCuts && (t1[det] < -40 || t1[det] > 40) && !fSimu && fB != -1710)
                continue;
            if (fCuts && (x1[det] < -0.3 || x1[det] > 0.3))
                continue;
            if (fCuts && !fPairs && q1[det] < cutQ)
                continue;
            if (fCuts && (tof < -50 || tof > 50) && !fSimu && fB == -1710)
                continue;
            if (fCuts && (tof < 20 || tof > 90) && !fSimu && fB != -1710)
                continue;

            hits30++;

            xFi30[mult30] = x1[det];
            yFi30[mult30] = y1[det];
            if (q1[det] > 7.2 && (q1[tofd1l] > 5.5 || q1[tofd2l] > 5.5))
                qFi30[mult30] = 6.;
            else
                qFi30[mult30] = 0;
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
                countdet++;
            }

            if (debug3)
                cout << "Fi30: " << ihit30 << " x1: " << x1[det] << " y1: " << y1[det] << " q1: " << q1[det]
                     << " t1: " << t1[det] << endl;
        }

        if (mult30 > 0 && maxWerte)
        {
            detector[countdet] = fi30;
            xdet[countdet] = xMax[fi30];
            ydet[countdet] = yMax[fi30];
            zdet[countdet] = zMax[fi30];
            qdet[countdet] = qMax[fi30];
            if (fSimu)
                qdet[countdet] = qdet[countdet] * 8. / 2.5;
            tdet[countdet] = tMax[fi30];

            fh_Fib_ToF_ac[det]->Fill(xdet[countdet] * 100., tStart - tdet[countdet]);
            fh_xy_Fib_ac[det]->Fill(xdet[countdet] * 100., ydet[countdet] * 100.);
            fh_mult_Fib_ac[det]->Fill(mult33);
            fh_ToT_Fib_ac[det]->Fill(xdet[countdet] * 100., qdet[countdet]);
            fh_Fibs_vs_Tofd_ac[det]->Fill(x1[tofd1r] * 100. + randx, xdet[countdet] * 100.);
            fh_Fib_vs_Events_ac[det]->Fill(fNEvents, xdet[countdet] * 100.);
            fh_ToF_vs_Events_ac[det]->Fill(fNEvents, tStart - tdet[countdet]);
            fh_Fib_Time_ac[det]->Fill(xdet[countdet] * 100., tdet[countdet]);

            countdet++;
        }

        // loop over fiber 23a
        Int_t mult23a = 0;
        auto detHit23a = fHitItems.at(DET_FI23A);
        Int_t nHits23a = detHit23a->GetEntriesFast();
        LOG(DEBUG) << "Fi23a hits: " << nHits23a << endl;
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

            // Fill histograms before cuts
            fh_Fib_ToF[det]->Fill(x1[det] * 100., tof);
            fh_xy_Fib[det]->Fill(x1[det] * 100., y1[det] * 100.);
            fh_mult_Fib[det]->Fill(nHits23a);
            fh_ToT_Fib[det]->Fill(x1[det] * 100., q1[det]);
            fh_Fibs_vs_Tofd[det]->Fill(x1[tofd1r] * 100. + randx, x1[det] * 100.);
            fh_Fib_vs_Events[det]->Fill(fNEvents, x1[det] * 100.);
            fh_ToF_vs_Events[det]->Fill(fNEvents, tof);
            fh_Fib_Time[det]->Fill(x1[det] * 100., t1[det]);

            // Cuts on Fi23a
            if (fCuts && (t1[det] < -30 || t1[det] > 30) && !fSimu && fB == -1710)
                continue;
            if (fCuts && (t1[det] < -40 || t1[det] > 40) && !fSimu && fB != -1710)
                continue;
            if (fCuts && (x1[det] < -0.06 || x1[det] > 0.06))
                continue;
            if (fCuts && (y1[det] < -0.06 || y1[det] > 0.06))
                continue;
            if (fCuts && !fPairs && q1[det] < cutQ)
                continue;
            if (fCuts && (tof < -50 || tof > 50) && !fSimu && fB == -1710)
                continue;
            if (fCuts && (tof < 20 || tof > 90) && !fSimu && fB != -1710)
                continue;

            xFi23a[mult23a] = x1[det];
            yFi23a[mult23a] = y1[det];
            qFi23a[mult23a] = 0.; // q1[det];
            tFi23a[mult23a] = tof;
            timeFi23a[mult33] = t1[det];
            mult23a++;
            if (mult23a > 100)
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
                countdet++;
            }

            if (debug3)
                cout << "Fi23a " << ihit23a << " x1: " << x1[det] << " y1: " << y1[det] << " q1: " << q1[det]
                     << " t1: " << tof << endl;
        }

        if (mult23a > 0 && maxWerte)
        {
            detector[countdet] = fi23a;
            xdet[countdet] = xMax[fi23a];
            ydet[countdet] = yMax[fi23a];
            zdet[countdet] = zMax[fi23a];
            qdet[countdet] = 0;
            tdet[countdet] = tMax[fi23a];

            fh_Fib_ToF_ac[det]->Fill(xdet[countdet] * 100., tStart - tdet[countdet]);
            fh_xy_Fib_ac[det]->Fill(xdet[countdet] * 100., ydet[countdet] * 100.);
            fh_mult_Fib_ac[det]->Fill(mult33);
            fh_ToT_Fib_ac[det]->Fill(xdet[countdet] * 100., qdet[countdet]);
            fh_Fibs_vs_Tofd_ac[det]->Fill(x1[tofd1r] * 100. + randx, xdet[countdet] * 100.);
            fh_Fib_vs_Events_ac[det]->Fill(fNEvents, xdet[countdet] * 100.);
            fh_ToF_vs_Events_ac[det]->Fill(fNEvents, tStart - tdet[countdet]);
            fh_Fib_Time_ac[det]->Fill(xdet[countdet] * 100., tdet[countdet]);

            countdet++;
        }

        // loop over fiber 23b
        Int_t mult23b = 0;
        auto detHit23b = fHitItems.at(DET_FI23B);
        Int_t nHits23b = detHit23b->GetEntriesFast();
        LOG(DEBUG) << "Fi23b hits: " << nHits23b << endl;
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
            fh_Fibs_vs_Tofd[det]->Fill(y1[tofd1l] * 100., y1[det] * 100.);
            fh_Fib_vs_Events[det]->Fill(fNEvents, y1[det] * 100.);
            fh_ToF_vs_Events[det]->Fill(fNEvents, tof);
            fh_Fib_Time[det]->Fill(y1[det] * 100., t1[det]);

            // Cuts on Fi23b
            if (fCuts && (t1[det] < -30 || t1[det] > 30) && !fSimu && fB == -1710)
                continue;
            if (fCuts && (t1[det] < -40 || t1[det] > 40) && !fSimu && fB != -1710)
                continue;
            if (fCuts && (y1[det] < -0.06 || y1[det] > 0.06))
                continue;
            if (fCuts && (x1[det] < -0.06 || x1[det] > 0.06))
                continue;
            if (fCuts && !fPairs && q1[det] < cutQ)
                continue;
            if (fCuts && (tof < -50 || tof > 50) && !fSimu && fB == -1710)
                continue;
            if (fCuts && (tof < 20 || tof > 90) && !fSimu && fB != -1710)
                continue;

            xFi23b[mult23b] = x1[det];
            yFi23b[mult23b] = y1[det];
            qFi23b[mult23b] = 0.; // q1[det];
            tFi23b[mult23b] = tof;
            timeFi23b[mult33] = t1[det];
            mult23b++;

            if (mult23b > 100)
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
                countdet++;
            }

            if (debug3)
                cout << "Fi23b " << ihit23b << " x1: " << x1[det] << " y1: " << y1[det] << " q1: " << q1[det]
                     << " t1: " << tof << endl;
        }

        if (mult23b > 0 && maxWerte)
        {
            detector[countdet] = fi23b;
            xdet[countdet] = xMax[fi23b];
            ydet[countdet] = yMax[fi23b];
            zdet[countdet] = zMax[fi23b];
            qdet[countdet] = 0;
            tdet[countdet] = tMax[fi23b];

            fh_Fib_ToF_ac[det]->Fill(xdet[countdet] * 100., tStart - tdet[countdet]);
            fh_xy_Fib_ac[det]->Fill(xdet[countdet] * 100., ydet[countdet] * 100.);
            fh_mult_Fib_ac[det]->Fill(mult33);
            fh_ToT_Fib_ac[det]->Fill(xdet[countdet] * 100., qdet[countdet]);
            fh_Fibs_vs_Tofd_ac[det]->Fill(x1[tofd1r] * 100. + randx, xdet[countdet] * 100.);
            fh_Fib_vs_Events_ac[det]->Fill(fNEvents, xdet[countdet] * 100.);
            fh_ToF_vs_Events_ac[det]->Fill(fNEvents, tStart - tdet[countdet]);
            fh_Fib_Time_ac[det]->Fill(xdet[countdet] * 100., tdet[countdet]);

            countdet++;
        }

        // if we have magnetic field runs we do not have hits in Fib23
        if (mult23a == 0 && mult23b == 0 && !fPairs && fB != -1710)
        {
            det = fi23a;
            detector[countdet] = det;
            xdet[countdet] = 0.;
            ydet[countdet] = 0.;
            zdet[countdet] = 0.;
            qdet[countdet] = 8.;
            tdet[countdet] = 0.;
            // mult23a++;
            countdet++;

            det = fi23b;
            detector[countdet] = det;
            xdet[countdet] = 0.;
            ydet[countdet] = 0.;
            zdet[countdet] = 0.;
            qdet[countdet] = 8.;
            tdet[countdet] = 0.;
            // mult23b++;
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

        // FIBCORREL
        if (fibCuts)
        {
            // Plots of correlations of Fiber detectors and register events for tracker
            for (Int_t i = 0; i < mult31; i++)
            {
                if (debug3)
                    cout << "Fib31: " << i << " x: " << xFi31[i] << " q: " << qFi31[i] << endl;

                Double_t xtemp31 = -0.732 * x1[tofd1r] * 100. - 29.364;
                if (std::abs(xtemp31 - xFi31[i] * 100.) > dx3)
                    continue;

                fh_Fib_ToF_ac[fi31]->Fill(xFi31[i] * 100., tFi31[i]);
                fh_xy_Fib_ac[fi31]->Fill(xFi31[i] * 100., yFi31[i] * 100.);
                fh_ToT_Fib_ac[fi31]->Fill(xFi31[i] * 100., qFi31[i]);
                fh_Fibs_vs_Tofd_ac[fi31]->Fill(x1[tofd1r] * 100. + randx, xFi31[i] * 100.);
                fh_Fib_Time_ac[fi31]->Fill(xFi31[i] * 100., tStart - tFi31[i]);

                for (Int_t j = 0; j < mult33; j++)
                {
                    if (debug3)
                        cout << "Fib33: " << j << " x: " << xFi33[j] << " q: " << qFi33[j] << endl;
                    Double_t x31 = 0.866 * xFi33[j] * 100. - 4.637; // 1.1165 * xFi31[i] * 100. + 2.53856;
                    // if(xFi33[j] > -100 && xFi31[i] > -100)
                    tempFi33 = false;
                    if (fSimu)
                        tempFi33 = true;
                    if (!fSimu && abs(timeFi31[i] - timeFi33[j]) < dtft)
                        tempFi33 = true;
                    if (abs(xFi31[i] * 100. - x31) < dx1 && xFi31[i] * 100. > -30. && xFi33[j] * 100. > -30. &&
                        tempFi33)
                    {
                        fh_Fib33_vs_Fib31->Fill(xFi31[i] * 100., xFi33[j] * 100.);
                        fh_Fib33_vs_Fib31_dx->Fill(xFi31[i] * 100., xFi33[j] * 100. - xFi31[i] * 100.);

                        fh_Fib_ToF_ac[fi33]->Fill(xFi33[j] * 100., tFi33[j]);
                        fh_xy_Fib_ac[fi33]->Fill(xFi33[j] * 100., yFi33[j] * 100.);
                        fh_ToT_Fib_ac[fi33]->Fill(xFi33[j] * 100., qFi33[j]);
                        fh_Fibs_vs_Tofd_ac[fi33]->Fill(x1[tofd1r] * 100. + randx, xFi33[j] * 100.);
                        fh_Fib_Time_ac[fi33]->Fill(xFi33[j] * 100., tStart - tFi33[j]);

                        cond1 = kTRUE;
                        if (debug3)
                            cout << "cond1" << endl;

                        if (!maxWerte)
                        {
                            if (!fFi31[i])
                            {
                                detector[countdet] = fi31;
                                xdet[countdet] = xFi31[i];
                                ydet[countdet] = yFi31[i];
                                zdet[countdet] = 0.;
                                qdet[countdet] = qFi31[i];
                                tdet[countdet] = tFi31[i];
                                countdet++;
                                fFi31[i] = true;
                            }
                            if (!fFi33[j])
                            {
                                detector[countdet] = fi33;
                                xdet[countdet] = xFi33[j];
                                ydet[countdet] = yFi33[j];
                                zdet[countdet] = 0.;
                                qdet[countdet] = qFi33[j];
                                tdet[countdet] = tFi33[j];
                                countdet++;
                                fFi33[j] = true;
                            }
                        }
                    }
                }
                for (Int_t j = 0; j < mult23a; j++)
                {
                    if (debug3)
                        cout << "Fib23a: " << j << " x: " << xFi23a[j] << " q: " << qFi23a[j] << endl;
                    Double_t x31 = -6.625 * xFi23a[j] * 100. - 29.29; //-6.8672 * xFi23a[j] * 100. - 27.3507;
                    // if (abs(xFi31[i] * 100. - x31) < dx2 && xFi31[i] > -100 && xFi23a[j] > -100)
                    tempFi23ar = false;
                    if (fSimu)
                        tempFi23ar = true;
                    if (!fSimu && abs(timeFi31[i] - timeFi23a[j]) < dtft)
                        tempFi23ar = true;
                    // if (abs(xFi31[i] * 100. - x31) < dx2 && xFi31[i]*100. > -30 && xFi23a[j]*100. > -6 &&
                    if (xFi23a[j] * 100. < -0.05 && tempFi23ar)
                    {
                        // if (fGraphCuts && !cut_fi31_fi23a->IsInside(xFi23a[j] * 100., xFi31[i] * 100.))
                        //	continue;
                        fh_Fib31_vs_Fib23a->Fill(xFi23a[j] * 100., xFi31[i] * 100.);
                        fh_Fib31_vs_Fib23a_dx->Fill(xFi23a[j] * 100., xFi31[i] * 100. - xFi23a[j] * 100.);

                        fh_Fib_ToF_ac[fi23a]->Fill(xFi23a[j] * 100., tFi23a[j]);
                        fh_xy_Fib_ac[fi23a]->Fill(xFi23a[j] * 100., yFi23a[j] * 100.);
                        fh_ToT_Fib_ac[fi23a]->Fill(xFi23a[j] * 100., qFi23a[j]);
                        fh_Fibs_vs_Tofd_ac[fi23a]->Fill(x1[tofd1r] * 100. + randx, xFi23a[j] * 100.);
                        fh_Fib_Time_ac[fi23a]->Fill(xFi23a[j] * 100., tStart - tFi23a[j]);

                        cond2 = kTRUE;
                        if (debug3)
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
                                tdet[countdet] = tFi23a[j];
                                countdet++;
                                fFi23a[j] = true;
                            }
                        }
                    }
                }
                for (Int_t j = 0; j < mult23b; j++)
                {
                    if (debug3)
                        cout << "Fib23b: " << j << " y: " << yFi23b[j] << " q: " << qFi23b[j] << endl;
                    // Double_t x31 = -6.8672 * xFi23b[j] * 100. - 27.3507;
                    // if (abs(xFi31[i] * 100. - x31) < dx2 && xFi31[i] > -100 && xFi23b[j] > -100)
                    tempFi23br = false;
                    if (fSimu)
                        tempFi23br = true;
                    if (!fSimu && abs(timeFi31[i] - timeFi23b[j]) < dtft)
                        tempFi23br = true;
                    if (xFi23b[j] * 100. > -6. && tempFi23br)
                    {
                        // if (fGraphCuts && !cut_fi30_fi23b->IsInside(yFi23b[j] * 100., xFi30[i] * 100.))
                        //	continue;
                        fh_Fib31_vs_Fib23b->Fill(xFi23b[j] * 100., xFi31[i] * 100.);
                        fh_Fib31_vs_Fib23b_dx->Fill(xFi23b[j] * 100., xFi31[i] * 100. - xFi23b[j] * 100.);

                        fh_Fib_ToF_ac[fi23b]->Fill(xFi23b[j] * 100., tFi23b[j]);
                        fh_xy_Fib_ac[fi23b]->Fill(xFi23b[j] * 100., yFi23b[j] * 100.);
                        fh_ToT_Fib_ac[fi23b]->Fill(xFi23b[j] * 100., qFi23b[j]);
                        fh_Fibs_vs_Tofd_ac[fi23b]->Fill(x1[tofd1r] * 100. + randx, xFi23b[j] * 100.);
                        fh_Fib_Time_ac[fi23b]->Fill(xFi23b[j] * 100., tStart - tFi23b[j]);

                        cond3 = kTRUE;
                        if (debug3)
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
                                tdet[countdet] = tFi23b[j];
                                countdet++;
                                fFi23b[j] = true;
                            }
                        }
                    }
                }
            }

            for (Int_t i = 0; i < mult30; i++)
            {
                if (debug3)
                    cout << "Fib30: " << i << " x: " << xFi30[i] << " q: " << qFi30[i] << endl;
                Double_t xtemp30 = 0.7796 * x1[tofd1l] * 100. - 28.669;
                if (std::abs(xtemp30 - xFi30[i] * 100.) > dx3)
                    continue;

                fh_Fib_ToF_ac[fi30]->Fill(xFi30[i] * 100., tStart - tFi30[i]);
                fh_xy_Fib_ac[fi30]->Fill(xFi30[i] * 100., yFi30[i] * 100.);
                fh_ToT_Fib_ac[fi30]->Fill(xFi30[i] * 100., qFi30[i]);
                fh_Fibs_vs_Tofd_ac[fi30]->Fill(x1[tofd1l] * 100. + randx, xFi30[i] * 100.);
                fh_Fib_Time_ac[fi30]->Fill(xFi30[i] * 100., tFi30[i]);

                for (Int_t j = 0; j < mult32; j++)
                {
                    if (debug3)
                        cout << "Fib32: " << j << " x: " << xFi32[j] << " q: " << qFi32[j] << endl;
                    Double_t x30 = 0.873 * xFi32[j] * 100. - 3.302; // 1.10926 * xFi30[i] * 100. + 2.8943;
                    tempFi32 = false;
                    if (fSimu)
                        tempFi32 = true;
                    if (!fSimu && abs(timeFi30[i] - timeFi32[j]) < dtft)
                        tempFi32 = true;

                    if ((abs(xFi30[i] * 100. - x30) < dx1 && xFi32[j] * 100. > -30 && xFi30[i] * 100. > -30. &&
                         tempFi32))
                    // if (xFi30[i] > -100 && xFi32[j] > -100)
                    {
                        fh_Fib32_vs_Fib30->Fill(xFi30[i] * 100., xFi32[j] * 100.);
                        fh_Fib32_vs_Fib30_dx->Fill(xFi30[i] * 100., xFi32[j] * 100. - xFi30[i] * 100.);

                        fh_Fib_ToF_ac[fi32]->Fill(xFi32[j] * 100., tFi32[j]);
                        fh_xy_Fib_ac[fi32]->Fill(xFi32[j] * 100., yFi32[j] * 100.);
                        fh_ToT_Fib_ac[fi32]->Fill(xFi32[j] * 100., qFi32[j]);
                        fh_Fibs_vs_Tofd_ac[fi32]->Fill(x1[tofd1l] * 100. + randx, xFi32[j] * 100.);
                        fh_Fib_Time_ac[fi32]->Fill(xFi32[j] * 100., tStart - tFi32[j]);

                        cond4 = kTRUE;
                        if (debug3)
                            cout << "cond4" << endl;

                        if (!maxWerte)
                        {
                            if (!fFi30[i])
                            {
                                detector[countdet] = fi30;
                                xdet[countdet] = xFi30[i];
                                ydet[countdet] = yFi30[i];
                                zdet[countdet] = 0.;
                                qdet[countdet] = qFi30[i];
                                tdet[countdet] = tFi30[i];
                                countdet++;
                                fFi30[i] = true;
                            }
                            if (!fFi32[j])
                            {
                                detector[countdet] = fi32;
                                xdet[countdet] = xFi32[j];
                                ydet[countdet] = yFi32[j];
                                zdet[countdet] = 0.;
                                qdet[countdet] = qFi32[j];
                                tdet[countdet] = tFi32[j];
                                countdet++;
                                fFi32[j] = true;
                            }
                        }
                    }
                }
                for (Int_t j = 0; j < mult23a; j++)
                {
                    if (debug3)
                        cout << "Fib23a: " << j << " x: " << xFi23a[j] << " q: " << qFi23a[j] << endl;
                    Double_t x30 = 6.498 * xFi23a[j] * 100. - 28.93;
                    // if (abs(xFi30[i] * 100. - x30) < dx2 && xFi30[i] > -100 && xFi23a[j] > -100)
                    tempFi23al = false;
                    if (fSimu)
                        tempFi23al = true;
                    if (!fSimu && abs(timeFi30[i] - timeFi23a[j]) < dtft)
                        tempFi23al = true;
                    // if (abs(xFi30[i] * 100. - x30) < dx2 && xFi30[i]*100. > -30 && xFi23a[j]*100. > -6 &&
                    if (xFi23a[j] * 100. > 0.05 && tempFi23al)
                    {
                        // if (fGraphCuts && !cut_fi30_fi23b->IsInside(yFi23b[j] * 100., xFi30[i] * 100.))
                        //	continue;
                        fh_Fib30_vs_Fib23a->Fill(xFi23a[j] * 100., xFi30[i] * 100.);
                        fh_Fib30_vs_Fib23a_dx->Fill(xFi23a[j] * 100., xFi30[i] * 100. - xFi23a[j] * 100.);

                        fh_Fib_ToF_ac[fi23a]->Fill(xFi23a[j] * 100., tFi23a[j]);
                        fh_xy_Fib_ac[fi23a]->Fill(xFi23a[j] * 100., yFi23a[j] * 100.);
                        fh_ToT_Fib_ac[fi23a]->Fill(xFi23a[j] * 100., qFi23a[j]);
                        fh_Fibs_vs_Tofd_ac[fi23a]->Fill(x1[tofd1l] * 100. + randx, xFi23a[j] * 100.);
                        fh_Fib_Time_ac[fi23a]->Fill(xFi23a[j] * 100., tStart - tFi23a[j]);

                        cond5 = kTRUE;
                        if (debug3)
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
                                tdet[countdet] = tFi23a[j];
                                countdet++;
                                fFi23a[j] = true;
                            }
                        }
                    }
                }
                for (Int_t j = 0; j < mult23b; j++)
                {
                    if (debug3)
                        cout << "Fib23b: " << j << " y: " << yFi23b[j] << " q: " << qFi23b[j] << endl;
                    // Double_t x30 = 6.98386 * xFi23b[j] * 100. - 27.39897;
                    // if (abs(xFi30[i] * 100. - x30) < dx2 && xFi30[i] > -100 && xFi23b[j] > -100)
                    tempFi23bl = false;
                    if (fSimu)
                        tempFi23bl = true;
                    if (!fSimu && abs(timeFi30[i] - timeFi23b[j]) < dtft)
                        tempFi23bl = true;
                    if ((yFi23b[j] * 100. > -6. && tempFi23bl))
                    {
                        // if (fGraphCuts && !cut_fi30_fi23b->IsInside(yFi23b[j] * 100., xFi30[i] * 100.))
                        //	continue;
                        fh_Fib30_vs_Fib23b->Fill(xFi23b[j] * 100., xFi30[i] * 100.);
                        fh_Fib30_vs_Fib23b_dx->Fill(xFi23b[j] * 100., xFi30[i] * 100. - xFi23b[j] * 100.);

                        fh_Fib_ToF_ac[fi23b]->Fill(xFi23b[j] * 100., tFi23b[j]);
                        fh_xy_Fib_ac[fi23b]->Fill(xFi23b[j] * 100., yFi23b[j] * 100.);
                        fh_ToT_Fib_ac[fi23b]->Fill(xFi23b[j] * 100., qFi23b[j]);
                        fh_Fibs_vs_Tofd_ac[fi23b]->Fill(x1[tofd1l] * 100. + randx, xFi23b[j] * 100.);
                        fh_Fib_Time_ac[fi23b]->Fill(xFi23b[j] * 100., tStart - tFi23b[j]);

                        cond6 = kTRUE;
                        if (debug3)
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
                                tdet[countdet] = tFi23b[j];
                                countdet++;
                                fFi23b[j] = true;
                            }
                        }
                    }
                }
            }
        }

        if (countdet > 50)
        {
            countdet50 += 1;
            if (debug3)
                cout << "Too many hits!!!" << endl;
            continue;
        }

        // cout << "Test: " << multTofd << endl;
        Bool_t temp_cond = false;
        if ((!fibCuts && ((mult30 > 0 && mult32 > 0 && mult23a > 0 && mult23b > 0) ||
                          (mult31 > 0 && mult33 > 0 && mult23a > 0 && mult23b > 0))) ||
            (fibCuts && ((cond1 && cond2 && cond3) || (cond4 && cond5 && cond6))))
            temp_cond = true;

        if (fB != -1710 && !fPairs)
        {
            Int_t ncount[10] = { 0 };

            // cout<<"Before: "<<endl;
            for (Int_t i = 0; i < countdet; i++)
            {
                ncount[detector[i]] += 1;

                // cout<<"Count: "<< i<<", det: "<<detector[i]<<", qdet: "<<qdet[i]<<endl;
            }

            countdet_s = 0;

            Double_t qdet_s7 = 0., qdet_s6 = 0., qdet_s8 = 0., qdet_s9 = 0.;
            if ((ncount[9] > 0 && ncount[7] > 0 && ncount[2] > 1 && ncount[4] > 1 && ncount[0] > 0 && ncount[1] > 0) ||
                (ncount[6] > 0 && ncount[8] > 0 && ncount[3] > 0 && ncount[5] > 0 && ncount[0] > 0 && ncount[1] > 0))
            // if(1 == 1)
            {
                Bool_t goodQ = true;

                for (Int_t i = 0; i < countdet; i++)
                {
                    Bool_t iffib33 = true;
                    if (detector[i] > 5 && goodQ)
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

                        if (debug)
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
                        qdet_s[countdet_s] = 0.; // qdet[i];
                        tdet_s[countdet_s] = tdet[i];
                        detector_s[countdet_s] = detector[i];

                        nsum[detector_s[countdet_s]] += 1;
                        xdet_sum[detector_s[countdet_s]] += xdet[i];
                        ydet_sum[detector_s[countdet_s]] += ydet[i];
                        zdet_sum[detector_s[countdet_s]] += zdet[i];
                        qdet_sum[detector_s[countdet_s]] += 0.; // qdet[i];
                        tdet_sum[detector_s[countdet_s]] += tdet[i];

                        if (debug)
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
                        qdet_s[countdet_s] = 0.; // qdet[i];
                        tdet_s[countdet_s] = tdet[i];
                        detector_s[countdet_s] = detector[i];

                        nsum[detector_s[countdet_s]] += 1;
                        xdet_sum[detector_s[countdet_s]] += xdet[i];
                        ydet_sum[detector_s[countdet_s]] += ydet[i];
                        zdet_sum[detector_s[countdet_s]] += zdet[i];
                        qdet_sum[detector_s[countdet_s]] += 0.; // qdet[i];
                        tdet_sum[detector_s[countdet_s]] += tdet[i];

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
                // cout<<"************* counter1: "<<counter1<<endl;
                if (counter1 < 1001 && !fAverage)
                {
                    for (Int_t i = 0; i < countdet_s; i++)
                    {
                        // if (debug)
                        cout << "counter1 :" << counter1 << " #" << i << " Det: " << detector_s[i]
                             << " x: " << xdet_s[i] * 100. << " y: " << ydet_s[i] * 100. << " q: " << qdet_s[i] << endl;

                        if (detector_s[i] == 0)
                        {
                            // write fiber detector hits
                            new ((*fFi23aHitItems)[fNofFi23aHitItems++]) R3BFiberMAPMTHitData(
                                0, xdet_s[i] * 100., ydet_s[i] * 100., qdet_s[i], tdet_s[i], 0, 0, 0., 0, 0.);
                            fh_check_QvsX[detector_s[i]]->Fill(xdet_s[i] * 100, qdet_s[i]);
                            fh_check_TvsX[detector_s[i]]->Fill(xdet_s[i] * 100, tdet_s[i]);
                            fh_check_XvsY[detector_s[i]]->Fill(xdet_s[i] * 100, ydet_s[i] * 100);
                        }
                        if (detector_s[i] == 1)
                        {
                            // write fiber detector hits
                            new ((*fFi23bHitItems)[fNofFi23bHitItems++]) R3BFiberMAPMTHitData(
                                0, xdet_s[i] * 100., ydet_s[i] * 100., qdet_s[i], tdet_s[i], 0, 0, 0., 0, 0.);
                            fh_check_QvsX[detector_s[i]]->Fill(ydet_s[i] * 100, qdet_s[i]);
                            fh_check_TvsX[detector_s[i]]->Fill(ydet_s[i] * 100, tdet_s[i]);
                            fh_check_XvsY[detector_s[i]]->Fill(xdet_s[i] * 100, ydet_s[i] * 100);
                        }
                        if (detector_s[i] == 2 && abs(fB) < 1710.)
                        {
                            // write fiber detector hits
                            new ((*fFi30HitItems)[fNofFi30HitItems++]) R3BFiberMAPMTHitData(
                                0, xdet_s[i] * 100., ydet_s[i] * 100., qdet_s[i], tdet_s[i], 0, 0, 0., 0, 0.);
                            fh_check_QvsX[detector_s[i]]->Fill(xdet_s[i] * 100, qdet_s[i]);
                            fh_check_TvsX[detector_s[i]]->Fill(xdet_s[i] * 100, tdet_s[i]);
                            fh_check_XvsY[detector_s[i]]->Fill(xdet_s[i] * 100, ydet_s[i] * 100);
                        }
                        if (detector_s[i] == 3 && abs(fB) > 1710.)
                        {
                            // write fiber detector hits
                            new ((*fFi31HitItems)[fNofFi31HitItems++]) R3BFiberMAPMTHitData(
                                0, xdet_s[i] * 100., ydet_s[i] * 100., qdet_s[i], tdet_s[i], 0, 0, 0., 0, 0.);
                            fh_check_QvsX[detector_s[i]]->Fill(xdet_s[i] * 100, qdet_s[i]);
                            fh_check_TvsX[detector_s[i]]->Fill(xdet_s[i] * 100, tdet_s[i]);
                            fh_check_XvsY[detector_s[i]]->Fill(xdet_s[i] * 100, ydet_s[i] * 100);
                        }
                        if (detector_s[i] == 4 && abs(fB) < 1710.)
                        {
                            // write fiber detector hits
                            new ((*fFi32HitItems)[fNofFi32HitItems++]) R3BFiberMAPMTHitData(
                                0, xdet_s[i] * 100., ydet_s[i] * 100., qdet_s[i], tdet_s[i], 0, 0, 0., 0, 0.);
                            fh_check_QvsX[detector_s[i]]->Fill(xdet_s[i] * 100, qdet_s[i]);
                            fh_check_TvsX[detector_s[i]]->Fill(xdet_s[i] * 100, tdet_s[i]);
                            fh_check_XvsY[detector_s[i]]->Fill(xdet_s[i] * 100, ydet_s[i] * 100);
                        }
                        if (detector_s[i] == 5 && abs(fB) > 1710.)
                        {
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
                }
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

        if (fPairs && temp_cond)
        {
            Int_t ncount[10] = { 0 };
            // cout<<"Before: "<<endl;
            for (Int_t i = 0; i < countdet; i++)
            {
                ncount[detector[i]] += 1;
                // cout<<"Count: "<< i<<", det: "<<detector_s[i]<<", qdet: "<<qdet_s[i]<<endl;
            }

            countdet_s = 0;

            Double_t qdet_s7 = 0., qdet_s6 = 0., qdet_s8 = 0., qdet_s9 = 0.;
            if (ncount[9] == ncount[7] && ncount[6] == ncount[8] && ncount[8] == ncount[9] && ncount[2] > 0 &&
                ncount[3] > 0 && ncount[4] > 0 && ncount[5] > 0 && ncount[0] > 0 && ncount[1] > 0 &&
                ncount[2] == ncount[4] && ncount[3] == ncount[5] && ncount[0] > ncount[2] && ncount[0] > ncount[3] &&
                ncount[1] > ncount[2] && ncount[1] > ncount[3])
            {
                Bool_t goodQ = false;
                for (Int_t i = 0; i < countdet; i++)
                {
                    if (detector[i] == 6)
                        qdet_s6 += qdet[i];
                    if (detector[i] == 7)
                        qdet_s7 += qdet[i];
                    if (detector[i] == 8)
                        qdet_s8 += qdet[i];
                    if (detector[i] == 9)
                        qdet_s9 += qdet[i];
                }
                if (qdet_s6 == qdet_s8 && qdet_s7 == qdet_s9)
                    goodQ = true;

                for (Int_t i = 0; i < countdet; i++)
                {
                    if (detector[i] > 5 && goodQ)
                    {
                        xdet_s[countdet_s] = xdet[i];
                        ydet_s[countdet_s] = ydet[i];
                        zdet_s[countdet_s] = zdet[i];
                        qdet_s[countdet_s] = qdet[i];
                        tdet_s[countdet_s] = tdet[i];
                        detector_s[countdet_s] = detector[i];

                        if (debug)
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
                        if (((detector[i] == 2 || detector[i] == 4) && abs(qdet_s7 - qdet[i] < 2.1)) ||
                            ((detector[i] == 3 || detector[i] == 5) && abs(qdet_s6 - qdet[i] < 2.1)))
                        {
                            xdet_s[countdet_s] = xdet[i];
                            ydet_s[countdet_s] = ydet[i];
                            zdet_s[countdet_s] = zdet[i];
                            qdet_s[countdet_s] = qdet[i];
                            tdet_s[countdet_s] = tdet[i];
                            detector_s[countdet_s] = detector[i];
                            if ((detector[i] == 2 || detector[i] == 4) && qdet_s7 == 2)
                                qdet_s[countdet_s] = 2;
                            if ((detector[i] == 3 || detector[i] == 5) && qdet_s6 == 2)
                                qdet_s[countdet_s] = 2;

                            if (debug)
                            {
                                cout << "fib  det:   " << detector[i] << ", count " << i << ", x: " << xdet[i]
                                     << ", q: " << qdet[i] << endl;
                                cout << "fib  det_s: " << detector_s[countdet_s] << ", count " << countdet_s
                                     << ", x: " << xdet_s[countdet_s] << ", q: " << qdet_s[countdet_s] << endl;
                            }
                            countdet_s++;
                        }
                    }
                    else if (detector[i] < 2 && goodQ)
                    {
                        xdet_s[countdet_s] = xdet[i];
                        ydet_s[countdet_s] = ydet[i];
                        zdet_s[countdet_s] = zdet[i];
                        qdet_s[countdet_s] = qdet[i];
                        tdet_s[countdet_s] = tdet[i];
                        detector_s[countdet_s] = detector[i];
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

                for (Int_t i = 0; i < countdet_s; i++)
                {
                    if (debug)
                        cout << "#" << i << " Det: " << detector_s[i] << " x: " << xdet_s[i] * 100.
                             << " y: " << ydet_s[i] * 100. << " q: " << qdet_s[i] << " t: " << tdet_s[i] << endl;

                    if (detector_s[i] == 0)
                    {
                        // write fiber detector hits
                        new ((*fFi23aHitItems)[fNofFi23aHitItems++]) R3BFiberMAPMTHitData(
                            0, xdet_s[i] * 100., ydet_s[i] * 100., qdet_s[i], tdet_s[i], 0, 0, 0., 0, 0.);
                        fh_check_QvsX[detector_s[i]]->Fill(xdet_s[i] * 100, qdet_s[i]);
                        fh_check_TvsX[detector_s[i]]->Fill(xdet_s[i] * 100, tdet_s[i]);
                        fh_check_XvsY[detector_s[i]]->Fill(xdet_s[i] * 100, ydet_s[i] * 100);
                    }
                    if (detector_s[i] == 1)
                    {
                        // write fiber detector hits
                        new ((*fFi23bHitItems)[fNofFi23bHitItems++]) R3BFiberMAPMTHitData(
                            0, xdet_s[i] * 100., ydet_s[i] * 100., qdet_s[i], tdet_s[i], 0, 0, 0., 0, 0.);
                        fh_check_QvsX[detector_s[i]]->Fill(ydet_s[i] * 100, qdet_s[i]);
                        fh_check_TvsX[detector_s[i]]->Fill(ydet_s[i] * 100, tdet_s[i]);
                        fh_check_XvsY[detector_s[i]]->Fill(xdet_s[i] * 100, ydet_s[i] * 100);
                    }
                    if (detector_s[i] == 2)
                    {
                        // write fiber detector hits
                        new ((*fFi30HitItems)[fNofFi30HitItems++]) R3BFiberMAPMTHitData(
                            0, xdet_s[i] * 100., ydet_s[i] * 100., qdet_s[i], tdet_s[i], 0, 0, 0., 0, 0.);
                        fh_check_QvsX[detector_s[i]]->Fill(xdet_s[i] * 100, qdet_s[i]);
                        fh_check_TvsX[detector_s[i]]->Fill(xdet_s[i] * 100, tdet_s[i]);
                        fh_check_XvsY[detector_s[i]]->Fill(xdet_s[i] * 100, ydet_s[i] * 100);
                    }
                    if (detector_s[i] == 3)
                    {
                        // write fiber detector hits
                        new ((*fFi31HitItems)[fNofFi31HitItems++]) R3BFiberMAPMTHitData(
                            0, xdet_s[i] * 100., ydet_s[i] * 100., qdet_s[i], tdet_s[i], 0, 0, 0., 0, 0.);
                        fh_check_QvsX[detector_s[i]]->Fill(xdet_s[i] * 100, qdet_s[i]);
                        fh_check_TvsX[detector_s[i]]->Fill(xdet_s[i] * 100, tdet_s[i]);
                        fh_check_XvsY[detector_s[i]]->Fill(xdet_s[i] * 100, ydet_s[i] * 100);
                    }
                    if (detector_s[i] == 4)
                    {
                        // write fiber detector hits
                        new ((*fFi32HitItems)[fNofFi32HitItems++]) R3BFiberMAPMTHitData(
                            0, xdet_s[i] * 100., ydet_s[i] * 100., qdet_s[i], tdet_s[i], 0, 0, 0., 0, 0.);
                        fh_check_QvsX[detector_s[i]]->Fill(xdet_s[i] * 100, qdet_s[i]);
                        fh_check_TvsX[detector_s[i]]->Fill(xdet_s[i] * 100, tdet_s[i]);
                        fh_check_XvsY[detector_s[i]]->Fill(xdet_s[i] * 100, ydet_s[i] * 100);
                    }
                    if (detector_s[i] == 5)
                    {
                        // write fiber detector hits
                        new ((*fFi33HitItems)[fNofFi33HitItems++]) R3BFiberMAPMTHitData(
                            0, xdet_s[i] * 100., ydet_s[i] * 100., qdet_s[i], tdet_s[i], 0, 0, 0., 0, 0.);
                        fh_check_QvsX[detector_s[i]]->Fill(xdet_s[i] * 100, qdet_s[i]);
                        fh_check_TvsX[detector_s[i]]->Fill(xdet_s[i] * 100, tdet_s[i]);
                        fh_check_XvsY[detector_s[i]]->Fill(xdet_s[i] * 100, ydet_s[i] * 100);
                    }
                    if ((detector_s[i] == 6 || detector_s[i] == 8 || detector_s[i] == 7 || detector_s[i] == 9))
                    {
                        if (ftrackerType == 0) // TofD data written out for each plane separately
                        {
                            Int_t ipl;
                            if (detector_s[i] == 6 || detector_s[i] == 7)
                                ipl = 1;
                            else
                                ipl = 2;
                            new ((*fTofdHitItems)[fNofTofdHitItems++]) R3BTofdHitData(
                                tdet_s[i], xdet_s[i] * 100., ydet_s[i] * 100., qdet_s[i], -5, qdet_s[i], ipl, 1, 0);

                            fh_check_QvsX[6]->Fill(xdet_s[i] * 100, qdet_s[i]);
                            fh_check_TvsX[6]->Fill(xdet_s[i] * 100, tdet_s[i]);
                            fh_check_XvsY[6]->Fill(xdet_s[i] * 100, ydet_s[i] * 100);
                        }
                        else // first and second plane will be written as one detector; this makes Dima's tracker faster
                        {
                            if (qdet_s[i] == 8)
                            {
                                ax += xdet_s[i];
                                ay += ydet_s[i];
                                at += tdet_s[i];
                                aq += qdet_s[i];
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

                    //    cout<<"TOFD C: " << " x: " << cx * 100.
                    //           << " y: " << cy * 100. << " q: " << cq << " t: " <<  endl;
                    //   cout<<"TOFD He: " << " x: " << ax * 100.
                    //         << " y: " << ay * 100. << " q: " << aq << " t: " <<  endl;

                    new ((*fTofdHitItems)[fNofTofdHitItems++])
                        R3BTofdHitData(at, ax * 100., ay * 100., aq, -5., aq, 10, 1, 0);

                    // Fill check spectra:
                    fh_check_QvsX[6]->Fill(cx * 100, cq);
                    fh_check_TvsX[6]->Fill(cx * 100, ct);
                    fh_check_XvsY[6]->Fill(cx * 100, cy * 100);
                    fh_check_QvsX[6]->Fill(ax * 100, aq);
                    fh_check_TvsX[6]->Fill(ax * 100, at);
                    fh_check_XvsY[6]->Fill(ax * 100, ay * 100);
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
            // if (debug)
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
    cout << "Events: " << fNEvents << endl;
    cout << "Events after selection: " << fNeventselect << endl;
    cout << "Wrong Trigger: " << counterWrongTrigger << endl;
    cout << "Wrong Tpat: " << counterWrongTpat << endl;
    cout << "ROLU veto: " << counterRolu << endl;
    cout << "Califa veto: " << counterCalifa << endl;
    cout << "TofD: " << counterTofd << endl;
    cout << "TofD multi: " << counterTofdMulti << endl;
    cout << "Selected events: " << counter1 << endl;
    cout << "Hits with countddet>50: " << countdet50 << endl;

    cout << "Hits TofD " << hits1 << endl;

    cout << "Eff. Fi30 min: " << hits30 << "  " << hits30 / hits1 << endl;
    cout << "Eff. Fi31 min: " << hits31 << "  " << hits31 / hits1 << endl;
    cout << "Eff. Fi32 min: " << hits32 << "  " << hits32 / hits1 << endl;
    cout << "Eff. Fi33 min: " << hits33 << "  " << hits33 / hits1 << endl;

    cout << "Eff. Fi30 max: " << hits30bc << "  " << hits30bc / hits1 << endl;
    cout << "Eff. Fi31 max: " << hits31bc << "  " << hits31bc / hits1 << endl;
    cout << "Eff. Fi32 max: " << hits32bc << "  " << hits32bc / hits1 << endl;
    cout << "Eff. Fi33 max: " << hits33bc << "  " << hits33bc / hits1 << endl;

    fh_Tpat->Write();
    fh_Trigger->Write();
    if (fMappedItems.at(DET_BMON))
    {
        fh_TOFDOR->Write();
        fh_SEE->Write();
        fh_IC->Write();
    }

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
        }
    }

    fh_Fib33_vs_Fib31->Write();
    fh_Fib33_vs_Fib31_dx->Write();
    fh_Fib31_vs_Fib23a->Write();
    fh_Fib31_vs_Fib23a_dx->Write();
    fh_Fib32_vs_Fib30->Write();
    fh_Fib32_vs_Fib30_dx->Write();
    fh_Fib30_vs_Fib23b->Write();
    fh_Fib30_vs_Fib23b_dx->Write();
    fh_Fib30_vs_Fib23a->Write();
    fh_Fib30_vs_Fib23a_dx->Write();
    fh_Fib31_vs_Fib23b->Write();
    fh_Fib31_vs_Fib23b_dx->Write();

    for (int i = 0; i < 7; i++)
    {
        fh_check_QvsX[i]->Write();
        fh_check_TvsX[i]->Write();
        fh_check_XvsY[i]->Write();
    }
}

ClassImp(R3BPreTrackS494)
