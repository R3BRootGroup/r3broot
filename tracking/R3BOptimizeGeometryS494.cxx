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

/******************************************************************************
 *****                   Task for fitting a fragment                      *****
 *****     Created Nov 2016 by D. Kresan                                  *****
 *****     Comments to the algorithm implementation from R. Plag          *****
 *****     Implemented ideal fitter, July 2017, D. Kresan                 *****
 ******************************************************************************/

#include "R3BOptimizeGeometryS494.h"
#include "FairIon.h"
#include "G4NistManager.hh"
#include "R3BFi4HitItem.h"
#include "R3BFragmentFitterGeneric.h"
#include "R3BGladFieldMap.h"
#include "R3BHit.h"
#include "R3BMCTrack.h"
#include "R3BPspPoint.h"
#include "R3BTGeoPar.h"
#include "R3BTPropagator.h"
#include "R3BTofdPoint.h"
#include "R3BTrackingDetector.h"
#include "R3BTrackingParticle.h"
#include "R3BTrackingSetup.h"

#include "R3BTofdCalData.h"
#include "R3BTofdHitData.h"
#include "R3BTofdMappedData.h"

#include "R3BFiberMAPMTCalData.h"
#include "R3BFiberMAPMTHitData.h"
#include "R3BFiberMAPMTMappedData.h"

#include "FairLogger.h"
#include "FairRootManager.h"
#include "FairRunAna.h"
#include "FairRuntimeDb.h"

#include "TArc.h"
#include "TClonesArray.h"
#include "TF1.h"
#include "TGraphErrors.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TMath.h"
#include "TRandom.h"

#include <fstream>
#include <iomanip>
#include <iostream>

using namespace std;

#define SPEED_OF_LIGHT 29.9792458 // cm/ns
//#define Amu 0.938272
//#define Fair_Amu 0.931494028
namespace GlobalChi2
{
    R3BTrackingSetup* gSetup;
    R3BOptimizeGeometryS494* gThisTask;
} // namespace GlobalChi2

using namespace GlobalChi2;

R3BOptimizeGeometryS494::R3BOptimizeGeometryS494(const char* name, Bool_t vis, Int_t verbose)
    : FairTask(name, verbose)
    , fFieldPar(NULL)
    , fPropagator(NULL)
    , fArrayMCTracks(NULL)
    , fDetectors(new R3BTrackingSetup())
    , fDetectorsLeft(new R3BTrackingSetup())
    , fDetectorsRight(new R3BTrackingSetup())
    , fArrayFragments(new TClonesArray("R3BTrackingParticle"))
    , candidate(new R3BTrackingParticle())
    , fNEvents(0)
    , fVis(vis)
    , fForward(kTRUE)
    , fBfield(-1710.)
    , fFitter(nullptr)
    , fEnergyLoss(kTRUE)
    , fSimu(kFALSE)
{
    // this is the list of detectors (active areas) we use for tracking
    /*
        fDetectorsLeft->AddDetector("target", kTarget, "TargetGeoPar");
        fDetectorsLeft->AddDetector("fi23a", kTargetGlad, "fi23aGeoPar", "Fi23aHit");
        fDetectorsLeft->AddDetector("fi23b", kTargetGlad, "fi23bGeoPar", "Fi23bHit");
        fDetectorsLeft->AddDetector("fi30", kAfterGlad, "fi30GeoPar", "Fi30Hit");
        fDetectorsLeft->AddDetector("fi32", kAfterGlad, "fi32GeoPar", "Fi32Hit");
        fDetectorsLeft->AddDetector("tofd", kTof, "tofdGeoPar", "TofdHit");

        fDetectorsRight->AddDetector("target", kTarget, "TargetGeoPar");
        fDetectorsRight->AddDetector("fi23a", kTargetGlad, "fi23aGeoPar", "Fi23aHit");
        fDetectorsRight->AddDetector("fi23b", kTargetGlad, "fi23bGeoPar", "Fi23bHit");
        fDetectorsRight->AddDetector("fi31", kAfterGlad, "fi31GeoPar", "Fi31Hit");
        fDetectorsRight->AddDetector("fi33", kAfterGlad, "fi33GeoPar", "Fi33Hit");
        fDetectorsRight->AddDetector("tofd", kTof, "tofdGeoPar", "TofdHit");
    */
    fDetectors->AddDetector("target", kTarget, "TargetGeoPar");
    fDetectors->AddDetector("fi23a", kTargetGlad, "fi23aGeoPar", "Fi23aHit");
    fDetectors->AddDetector("fi23b", kTargetGlad, "fi23bGeoPar", "Fi23bHit");
    fDetectors->AddDetector("fi31", kAfterGlad, "fi31GeoPar", "Fi31Hit");
    fDetectors->AddDetector("fi30", kAfterGlad, "fi30GeoPar", "Fi30Hit");
    fDetectors->AddDetector("fi33", kAfterGlad, "fi33GeoPar", "Fi33Hit");
    fDetectors->AddDetector("fi32", kAfterGlad, "fi32GeoPar", "Fi32Hit");
    fDetectors->AddDetector("tofd", kTof, "tofdGeoPar", "TofdHit");
}

R3BOptimizeGeometryS494::~R3BOptimizeGeometryS494() {}

/* For the tracking we use a user-defined list of TrackingDetectors,
 * stored in a TClonesArrays. The TrackingDetectors will provide
 * detector properties as well as hit coordinates to which we want to
 * fit the track.
 *
 * The tracker also needs a (measured) time-of-flight. Duh!
 * And a charge.
 * The charge is part of the TOF-wall hit level (or at least: will be)
 * The time-of-flight can be constructed from the hits of
 * start and stop detector, aka LOS and TOFD.
 *
 * Unfortunately, the HIT level data is not universal, hence, we cannot use the
 * same generic code for all detectors. We *could* have X,Y,C,T as basic structure
 * for ALL hit data. But not all detectors provide all variables.
 * We could inherit getters from a common base class and return NAN in the
 * overloaded members. Can fix this later. For now we can manually construct
 * a list of detectors with X,Y and additionaly fetch TOF and charge.
 *
 * Btw: the old tracker required the user to provide the measured positions
 * for the track fitting. That was actually quite convenient. Can this
 * be done in the macro? Nope. Because there is no event loop in the macro.
 * Pity.
 *
 * So we need to do it here in the Exec function. Like:
 * for each psp1_hit
 *   for each psp2_hit
 *     for each mtof_hit
 *          track_fragment()
 *
 * Also, we need to know the materials of the detectors and a way to        !! discuss !! ToDo
 * calc the energy loss (forward and backward).
 *
 */
InitStatus R3BOptimizeGeometryS494::Init()
{
    cout << setprecision(7);
    fSide = 3; // 0 = front, 1 = left, 2 = right
    FairRootManager* man = FairRootManager::Instance();

    fArrayMCTracks = (TClonesArray*)man->GetObject("MCTrack");
    if (NULL == fArrayMCTracks)
    {
        LOG(error) << "No MC Track array found in input file.";
        //        return kERROR;
    }

    assert(DET_MAX + 1 == sizeof(fDetectorNames) / sizeof(fDetectorNames[0]));
    printf("Have %d fiber detectors.\n", NOF_FIB_DET);
    for (int det = 0; det < DET_MAX; det++)
    {
        fArrayHits.push_back((TClonesArray*)man->GetObject(Form("%sHit", fDetectorNames[det])));
        // if (det == DET_MAX - 1)
        //    maxevent = man->CheckMaxEventNo();
        if (NULL == fArrayHits.at(det))
        {
            printf("Could not find hit data for '%s'.\n", fDetectorNames[det]);
        }
        else
        {
            printf("Found hit data for '%s'.\n", fDetectorNames[det]);
        }
    }

    cout << "[INFO] Finished with getting data" << endl;

    man->Register("TrackingParticle", "Tracking", fArrayFragments, kTRUE);

    cout << "[INFO] Finished with registring tracking particle" << endl;

    if (!InitPropagator())
    {
        LOG(error) << "!InitPropagator";
        return kERROR;
    }

    cout << "[INFO] Propagator init" << endl;

    // fDetectorsLeft->Init();
    // fDetectorsRight->Init();
    fDetectors->Init();

    cout << "[INFO] Finished with detectors init" << endl;

    // For the moment we have to change manually between the setup
    gSetup = fDetectors;
    // gSetup = fDetectorsLeft;
    // gSetup = fDetectorsRight;

    cout << "[INFO] Finished with gsetup" << endl;

    gThisTask = this;

    fFitter->Init(fPropagator, fEnergyLoss);

    fh_chi2 = new TH1F("h_chi2", "Chi2", 10000, 0., 1000.);
    fh_Erel = new TH1F("h_Erel", "Erel / MeV", 600, -10., 50);
    fh_psum = new TH1F("h_psum", "psum / MeV/c", 500, 16000., 18000.);

    cout << "[INFO] FInished with init" << endl;

    return kSUCCESS;
}

InitStatus R3BOptimizeGeometryS494::ReInit()
{
    SetParContainers();

    if (!InitPropagator())
    {
        return kERROR;
    }
    fFitter->Init(fPropagator, fEnergyLoss);
    return kSUCCESS;
}

void R3BOptimizeGeometryS494::SetParContainers()
{
    fFieldPar = (R3BFieldPar*)FairRuntimeDb::instance()->getContainer("R3BFieldPar");

    // fDetectorsLeft->SetParContainers();
    // fDetectorsRight->SetParContainers();
    fDetectors->SetParContainers();
}

void R3BOptimizeGeometryS494::Exec(const Option_t*)
{

    // For the moment one has to swap manually between the setups
    /*
        //R3BTrackingDetector* target = fDetectorsLeft->GetByName("target");
        R3BTrackingDetector* fi23a = fDetectorsLeft->GetByName("fi23a");
        R3BTrackingDetector* fi23b = fDetectorsLeft->GetByName("fi23b");
        R3BTrackingDetector* fi30 = fDetectorsLeft->GetByName("fi30");
        R3BTrackingDetector* fi31 = fDetectorsRight->GetByName("fi31");
        R3BTrackingDetector* fi32 = fDetectorsLeft->GetByName("fi32");
        R3BTrackingDetector* fi33 = fDetectorsRight->GetByName("fi33");
        R3BTrackingDetector* tof = fDetectorsLeft->GetByName("tofd");
    */

    R3BTrackingDetector* target = fDetectors->GetByName("target");
    R3BTrackingDetector* fi23a = fDetectors->GetByName("fi23a");
    R3BTrackingDetector* fi23b = fDetectors->GetByName("fi23b");
    R3BTrackingDetector* fi31 = fDetectors->GetByName("fi31");
    R3BTrackingDetector* fi30 = fDetectors->GetByName("fi30");
    R3BTrackingDetector* fi33 = fDetectors->GetByName("fi33");
    R3BTrackingDetector* fi32 = fDetectors->GetByName("fi32");
    R3BTrackingDetector* tof = fDetectors->GetByName("tofd");

    /*
        cout << "tof: " << tof->fArrayHits->GetEntriesFast() << endl;
        cout << "fi23a: " << fi23a->fArrayHits->GetEntriesFast() << endl;
        cout << "fi23b: " << fi23b->fArrayHits->GetEntriesFast() << endl;
        cout << "fi30: " << fi30->fArrayHits->GetEntriesFast() << endl;
        cout << "fi31: " << fi31->fArrayHits->GetEntriesFast() << endl;
        cout << "fi32: " << fi32->fArrayHits->GetEntriesFast() << endl;
        cout << "fi33: " << fi33->fArrayHits->GetEntriesFast() << endl;
    */

    if (tof->fArrayHits->GetEntriesFast() < 2)
        return;
    if (fi23a->fArrayHits->GetEntriesFast() < 2)
        return;
    if (fi23b->fArrayHits->GetEntriesFast() < 2)
        return;
    if (fi30->fArrayHits->GetEntriesFast() + fi32->fArrayHits->GetEntriesFast() + fi31->fArrayHits->GetEntriesFast() +
            fi33->fArrayHits->GetEntriesFast() <
        4)
        return;

    // fDetectorsLeft->CopyToBuffer();
    if (fNEvents > 1000)
        return;
    // if(fNEvents >  10) return;

    fDetectors->CopyToBuffer();
    cout << "Saved events: " << gThisTask->GetNEvents() << endl;

    fNEvents += 1;
}

double Chi2AllEvents(const double* xx)
{

    // For the moment ...
    Int_t fSide = 3;
    Double_t chi2_all = 0;
    Double_t fi23x, fi23y, fi23z, fi30x, fi30z, fi30a, fi31x, fi31z, fi31a, fi32x, fi32z, fi32a, fi33x, fi33z, fi33a,
        tofdx, tofdy, tofdz, zdif;

    Int_t nof = 0;
    cout.precision(7);
    /*
    cout << "new correction: " << xx[0] << "  " << xx[1] << "  " << xx[2] << "  " << xx[3] << "  " << xx[4] << "  "
         << xx[5] << "  " << xx[6] << "  " << xx[7] << "  " << xx[8] << "  " << xx[9] << "  " << xx[10] << "  "
         << xx[11] << "  " << xx[12] << "  " << xx[13] << "  " << xx[14] << "  " << xx[15] << "  " << endl;
    */
    fi23x = 0.00000000;
    fi23y = 0.0000000;
    // fi23z = 91.2;
    fi23z = 91.04;

    fi30x = -62.18069;
    fi30z = 570.59330;
    // fi30a = -13.68626;
    fi30a = -14.;

    fi31x = -112.67330;
    fi31z = 537.9038;
    // fi31a = -193.8265;
    fi31a = -194.;

    fi32x = -82.41713;
    fi32z = 632.9688;
    // fi32a = -14.02967;
    fi32a = -14.;

    fi33x = -131.8124;
    fi33z = 597.7741;
    // fi33a = -193.8699;
    fi33a = -194.;

    tofdx = -129.300;
    tofdy = 0.0;
    tofdz = 685.4;

    if (fSide == 3) // infront
    {
        fi23x = xx[0];
        fi23y = xx[1];
        fi23z = xx[2]; //*10.;
        fi30x = xx[3]; //*10.;
        fi30z = xx[4]; //*100.;
        fi32x = xx[5]; //*10.;
        fi31x = xx[6]; //*100.;
        fi31z = xx[7]; //*100.;
        fi33x = xx[8]; //*100.;
        tofdx = xx[9]; //*100.;
        tofdy = xx[10];
        tofdz = xx[11]; //*100.;
        fi30a = xx[12]; //*10.;
        fi31a = xx[13]; //*100.;
        fi32a = xx[14]; //*10.;
        fi33a = xx[15]; //*100.;
        // fi32a = xx[12];//*10.;
        // fi33a = xx[13];//*100.;
    }

    // if(!fSimu)
    {
        fi32z = fi30z + 62.77813349;
        fi33z = fi31z + 62.77813349;
    }
    // else
    {
        // simu:
        // fi32z = fi30z + 632.9688-570.59330;
        // fi33z = fi31z + 597.7741-537.9038;
    }

    // for break-up run we need this line
    //((R3BGladFieldMap*)FairRunAna::Instance()->GetField())->SetTrackerCorrectionScale(xx[6]);

    //((R3BGladFieldMap*)FairRunAna::Instance()->GetField())->SetTrackerCorrectionScale(xx[15]);
    // Double_t scale = ((R3BGladFieldMap *) FairRunAna::Instance()->GetField())->GetScale();
    // Double_t field = ((R3BGladFieldMap *) FairRunAna::Instance()->GetField())->GetBy(0.,0.,240.);
    // cout << "Field:" << field << " scale: " << scale << endl;

    R3BTrackingDetector* fi23a = gSetup->GetByName("fi23a");

    fi23a->pos0 = TVector3(0., 0., 0.);
    fi23a->pos1 = TVector3(5.42976, 5., 0.);
    fi23a->pos2 = TVector3(-5.42976, 5., 0.);

    fi23a->pos0.RotateY(0. * TMath::DegToRad());
    fi23a->pos1.RotateY(0. * TMath::DegToRad());
    fi23a->pos2.RotateY(0. * TMath::DegToRad());

    TVector3 trans3a(fi23x, 0., fi23z - 0.015);

    fi23a->pos0 += trans3a;
    fi23a->pos1 += trans3a;
    fi23a->pos2 += trans3a;
    fi23a->norm = ((fi23a->pos1 - fi23a->pos0).Cross(fi23a->pos2 - fi23a->pos0)).Unit();

    // fi23a->pos0.Print();
    //  fi23a->pos1.Print();
    // fi23a->pos2.Print();

    R3BTrackingDetector* fi23b = gSetup->GetByName("fi23b");

    fi23b->pos0 = TVector3(0., 0., 0.);
    fi23b->pos1 = TVector3(5., 5.42976, 0.);
    fi23b->pos2 = TVector3(-5., 5.42976, 0.);

    fi23b->pos0.RotateY(0. * TMath::DegToRad());
    fi23b->pos1.RotateY(0. * TMath::DegToRad());
    fi23b->pos2.RotateY(0. * TMath::DegToRad());

    TVector3 trans3b(0., fi23y, fi23z + 0.015);

    fi23b->pos0 += trans3b;
    fi23b->pos1 += trans3b;
    fi23b->pos2 += trans3b;
    fi23b->norm = ((fi23b->pos1 - fi23b->pos0).Cross(fi23b->pos2 - fi23b->pos0)).Unit();

    // fi23b->pos0.Print();
    // fi23b->pos1.Print();
    // fi23b->pos2.Print();

    R3BTrackingDetector* fi30 = gSetup->GetByName("fi30");
    fi30->pos0 = TVector3(0., 0., 0.);
    fi30->pos1 = TVector3(26.7351, 25., 0.);
    fi30->pos2 = TVector3(-26.7351, 25., 0.);

    fi30->pos0.RotateY(fi30a * TMath::DegToRad());
    fi30->pos1.RotateY(fi30a * TMath::DegToRad());
    fi30->pos2.RotateY(fi30a * TMath::DegToRad());

    TVector3 trans30(fi30x, 0., fi30z);

    fi30->pos0 += trans30;
    fi30->pos1 += trans30;
    fi30->pos2 += trans30;
    fi30->norm = ((fi30->pos1 - fi30->pos0).Cross(fi30->pos2 - fi30->pos0)).Unit();

    // fi30->pos0.Print();
    // fi30->pos1.Print();
    // fi30->pos2.Print();

    R3BTrackingDetector* fi31 = gSetup->GetByName("fi31");

    fi31->pos0 = TVector3(0., 0., 0.);
    fi31->pos1 = TVector3(26.7092, 25., 0.);
    fi31->pos2 = TVector3(-26.7092, 25., 0.);

    fi31->pos0.RotateY(fi31a * TMath::DegToRad());
    fi31->pos1.RotateY(fi31a * TMath::DegToRad());
    fi31->pos2.RotateY(fi31a * TMath::DegToRad());
    TVector3 trans31(fi31x, 0., fi31z);

    fi31->pos0 += trans31;
    fi31->pos1 += trans31;
    fi31->pos2 += trans31;
    fi31->norm = ((fi31->pos1 - fi31->pos0).Cross(fi31->pos2 - fi31->pos0)).Unit();

    // fi31->pos0.Print();
    // fi31->pos1.Print();
    // fi31->pos2.Print();

    R3BTrackingDetector* fi32 = gSetup->GetByName("fi32");

    fi32->pos0 = TVector3(0., 0., 0.);
    fi32->pos1 = TVector3(26.4765, 25., 0.);
    fi32->pos2 = TVector3(-26.4765, 25., 0.);

    fi32->pos0.RotateY(fi32a * TMath::DegToRad());
    fi32->pos1.RotateY(fi32a * TMath::DegToRad());
    fi32->pos2.RotateY(fi32a * TMath::DegToRad());

    TVector3 trans32(fi32x, 0., fi32z);

    fi32->pos0 += trans32;
    fi32->pos1 += trans32;
    fi32->pos2 += trans32;
    fi32->norm = ((fi32->pos1 - fi32->pos0).Cross(fi32->pos2 - fi32->pos0)).Unit();

    // fi32->pos0.Print();
    // fi32->pos1.Print();
    // fi32->pos2.Print();

    R3BTrackingDetector* fi33 = gSetup->GetByName("fi33");
    fi33->pos0 = TVector3(0., 0., 0.);
    fi33->pos1 = TVector3(26.5024, 25., 0.);
    fi33->pos2 = TVector3(-26.5024, 25., 0.);

    fi33->pos0.RotateY(fi33a * TMath::DegToRad());
    fi33->pos1.RotateY(fi33a * TMath::DegToRad());
    fi33->pos2.RotateY(fi33a * TMath::DegToRad());

    TVector3 trans33(fi33x, 0., fi33z);

    fi33->pos0 += trans33;
    fi33->pos1 += trans33;
    fi33->pos2 += trans33;
    fi33->norm = ((fi33->pos1 - fi33->pos0).Cross(fi33->pos2 - fi33->pos0)).Unit();

    // fi33->pos0.Print();
    // fi33->pos1.Print();
    // fi33->pos2.Print();

    R3BTrackingDetector* tofd = gSetup->GetByName("tofd");
    tofd->pos0 = TVector3(0., 0., 0.);
    tofd->pos1 = TVector3(60., 40., 0.);
    tofd->pos2 = TVector3(-60., 40., 0.);

    tofd->pos0.RotateY(-18. * TMath::DegToRad());
    tofd->pos1.RotateY(-18. * TMath::DegToRad());
    tofd->pos2.RotateY(-18. * TMath::DegToRad());

    TVector3 transtofd(tofdx, tofdy, tofdz);

    tofd->pos0 += transtofd;
    tofd->pos1 += transtofd;
    tofd->pos2 += transtofd;
    tofd->norm = ((tofd->pos1 - tofd->pos0).Cross(tofd->pos2 - tofd->pos0)).Unit();

    // tofd->pos0.Print();
    // tofd->pos1.Print();
    // tofd->pos2.Print();

    for (Int_t iev = 0; iev < gThisTask->GetNEvents(); iev++)
    {
        gSetup->TakeHitsFromBuffer(iev);
        Double_t fieldScale;

        // für break-up events wieder rein

        //((R3BGladFieldMap*)FairRunAna::Instance()->GetField())->SetTrackerCorrectionScale(xx[15]);
        // Double_t scale = ((R3BGladFieldMap *) FairRunAna::Instance()->GetField())->GetScale();
        // Double_t field = ((R3BGladFieldMap *) FairRunAna::Instance()->GetField())->GetBy(0.,0.,240.);
        // cout << "Field:" << field << " scale: " << scale << endl;

        Double_t chi2 = gThisTask->Chi2();
        if (chi2 > 0)
        {
            // cout << "For event: "<<iev+1<<" Chi2 value: " << chi2 << endl;
            // cout<<"  "<<endl;
            chi2_all += chi2;
            nof++;
        }
    }

    // cout << "Test: " << nof << endl;
    chi2_all = chi2_all / (double)nof;

    if (nof < 1)
        chi2_all = 1.e6;
    cout << " " << endl;
    cout << "With parameters: " << xx[0] << ", " << xx[1] << ", " << xx[2] << ", " << xx[3] << ", " << xx[4] << ", "
         << xx[5] << ", " << xx[6] << ", " << xx[7] << ", " << xx[8] << ", " << xx[9] << ", " << xx[10] << ", "
         << xx[11] << ", " << xx[12] << ", " << xx[13] << ", " << xx[14] << ", " << xx[15] << endl;
    cout << "Chi2 all: " << chi2_all << " for nof: " << nof << endl;
    return chi2_all;
}

double R3BOptimizeGeometryS494::Chi2()
{

    fArrayFragments->Clear();
    if (fFragments.size() > 0)
    {
        for (auto const& x : fFragments)
        {
            delete x;
        }
        fFragments.clear();
    }
    Bool_t debug = false;

    /* this part needs to be adopted to each experiment / setup
     *
     * First thing is to tell each detector which hit we want to use
     * for the track fitting. This means: manually copy the hit we want to
     * use for tracking from the TClonesArray into the single hit structure.
     *
     * Then call fit_fragment() to optimize the track and calculate
     * particle properties.
     */

    R3BTrackingDetector* target = gSetup->GetByName("target");
    R3BTrackingDetector* fi23a = gSetup->GetByName("fi23a");
    R3BTrackingDetector* fi23b = gSetup->GetByName("fi23b");
    R3BTrackingDetector* fi31 = gSetup->GetByName("fi31");
    R3BTrackingDetector* fi30 = gSetup->GetByName("fi30");
    R3BTrackingDetector* fi33 = gSetup->GetByName("fi33");
    R3BTrackingDetector* fi32 = gSetup->GetByName("fi32");
    R3BTrackingDetector* tof = gSetup->GetByName("tofd");

    //  if (target->hits.size() < 1)
    //    target->hits.push_back(new R3BHit(0, 0.0, 0.0, 0., 0., 0));

    if (tof->hits.size() > 0 && debug)
    {
        cout << "*************** NEW EVENT ****" << endl;
        cout << "Hits ToFD: " << tof->hits.size() << endl;
        cout << "Hits right: " << fi23a->hits.size() << "  " << fi23b->hits.size() << "  " << fi31->hits.size() << "  "
             << fi33->hits.size() << endl;
        cout << "Hits left: " << fi23a->hits.size() << "  " << fi23b->hits.size() << "  " << fi32->hits.size() << "  "
             << fi30->hits.size() << endl;

        for (Int_t i = 0; i < fi23a->hits.size(); i++)
        {
            cout << "Fib23a hits: " << i << ", " << fi23a->hits.at(i)->GetX() << endl;
        }
        for (Int_t i = 0; i < fi23b->hits.size(); i++)
        {
            cout << "Fib23b hits: " << i << ", " << fi23b->hits.at(i)->GetY() << endl;
        }
        for (Int_t i = 0; i < fi30->hits.size(); i++)
        {
            cout << "Fib30 hits: " << i << ", " << fi30->hits.at(i)->GetX() << endl;
        }
        for (Int_t i = 0; i < fi31->hits.size(); i++)
        {
            cout << "Fib31 hits: " << i << ", " << fi31->hits.at(i)->GetX() << endl;
        }
        for (Int_t i = 0; i < fi32->hits.size(); i++)
        {
            cout << "Fib32 hits: " << i << ", " << fi32->hits.at(i)->GetX() << endl;
        }
        for (Int_t i = 0; i < fi33->hits.size(); i++)
        {
            cout << "Fib33 hits: " << i << ", " << fi33->hits.at(i)->GetX() << endl;
        }
        for (Int_t i = 0; i < tof->hits.size(); i++)
        {
            cout << "Tofd hits  : " << i << ", " << tof->hits.at(i)->GetX() << ", " << tof->hits.at(i)->GetY() << endl;
        }
    }

    // Start values
    Double_t beta = 0.0;
    Double_t beta0 = 0.0;
    Double_t x0 = 0.0;
    Double_t y0 = 0.0;
    Double_t z0 = 0.0;
    Double_t px0 = 0.0;
    Double_t py0 = 0.0;
    Double_t pz0 = 0.0;
    Double_t p = 0.0;
    Double_t mass = 0.0;

    Double_t betaHe = 0.0;
    Double_t x0He = 0.0;
    Double_t y0He = 0.0;
    Double_t z0He = 0.0;
    Double_t px0He = 0.0;
    Double_t py0He = 0.0;
    Double_t pz0He = 0.0;
    Double_t pHe = 0.0;
    Double_t massHe = 0.0;

    Double_t betaC = 0.0;
    Double_t x0C = 0.0;
    Double_t y0C = 0.0;
    Double_t z0C = 0.0;
    Double_t px0C = 0.0;
    Double_t py0C = 0.0;
    Double_t pz0C = 0.0;
    Double_t pC = 0.0;
    Double_t massC = 0.0;

    Double_t betaO = 0.0;
    Double_t x0O = 0.0;
    Double_t y0O = 0.0;
    Double_t z0O = 0.0;
    Double_t px0O = 0.0;
    Double_t py0O = 0.0;
    Double_t pz0O = 0.0;
    Double_t pO = 0.0;
    Double_t massO = 0.0;

    Double_t pBeamz = 17.3915;
    Double_t px0_cand = 0.0;
    Double_t py0_cand = 0.0;
    Double_t pz0_cand = 0.0;
    Double_t beta0_cand = 0.0;
    Double_t pxmem = 0.0;
    Double_t pymem = 0.0;
    Double_t pzmem = 0.0;
    Double_t xmem = 0.0;
    Double_t ymem = 0.0;
    Double_t zmem = 0.0;
    Double_t yC = -1000.0;
    Double_t yCexp = -1000.0;
    Double_t yCfib = -1000.0;
    Double_t yCfibexp = -1000.0;
    Double_t masmem = -1000.0;
    Double_t chargemem = -1000.0;
    Double_t AoverZmem = -1000.0;
    Double_t pCx = -1000.0, pCy = -1000.0, pCz = -10000.0;
    Double_t pHex = -1000.0, pHey = -1000.0, pHez = -10000.0;
    Int_t ixtmax;
    Double_t Erel = 0.0;
    Double_t psum = 0.0;
    Double_t pChi2 = 1e10;
    Double_t eChi2 = 1e10;

    Double_t amu;
    if (!fSimu) // these are atima values in GeV/c2
    {
        amu = 0.93149403;
    }
    else // these are geant3/4 values in GeV/c2
    {
        amu = 0.931494028;
    }

    Double_t AoverZ = 2.0;
    TVector3 pos23a;
    TVector3 pos23b;
    TVector3 postofd;

    // try to fit all possible combination of hits.

    fPropagator->SetVis(kFALSE);

    Int_t nCand = 0;

    Double_t m0Chi2 = -1.0;
    Double_t m1Chi2 = -1.0;
    Double_t yChi2 = -1.0;
    Double_t zChi2 = -1.0;
    minChi2_12C = 1.e6;
    minChi2_4He = 1.e6;

    Int_t charge_mem = 0;
    Bool_t alpha = kFALSE;
    Bool_t carbon = kFALSE;
    Bool_t oxygen = kFALSE;

    // The idea is to loop twice over the ToF wall hits.
    // First we want to look for 12C particle and then for 4He

    Int_t Icountleft = 0;
    Int_t Icountright = 0;

    Int_t lmin;
    Int_t lmax;
    if (fBfield != -1710.0)
    {
        lmin = 0;
        lmax = 1;
    }
    else
    {
        lmin = 1;
        lmax = 3;
    }
    // cout << "Test l: " << lmin << "  " << lmax << "  " << fBfield << endl;

    Int_t iretrack_max = 1;
    Double_t psum_mem = -10000.;

    for (Int_t iretrack = 0; iretrack < iretrack_max + 1; iretrack++)
    {
        Int_t ifi23a = 0;
        Int_t ifi23b = 0;
        Int_t ifi30 = 0;
        Int_t ifi31 = 0;
        Int_t ifi32 = 0;
        Int_t ifi33 = 0;
        Int_t itof = 0;
        if (0 == fi23a->hits.size())
        {
            ifi23a = -1;
        }
        if (0 == fi23b->hits.size())
        {
            ifi23b = -1;
        }
        if (0 == fi30->hits.size())
        {
            ifi30 = -1;
        }
        if (0 == fi31->hits.size())
        {
            ifi31 = -1;
        }
        if (0 == fi32->hits.size())
        {
            ifi32 = -1;
        }
        if (0 == fi33->hits.size())
        {
            ifi33 = -1;
        }
        if (0 == tof->hits.size())
        {
            itof = -1;
        }

        alpha = kFALSE;
        carbon = kFALSE;

        for (Int_t i = 0; i < fi23a->hits.size(); i++)
        {
            fi23a->free_hit[i] = true;
        }
        for (Int_t i = 0; i < fi23b->hits.size(); i++)
        {
            fi23b->free_hit[i] = true;
        }
        for (Int_t i = 0; i < fi30->hits.size(); i++)
        {
            fi30->free_hit[i] = true;
        }
        for (Int_t i = 0; i < fi32->hits.size(); i++)
        {
            fi32->free_hit[i] = true;
        }
        for (Int_t i = 0; i < fi31->hits.size(); i++)
        {
            fi31->free_hit[i] = true;
        }
        for (Int_t i = 0; i < fi33->hits.size(); i++)
        {
            fi33->free_hit[i] = true;
        }

        for (Int_t i = 0; i < tof->hits.size(); i++)
        {
            tof->free_hit[i] = true;
        }

        for (Int_t l = lmin; l < lmax; l++)
        {
            // l = 0: 16O
            // l = 1: 12C
            // l = 1: 4He

            Int_t charge = 0;
            Int_t charge_requested = 0;
            Double_t Charge = 0.0;
            Double_t m0 = 0.0;
            Double_t p0 = 0.0;

            if (l == 0)
            {
                charge_requested = 8;
            }
            else if (l == 1)
            {
                charge_requested = 6;
            }
            else if (l == 2)
            {
                charge_requested = 2;
            }

            if (iretrack == 1 && psum_mem < 0.)
                continue;

            if (iretrack == 0)
            {
                x0 = 0.0;
            }
            if (iretrack == 1 && l < 2 && psum_mem > 0.)
            {
                x0 = xmem - (psum_mem - ps) / 978.518;
            }
            if (l == 2)
            {
                x0 = xmem;
                y0 = ymem;
            }

            z0 = 0.0;

            // Loop over all combinations of hits
            for (Int_t i = 0; i < tof->hits.size(); i++) // loop over all ToFD hits
            {
                if (fSimu)
                {
                    // For tracking of simulations:
                    //   charge = sqrt(tof->hits.at(i)->GetEloss()) * 26.76 + 0.5;
                    //   Charge = sqrt(tof->hits.at(i)->GetEloss()) * 26.76;

                    // if digi used
                    // Charge = tof->hits.at(i)->GetEloss();
                    // charge = (int)(Charge + 0.5);

                    // if digiHit used
                    charge = (int)(tof->hits.at(i)->GetEloss() + 0.5);
                    Charge = tof->hits.at(i)->GetEloss();
                }
                else
                {
                    // For tracking of exp. data:
                    charge = (int)(tof->hits.at(i)->GetEloss() + 0.5);
                    Charge = tof->hits.at(i)->GetEloss();
                }
                if (debug)
                    cout << "Charge: " << charge << " requested charge: " << charge_requested << endl;

                if (charge != charge_requested)
                    continue;

                beta0 = 0.7593; // velocity could eventually be calculated from ToF
                tof->res_t = 0.03;
                // Double_t m0 = charge * 2. * 0.931494028; // First guess of mass

                // Masse wird nicht getreckt, nur momentum
                if (charge == 8)
                {
                    m0 = 14.895085;
                    if (fSimu)
                        m0 = mO;  // has to have the same value as what geant uses
                    p0 = 17.3915; // in GeV/c2
                }
                if (charge == 6)
                {
                    m0 = 11.174862;
                    if (fSimu)
                        m0 = mC;
                    p0 = 13.043625;
                    if (!fSimu)
                        pC = p0;
                    massC = m0 * 1.e3;
                }
                if (charge == 2)
                {
                    m0 = 3.7273791;
                    if (fSimu)
                        m0 = mHe;
                    p0 = 4.347875;
                    if (!fSimu)
                        pHe = p0;
                    massHe = m0 * 1.e3;
                    if (fSimu)
                    {
                        //  x0 = x0He;
                        // y0 = y0He;
                    }
                }

                if (debug)
                {
                    cout << "Mass: " << m0 << endl;
                    cout << "Position on TofD: " << tof->hits.at(i)->GetX() << endl;
                }

                if (!tof->free_hit[i]) // if the hit was used already, continue
                {
                    if (debug)
                        cout << "ToFD hit already used" << endl;
                    continue;
                }

                tof->LocalToGlobal(postofd, tof->hits.at(i)->GetX(), tof->hits.at(i)->GetY());
                Double_t ltofd = sqrt((postofd.Z() - z_tp) * (postofd.Z() - z_tp) + postofd.X() * postofd.X());

                Double_t foffset = 0.441;
                Double_t fslope = 1.0;
                if (tof->hits.at(i)->GetX() > 0.)
                    fslope = 0.883; // fi30/32, tofdx_loc>0
                if (tof->hits.at(i)->GetX() < 0.)
                    fslope = 0.901; // fi31/33, tofdx_loc<0

                if (tof->hits.at(i)->GetX() > 0.0 && tof->hits.at(i)->GetX() != 4.11 && fi30->hits.size() > 0 &&
                    fi32->hits.size() > 0 && fi23a->hits.size() > 0 && fi23b->hits.size() > 0)
                {
                    // left branch in beam direction, don't consider hits in the detectors of the other side

                    if (!fSimu)
                    {
                        //( (R3BGladFieldMap*) FairRunAna::Instance()->GetField())->SetTrackerCorrectionScale(1.002);
                    }
                    // simu from 01/02/22:
                    if (fSimu)
                    {
                        /*   ((R3BGladFieldMap*)FairRunAna::Instance()->GetField())->SetTrackerCorrectionAngleY(-14.);
                           ((R3BGladFieldMap*)FairRunAna::Instance()->GetField())->SetTrackerCorrectionZ(174.95);
                           ((R3BGladFieldMap*)FairRunAna::Instance()->GetField())->SetTrackerCorrectionY(1.08);
                           ((R3BGladFieldMap*)FairRunAna::Instance()->GetField())->SetTrackerCorrectionAngleZ(0.);
                           ((R3BGladFieldMap*)FairRunAna::Instance()->GetField())->SetTrackerCorrectionScale(1.00);*/
                    }

                    // R3BTrackingDetector* target = fDetectorsLeft->GetByName("target");
                    if (fNEventsLeft == 0)
                    {
                        target->hits.push_back(new R3BHit(0, 0.0, 0.0, 0., 0., 0));
                    }

                    do // fi32
                    {

                        if ((ifi32 >= 0 && !fi32->free_hit[ifi32]) ||
                            charge_requested !=
                                fi32->hits.at(ifi32)->GetEloss()) // if the hit was used already, continue
                        {

                            ifi32 += 1;
                            continue;
                        }

                        do // fi30
                        {
                            if ((ifi30 >= 0 && !fi30->free_hit[ifi30]) ||
                                charge_requested !=
                                    fi30->hits.at(ifi30)->GetEloss()) // if the hit was used already, continue
                            {
                                ifi30 += 1;
                                continue;
                            }
                            do // fi23b
                            {
                                if ((l > 0 && abs(fi23b->hits.at(ifi23b)->GetY()) < cut_yfib23) ||
                                    (l < 3 && ifi23b >= 0 &&
                                     !fi23b->free_hit[ifi23b])) // if the hit was used already, continue
                                {

                                    ifi23b += 1;
                                    continue;
                                }
                                fi23b->LocalToGlobal(pos23b, 0., fi23b->hits.at(ifi23b)->GetY());
                                Double_t y_tp =
                                    (postofd.Y() - foffset + fslope * ltofd * pos23b.Y() / (z_tp - pos23b.Z())) /
                                    (1. + fslope * ltofd / (z_tp - pos23b.Z()));
                                y0 = pos23b.Y() - pos23b.Z() * (y_tp - pos23b.Y()) / (z_tp - pos23b.Z());

                                // reject fib23b hits that don't correspond to tofdy:
                                /*  if(abs(y0) > 1.5)
                                  {
                                      ifi23b += 1;
                                      continue;
                                  }*/
                                do // fi23a
                                {

                                    if ((l > 0 && abs(fi23a->hits.at(ifi23a)->GetX()) < cut_xfib23) ||
                                        (l < 3 && ifi23a >= 0 && !fi23a->free_hit[ifi23a]))
                                    {

                                        ifi23a += 1;
                                        continue;
                                    }

                                    fi23a->LocalToGlobal(pos23a, fi23a->hits.at(ifi23a)->GetX(), 0.);
                                    if (l < 2)
                                    {
                                        /*
                                         py0_cand =
                                             (pos23b.Y() - y0 + 0.000365) / pos23b.Z() * p0;
                                         px0_cand =
                                             (pos23a.X() - x0 - 0.0093 * 1.) / pos23a.Z() * p0;
                                         pz0_cand = sqrt(p0*p0-py0_cand*py0_cand-px0_cand*px0_cand);
                                         */
                                        px0_cand = 0.0;
                                        py0_cand = 0.0;
                                        pz0_cand = p0;
                                        Double_t ptot_cand =
                                            sqrt(px0_cand * px0_cand + py0_cand * py0_cand + pz0_cand * pz0_cand);
                                        beta0_cand = sqrt(1.0 / (1.0 + (m0 / ptot_cand) * (m0 / ptot_cand)));

                                        candidate = new R3BTrackingParticle(
                                            charge, x0, y0, z0, px0_cand, py0_cand, pz0_cand, beta0_cand, m0);
                                    }
                                    else if (l == 2)
                                    {
                                        px0_cand = -pxmem;
                                        py0_cand = -pymem;
                                        pz0_cand = sqrt(p0 * p0 - pxmem * pxmem - pymem * pymem);
                                        Double_t ptot_cand =
                                            sqrt(px0_cand * px0_cand + py0_cand * py0_cand + pz0_cand * pz0_cand);
                                        beta0_cand = sqrt(1.0 / (1.0 + (m0 / ptot_cand) * (m0 / ptot_cand)));

                                        candidate = new R3BTrackingParticle(
                                            charge, xmem, ymem, zmem, px0_cand, py0_cand, pz0_cand, beta0_cand, m0);
                                    }

                                    if (debug)
                                    {
                                        cout << "left side of setup" << endl;
                                        cout << "Charge requested: " << charge_requested << endl;
                                        cout << "Start values to fit, x0: " << x0 << " y0: " << y0 << " z0: " << z0
                                             << " p0: " << p0 << " beta0: " << beta0 << " m0: " << m0 << endl;
                                        cout << "Hit Tofd # " << i << " x: " << tof->hits.at(i)->GetX()
                                             << " y: " << tof->hits.at(i)->GetY() << endl;
                                    }
                                    if (ifi23a > -1 && debug)
                                        cout << " Fi23a left # " << ifi23a << " x: " << fi23a->hits.at(ifi23a)->GetX()
                                             << endl;
                                    if (ifi23b > -1 && debug)
                                        cout << " left Fi23b # " << ifi23b << " y: " << fi23b->hits.at(ifi23b)->GetY()
                                             << endl;
                                    if (ifi30 > -1 && debug)
                                        cout << " fi30 # " << ifi30 << " x: " << fi30->hits.at(ifi30)->GetX() << endl;
                                    if (ifi32 > -1 && debug)
                                        cout << " fi32 # " << ifi32 << " x: " << fi32->hits.at(ifi32)->GetX() << endl;
                                    // add points through which tracker has to go:
                                    candidate->AddHit("target", 0);
                                    candidate->AddHit("tofd", i);
                                    candidate->AddHit("fi23a", ifi23a);
                                    candidate->AddHit("fi23b", ifi23b);
                                    candidate->AddHit("fi32", ifi32);
                                    candidate->AddHit("fi30", ifi30);

                                    Int_t status = 10;
                                    if (fForward)
                                    {
                                        status = fFitter->FitTrackMomentumForward(candidate, fDetectors);
                                    }
                                    else
                                    {
                                        // status = fFitter->FitTrackBackward2D(candidate, fDetectors);
                                        status = fFitter->FitTrackMomentumBackward(candidate, fDetectors);
                                    }
                                    if (debug)
                                        cout << " Chi: " << candidate->GetChi2() << "  "
                                             << candidate->GetStartMomentum().Mag() << "  "
                                             << 1000. * (candidate->GetStartMomentum().Mag() - p0) *
                                                    (candidate->GetStartMomentum().Mag() - p0)
                                             << endl;
                                    if (debug)
                                        cout << "--------------------------------" << endl;
                                    nCand += 1;
                                    Icountleft += 1;

                                    //     cout << fNEvents<< ", LEFT SIDE: Charge: "<< charge<<", Momentum: " <<
                                    //     candidate->GetMomentum().Mag()<<
                                    //  ", Momentum Z: "<<candidate->GetMomentum().Z() <<
                                    //   ", Momentum X: "<<candidate->GetMomentum().X()<< endl;

                                    if (TMath::IsNaN(candidate->GetMomentum().Z()))
                                    {
                                        delete candidate;
                                        continue;
                                    }

                                    if (10 > status)
                                    {
                                        if (fForward)
                                        {
                                            candidate->Reset();
                                        }
                                        else
                                        {
                                            candidate->SetStartPosition(candidate->GetPosition());
                                            candidate->SetStartMomentum(-1. * candidate->GetMomentum());
                                            candidate->SetStartBeta(beta0);
                                            candidate->UpdateMomentum();

                                            // candidate->SetStartPosition(candidate->GetPosition()); // @target
                                            // candidate->SetStartMomentum(-1. * candidate->GetMomentum());
                                            // candidate->SetStartBeta(0.8328);
                                            // candidate->SetStartBeta(beta0);
                                            // candidate->UpdateMomentum();
                                            candidate->Reset();

                                            // candidate->GetStartPosition().Print();
                                            // candidate->GetStartMomentum().Print();
                                            // cout << "chi2: " << candidate->GetChi2() << endl;
                                            // status = FitFragment(candidate);
                                        }
                                        if (10 > status)
                                        {
                                            // if(candidate->GetChi2() < 3.)
                                            {
                                                fFragments.push_back(candidate);
                                            }
                                        }
                                        else
                                        {
                                            delete candidate;
                                        }
                                    }
                                    else
                                    {
                                        delete candidate;
                                    }

                                    // return;
                                    ifi23a += 1;
                                } while (ifi23a < fi23a->hits.size());
                                ifi23a = 0;
                                if (0 == fi23a->hits.size())
                                    ifi23a = -1;

                                ifi23b += 1;
                            } while (ifi23b < fi23b->hits.size());
                            ifi23b = 0;
                            if (0 == fi23b->hits.size())
                                ifi23b = -1;

                            ifi30 += 1;
                        } while (ifi30 < fi30->hits.size());
                        ifi30 = 0;
                        if (0 == fi30->hits.size())
                            ifi30 = -1;

                        ifi32 += 1;
                    } while (ifi32 < fi32->hits.size());
                    ifi32 = 0;
                    if (0 == fi32->hits.size())
                        ifi32 = -1;
                } // end if left branch

                if (tof->hits.at(i)->GetX() < 0.0 && fi31->hits.size() > 0 && fi33->hits.size() > 0 &&
                    fi23a->hits.size() > 0 && fi23b->hits.size() > 0)
                {
                    // right branch in beam direction, don't consider hits in the detectors of the other side

                    if (!fSimu)
                    {
                        //( (R3BGladFieldMap*) FairRunAna::Instance()->GetField())->SetTrackerCorrectionScale(1.004);
                    }
                    // simu from 01/02/22:
                    if (fSimu)
                    {
                        /*  ((R3BGladFieldMap*)FairRunAna::Instance()->GetField())->SetTrackerCorrectionAngleY(-14.);
                          ((R3BGladFieldMap*)FairRunAna::Instance()->GetField())->SetTrackerCorrectionZ(174.95);
                          ((R3BGladFieldMap*)FairRunAna::Instance()->GetField())->SetTrackerCorrectionY(1.08);
                          ((R3BGladFieldMap*)FairRunAna::Instance()->GetField())->SetTrackerCorrectionAngleZ(0.);
                          ((R3BGladFieldMap*)FairRunAna::Instance()->GetField())->SetTrackerCorrectionScale(1.00);*/
                    }

                    if (fNEventsRight == 0)
                    {
                        // R3BHit(detId, x,  y, eloss,  time, hitId = -1);
                        target->hits.push_back(new R3BHit(0, 0.0, 0.0, 0., 0., 0));
                    }

                    do // fi33
                    {
                        if ((ifi33 >= 0 && !fi33->free_hit[ifi33]) ||
                            charge_requested !=
                                fi33->hits.at(ifi33)->GetEloss()) // if the hit was used already, continue
                        {

                            ifi33 += 1;
                            continue;
                        }
                        do // fi31
                        {
                            if ((ifi31 >= 0 && !fi31->free_hit[ifi31]) ||
                                charge_requested !=
                                    fi31->hits.at(ifi31)->GetEloss()) // if the hit was used already, continue
                            {

                                ifi31 += 1;
                                continue;
                            }
                            do // fi23b
                            {
                                if ((l > 0 && abs(fi23b->hits.at(ifi23b)->GetY()) < cut_yfib23) ||
                                    (l < 3 && ifi23b >= 0 &&
                                     !fi23b->free_hit[ifi23b])) // if the hit was used already, continue
                                {
                                    ifi23b += 1;
                                    continue;
                                }
                                fi23b->LocalToGlobal(pos23b, 0., fi23b->hits.at(ifi23b)->GetY());
                                Double_t y_tp =
                                    (postofd.Y() - foffset + fslope * ltofd * pos23b.Y() / (z_tp - pos23b.Z())) /
                                    (1. + fslope * ltofd / (z_tp - pos23b.Z()));
                                y0 = pos23b.Y() - pos23b.Z() * (y_tp - pos23b.Y()) / (z_tp - pos23b.Z());
                                // y0 = 0.0;
                                // reject fib23b hits that don't correspond to tofdy:
                                /* if(abs(y0) > 1.5)
                                 {
                                     ifi23b += 1;
                                     continue;
                                 }    */
                                do // fi23a
                                {

                                    if ((l > 0 && abs(fi23a->hits.at(ifi23a)->GetX()) < cut_xfib23) ||
                                        (l < 3 && ifi23a >= 0 && !fi23a->free_hit[ifi23a]))
                                    {
                                        ifi23a += 1;
                                        continue;
                                    }

                                    fi23a->LocalToGlobal(pos23a, fi23a->hits.at(ifi23a)->GetX(), 0.);
                                    if (l < 2)
                                    {
                                        /*py0_cand =
                                            (pos23b.Y() - y0 + 0.000365) / pos23b.Z() * p0;
                                        px0_cand =
                                            (pos23a.X() - x0 - 0.0093 * 1.) / pos23a.Z() * p0;
                                        pz0_cand = sqrt(p0*p0-py0_cand*py0_cand-px0_cand*px0_cand);*/

                                        px0_cand = 0.0;
                                        py0_cand = 0.0;
                                        pz0_cand = p0;
                                        Double_t ptot_cand =
                                            sqrt(px0_cand * px0_cand + py0_cand * py0_cand + pz0_cand * pz0_cand);
                                        beta0_cand = sqrt(1.0 / (1.0 + (m0 / ptot_cand) * (m0 / ptot_cand)));

                                        candidate = new R3BTrackingParticle(
                                            charge, x0, y0, z0, px0_cand, py0_cand, pz0_cand, beta0_cand, m0);
                                    }
                                    else if (l == 2)
                                    {
                                        px0_cand = -pxmem;
                                        py0_cand = -pymem;
                                        pz0_cand = sqrt(p0 * p0 - pxmem * pxmem - pymem * pymem);
                                        Double_t ptot_cand =
                                            sqrt(px0_cand * px0_cand + py0_cand * py0_cand + pz0_cand * pz0_cand);
                                        beta0_cand = sqrt(1.0 / (1.0 + (m0 / ptot_cand) * (m0 / ptot_cand)));

                                        candidate = new R3BTrackingParticle(
                                            charge, xmem, ymem, zmem, px0_cand, py0_cand, pz0_cand, beta0_cand, m0);
                                    }

                                    if (debug)
                                    {
                                        cout << "right side of setup" << endl;
                                        cout << "Charge requested: " << charge_requested << endl;
                                        cout << "Start values to fit, x0: " << x0 << " y0: " << y0 << " z0: " << z0
                                             << " p0: " << p0 << " beta0: " << beta0 << " m0: " << m0 << endl;
                                        cout << "Hit Tofd # " << i << " x: " << tof->hits.at(i)->GetX()
                                             << " y: " << tof->hits.at(i)->GetY() << endl;
                                    }
                                    if (ifi23a > -1 && debug)
                                        cout << "Fi23a # " << ifi23a << " x: " << fi23a->hits.at(ifi23a)->GetX()
                                             << endl;
                                    if (ifi23b > -1 && debug)
                                        cout << "Fi23b # " << ifi23b << " y: " << fi23b->hits.at(ifi23b)->GetY()
                                             << endl;
                                    if (ifi33 > -1 && debug)
                                        cout << "Fi33 # " << ifi33 << " x: " << fi33->hits.at(ifi33)->GetX() << endl;
                                    if (ifi31 > -1 && debug)
                                        cout << "Fi31  # " << ifi31 << " x: " << fi31->hits.at(ifi31)->GetX() << endl;

                                    candidate->AddHit("target", 0);
                                    candidate->AddHit("tofd", i);
                                    candidate->AddHit("fi23a", ifi23a);
                                    candidate->AddHit("fi23b", ifi23b);
                                    candidate->AddHit("fi31", ifi31);
                                    candidate->AddHit("fi33", ifi33);

                                    Int_t status = 10;
                                    if (fForward)
                                    {
                                        status = fFitter->FitTrackMomentumForward(candidate, fDetectors);
                                    }
                                    else
                                    {
                                        // status = fFitter->FitTrackBackward2D(candidate, fDetectors);
                                        status = fFitter->FitTrackMomentumBackward(candidate, fDetectors);
                                    }
                                    if (debug)
                                        cout << "Chi: " << candidate->GetChi2() << "  pstart.Mag "
                                             << candidate->GetStartMomentum().Mag() << " dp.Mag "
                                             << 1000. * (candidate->GetStartMomentum().Mag() - p0) *
                                                    (candidate->GetStartMomentum().Mag() - p0)
                                             << endl;
                                    if (debug)
                                        cout << "--------------------------------" << endl;
                                    nCand += 1;
                                    Icountright += 1;

                                    //    cout <<fNEvents<<", RIGHT SIDE: Charge: "<< charge<<", Momentum: " <<
                                    //    candidate->GetMomentum().Mag()<<", Momentum Z:
                                    //    "<<candidate->GetMomentum().Z()<<
                                    //   ", Momentum X: "<<candidate->GetMomentum().X() << " status: "<<status<<endl;

                                    if (TMath::IsNaN(candidate->GetMomentum().Z()))
                                    {
                                        delete candidate;
                                        continue;
                                    }

                                    if (10 > status)
                                    {
                                        if (fForward)
                                        {
                                            candidate->Reset();
                                        }
                                        else
                                        {
                                            candidate->SetStartPosition(candidate->GetPosition());
                                            candidate->SetStartMomentum(-1. * candidate->GetMomentum());
                                            candidate->SetStartBeta(beta0);
                                            candidate->UpdateMomentum();

                                            // candidate->SetStartPosition(candidate->GetPosition());
                                            // candidate->SetStartMomentum(-1. * candidate->GetMomentum());
                                            // candidate->SetStartBeta(0.8328);
                                            // candidate->SetStartBeta(beta0);
                                            // candidate->UpdateMomentum();

                                            candidate->Reset();

                                            // candidate->GetStartPosition().Print();
                                            // candidate->GetStartMomentum().Print();
                                            // cout << "chi2: " << candidate->GetChi2() << endl;
                                            // status = FitFragment(candidate);
                                        }
                                        if (10 > status)
                                        {
                                            // if(candidate->GetChi2() < 3.)
                                            {
                                                fFragments.push_back(candidate);
                                            }
                                        }
                                        else
                                        {
                                            delete candidate;
                                        }
                                    }
                                    else
                                    {
                                        delete candidate;
                                    }

                                    // return;
                                    ifi23a += 1;
                                } while (ifi23a < fi23a->hits.size());
                                ifi23a = 0;
                                if (0 == fi23a->hits.size())
                                    ifi23a = -1;

                                ifi23b += 1;
                            } while (ifi23b < fi23b->hits.size());
                            ifi23b = 0;
                            if (0 == fi23b->hits.size())
                                ifi23b = -1;

                            ifi31 += 1;
                        } while (ifi31 < fi31->hits.size());
                        ifi31 = 0;
                        if (0 == fi31->hits.size())
                            ifi31 = -1;

                        ifi33 += 1;
                    } while (ifi33 < fi33->hits.size());
                    ifi33 = 0;
                    if (0 == fi33->hits.size())
                        ifi33 = -1;
                } // end if right branch

                charge_mem = charge;

            } // end for TofD

            charge = charge_mem;
            Charge = double(charge_mem);

            // if (candidate->GetSize() > 0 && !fSimu)
            {
                // candidate->Clear();
            }

            //   if(candidate) delete candidate;

            R3BTrackingParticle* bestcandidate;
            //	= new R3BTrackingParticle(-1000, -1000, -1000, -1000, -1000, -1000, -1000, -1000, -1000);

            minChi2 = 1e10;
            pChi2 = 1e10;
            eChi2 = 1e10;
            Double_t parChi2 = 1e10;
            Double_t xChi2 = 1e10;
            Double_t sigmap = 0.01 * p0;

            if (fFragments.size() > 0)
            {
                for (auto const& x : fFragments)
                {
                    // pChi2 = (x->GetStartMomentum().Mag() - p0) * (x->GetStartMomentum().Mag() - p0) /
                    //       (p0 * p0 * 0.01 * 0.01);

                    pChi2 = (x->GetStartMomentum().Mag() - p0) * (x->GetStartMomentum().Mag() - p0) / (sigmap * sigmap);

                    xChi2 = x->GetChi2();

                    // cout<<"For l = "<<l<<" xChi2 = "<<xChi2 <<", pCHi2 = "<<pChi2<<", p =
                    // "<<x->GetStartMomentum().Mag()<<endl;
                    // parChi2 = xChi2;
                    // parChi2 = sqrt(pChi2 * pChi2 + xChi2 * xChi2);

                    if (l < 2) // 12C or 16O
                    {
                        parChi2 = xChi2;
                        // parChi2 = sqrt(pChi2 * pChi2 + xChi2 * xChi2);
                        if (parChi2 < minChi2)
                        {
                            bestcandidate = x;
                            minChi2 = parChi2;

                            // cout << "New min chi2 for 12C: " << minChi2 << endl;
                            // cout << "Corresponding Mass   : " << x->GetMass() << endl;
                            // cout << "Corresponding Mass   : " << bestcandidate->GetMass() << endl;
                        }
                    }
                    if (l == 2) // 4He
                    {
                        // parChi2 = xChi2;
                        // parChi2 = sqrt(xChi2*xChi2 + eChi2*eChi2 + pChi2 * pChi2);
                        pHex = x->GetStartMomentum().X() * 1000.0;
                        pHey = x->GetStartMomentum().Y() * 1000.0;
                        pHez = x->GetStartMomentum().Z() * 1000.0;
                        alphaP.SetPxPyPzE(
                            pHex, pHey, pHez, sqrt(pow(pHex, 2) + pow(pHey, 2) + pow(pHez, 2) + pow(massHe, 2)));
                        p4He = alphaP.Vect();
                        Double_t m_inva = (alphaP + carbonP).M(); // invariant mass
                        Erel = m_inva - massHe - massC;           // relative Energy
                        psum = (p12C + p4He).Mag();
                        pChi2 = (psum - ps) * (psum - ps) / (ps * ps * 0.01 * 0.01);
                        Double_t sigmaErel = 0.01 * 4.36;
                        eChi2 = (Erel - 4.36) * (Erel - 4.36) / (sigmaErel * sigmaErel);

                        // parChi2 = sqrt(xChi2 * xChi2 + pChi2 * pChi2 + eChi2 * eChi2);

                        if (iretrack == 0)
                            parChi2 = xChi2;
                        if (iretrack == iretrack_max)
                            parChi2 = sqrt(xChi2 * xChi2 + eChi2 * eChi2);

                        // cout<<"New 4He candidate: "<<iretrack<<", "<<psum<<", "<<Erel<<", "<<eChi2<<", "<<xChi2<<",
                        // "<<parChi2<<endl; 	cout<<"selected: "<<bestcandidate->GetStartMomentum().Z()<<",
                        //"<<bestcandidate->GetStartPosition().X()<<", "<< 	bestcandidate->GetStartPosition().Y()<<",
                        //"<<minChi2<<", "<<Erel<<endl;

                        if (parChi2 < minChi2)
                        {
                            bestcandidate = x;
                            minChi2 = parChi2;

                            if (l == 2)
                                psum_mem = psum;

                            // cout << "For event: "<<fNEvents_nonull<<" new min chi2 for He: " << minChi2 << ", Erel:
                            // "<<Erel<<", "<<Erel_check<<endl;
                        }
                        //	cout<<"selected: "<<bestcandidate->GetStartMomentum().Z()<<",
                        //"<<bestcandidate->GetStartPosition().X()<<", "<< 	bestcandidate->GetStartPosition().Y()<<",
                        //"<<minChi2<<", "<<Erel<<endl;
                    }
                }

                if (l == 1)
                {
                    pCx = bestcandidate->GetStartMomentum().X() * 1000.0;
                    pCy = bestcandidate->GetStartMomentum().Y() * 1000.0;
                    pCz = bestcandidate->GetStartMomentum().Z() * 1000.0;
                    carbonP.SetPxPyPzE(pCx, pCy, pCz, sqrt(pow(pCx, 2) + pow(pCy, 2) + pow(pCz, 2) + pow(massC, 2)));
                    p12C = carbonP.Vect();
                    // cout << "For event: "<<fNEvents_nonull<<" new min chi2 for C: " << minChi2 <<endl;
                    minChi2_12C = minChi2;
                    pxmem = bestcandidate->GetStartMomentum().X();
                    pymem = bestcandidate->GetStartMomentum().Y();
                    pzmem = bestcandidate->GetStartMomentum().Z();
                    xmem = bestcandidate->GetStartPosition().X();
                    ymem = bestcandidate->GetStartPosition().Y();
                    zmem = bestcandidate->GetStartPosition().Z();
                }
                if (l == 2) // 4He
                {
                    minChi2_4He = minChi2;
                }

                // if (minChi2 > 1.e5) continue;

                if (bestcandidate->GetStartMomentum().X() < 0)
                {
                    fi23a->free_hit[bestcandidate->GetHitIndexByName("fi23a")] = false;
                    fi23b->free_hit[bestcandidate->GetHitIndexByName("fi23b")] = false;
                    fi31->free_hit[bestcandidate->GetHitIndexByName("fi31")] = false;
                    fi33->free_hit[bestcandidate->GetHitIndexByName("fi33")] = false;
                }
                else
                {
                    fi23a->free_hit[bestcandidate->GetHitIndexByName("fi23a")] = false;
                    fi23b->free_hit[bestcandidate->GetHitIndexByName("fi23b")] = false;
                    fi30->free_hit[bestcandidate->GetHitIndexByName("fi30")] = false;
                    fi32->free_hit[bestcandidate->GetHitIndexByName("fi32")] = false;
                }
                tof->free_hit[bestcandidate->GetHitIndexByName("tofd")] = false;

                Double_t x0soll = 0.0;
                Double_t y0soll = 0.0;
                Double_t z0soll = 0.0;
                Double_t psoll = 0.0;
                Double_t px0soll = 0.0;
                Double_t py0soll = 0.0;
                Double_t pz0soll = 0.0;
                Double_t beta0soll = 0.0;
                Double_t m0soll = 0.0;

                if (l == 0)
                {
                    if (debug)
                        cout << "16O" << endl;

                    oxygen = kTRUE;
                    if (fSimu)
                    {
                        x0soll = x0O;
                        y0soll = y0O;
                        z0soll = z0O;
                        px0soll = px0O;
                        py0soll = py0O;
                        pz0soll = pz0O;
                        psoll = pO;
                        m0soll = massO;
                        beta0soll = betaO;
                    }
                    else
                    {
                        x0soll = x0;
                        y0soll = y0;
                        z0soll = z0;
                        px0soll = 0.0;
                        py0soll = 0.0;
                        pz0soll = 17.3915;
                        psoll = 17.3915;
                        m0soll = m0;
                        beta0soll = 0.7593209;
                    }
                }

                if (l == 1)
                {
                    if (debug)
                        cout << "12C" << endl;

                    carbon = kTRUE;
                    x0soll = x0C;
                    y0soll = y0C;
                    z0soll = z0C;
                    px0soll = px0C;
                    py0soll = py0C;
                    pz0soll = pz0C;
                    psoll = pC;
                    m0soll = massC * 1.e-3;
                    beta0soll = betaC;
                    counterC += 1;
                    chargemem = Charge;
                }

                if (l == 2)
                {
                    if (debug)
                        cout << "4He" << endl;
                    alpha = kTRUE;
                    x0soll = xmem; // x0He;
                    y0soll = ymem; // y0He;
                    z0soll = zmem; // z0He;
                    px0soll = px0He;
                    py0soll = py0He;
                    pz0soll = pz0He;
                    psoll = pHe;
                    m0soll = massHe * 1.e-3;
                    beta0soll = betaHe;
                    counterHe += 1;
                }
                if (debug)
                {
                    cout << "Results after tracking :" << endl;
                    cout << "Charge   : " << charge << endl;
                    cout << "Position (soll) x: " << x0soll << " y: " << y0soll << " z: " << z0soll << endl;
                    cout << "Position (ist)  x: " << bestcandidate->GetStartPosition().X()
                         << " y: " << bestcandidate->GetStartPosition().Y()
                         << " z: " << bestcandidate->GetStartPosition().Z() << endl;

                    cout << "Momentum (soll): " << psoll << " px : " << px0soll << " py: " << py0soll
                         << " pz: " << pz0soll << endl;
                    cout << "Momentum (ist) : " << bestcandidate->GetStartMomentum().Mag()
                         << " px : " << bestcandidate->GetStartMomentum().X()
                         << " py: " << bestcandidate->GetStartMomentum().Y()
                         << " pz: " << bestcandidate->GetStartMomentum().Z() << endl;
                    cout << "chi2 " << minChi2 << endl;

                    //   cout << "Beta   : " << bestcandidate->GetStartBeta() << endl;
                }

                if (minChi2 < 1.e6 && l == 2 && (carbon && alpha))
                {
                    totalChi2Mass += (minChi2 + minChi2_12C);
                    totalEvents++;
                }

                if (l == 2)
                    totalChi2P += sqrt(minChi2_12C * minChi2_12C + minChi2 * minChi2);

                if (l == 2 && (carbon && alpha))
                {
                    if (iretrack == 0)
                        m0Chi2 = sqrt(minChi2_12C * minChi2_12C + minChi2_4He * minChi2_4He);
                    if (iretrack == iretrack_max)
                        m1Chi2 = sqrt(minChi2_12C * minChi2_12C + minChi2_4He * minChi2_4He);
                }

                fPropagator->SetVis(fVis);
                bestcandidate->Reset();
            }
            // delete all stored fragments
            fArrayFragments->Clear();
            if (fFragments.size() > 0)
            {
                for (auto const& x : fFragments)
                {
                    delete x;
                }
                fFragments.clear();
            }

        } // end for two particle (12C and 4He)
        /*
       if ((alpha && carbon ))
       {
           fh_Erel->Fill(Erel);
           fh_psum->Fill(psum);
           fh_chi2->Fill(minChi2_12C);
           fh_chi2->Fill(minChi2_4He);
       }
     */

        if (debug)
        {
            cout << "For iretrack = " << iretrack << endl;
            cout << "chi2 12C: " << minChi2_12C << endl;
            cout << "Momentum 12C: " << pCx << "  " << pCy << "  " << pCz << endl;
            cout << "chi2 4He: " << minChi2 << endl;
            cout << "Momentum 4He: " << pHex << "  " << pHey << "  " << pHez << endl;
            cout << "Found Tracks with chi2 He/C= " << minChi2 << " / " << minChi2_12C << ", Erel: " << Erel
                 << ", and psum: " << psum << ", x0: " << xmem << ", y0: " << ymem << endl;
        }

    } // end for iretrack

    if (m1Chi2 < 1000000.0)
    {
        return m1Chi2;
    }
    return -1.;
}

void R3BOptimizeGeometryS494::Finish()
{
    // Optimize global Chi2

    Int_t nvariables = 16; // 14

    // Place here Minuit part for
    // minimization of Chi2AllEvents()

    // ROOT::Math::Minimizer* mini = ROOT::Math::Factory::CreateMinimizer("Minuit2", "Minimize");
    // ROOT::Math::Minimizer* mini = ROOT::Math::Factory::CreateMinimizer("Minuit2", "Migrad");
    // ROOT::Math::Minimizer* mini = ROOT::Math::Factory::CreateMinimizer("Minuit", "MigradImproved");
    // ROOT::Math::Minimizer* mini = ROOT::Math::Factory::CreateMinimizer("Minuit2", "Minos");
    ROOT::Math::Minimizer* mini = ROOT::Math::Factory::CreateMinimizer("Minuit2", "Scan");
    // for SCAN: variable is changed in step/10 steps
    //   ROOT::Math::Minimizer* mini = ROOT::Math::Factory::CreateMinimizer("Minuit2", "Simplex");
    // Minimum = ROOT::Math::Factory::CreateMinimizer("Minuit2", "Migrad");
    // Migrad, Simplex, Scan,
    // set tolerance , etc...
    mini->SetMaxFunctionCalls(1); // for Minuit/Minuit2
    mini->SetMaxIterations(100);  // for GSL
    mini->SetTolerance(10);
    mini->SetStrategy(1);
    mini->SetPrintLevel(1);

    // create funciton wrapper for minmizer
    // a IMultiGenFunction type
    ROOT::Math::Functor* fm1 = new ROOT::Math::Functor(&Chi2AllEvents, nvariables);
    mini->SetFunction(*fm1);

    // Double_t step[16] = { 0.01, 0.013, 0.01, 0.5, 0.9, 0.5, 0.5, 0.9, 0.5, 1., 0.1, 1.0, 0.1, 0.1, 0.1, 0.1 };
    // Double_t step[16] = {0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1 };   // scan
    Double_t step[16] = { 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1 }; // scan

    // Double_t step[16] = { 0.01, 0.05, 0.01, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.2, 0.2, 0.2, 0.1, 0.1, 0.1, 0.1 };

    // simu
    // Double_t variable_true[16] = {0.0, 0.0, 91.04, -62.18069, 570.59330, -82.41713, -112.67330, 537.9038, -131.8124,
    // -130.487, 0, 685.4, -13.68626, -193.8265, -14.02967, -193.88699};
    // Double_t variable_mem[16] = {0.05242090, 0.05805132, 91.04061, -62.60678, 571.5501, -83.10044, -112.4452,
    // 538.5381, -132.5500, -130.8368, 2.263945, 685.7456, -13.85543, -193.8466, -13.95070, -193.9429};

    // default
    //   Double_t variable_default[16] =
    //   {0.,0.,91.2,-62.8664,570.7671,-82.8087,-111.1057,537.8183,-131.0650,-129.9365,0,685.4,-14.0,-194.0,-14.0,-194.0};
    //   Double_t variable_mem[16] =
    //   {0.,0.,91.2,-62.8664,570.7671,-82.8087,-111.1057,537.8183,-131.0650,-129.9365,0,685.4,-14.0,-194.0,-14.0,-194.0};
    // 23rd optim
    Double_t variable_default[16] = { 0.007419462, 0.007327982, 91.02052,  -62.59230, 571.7103, -83.07132,
                                      -112.4498,   538.1335,    -132.5532, -129.6732, 0.492008, 686.1278,
                                      -13.89355,   -193.8936,   -13.99355, -193.9936 }; // 23rd

    Double_t variable_mem[16] = { 0.007419462, 0.007327982, 91.02052,  -62.59230, 571.7103, -83.07132,
                                  -112.4498,   538.1335,    -132.5532, -129.6732, 0.492008, 686.1278,
                                  -13.89355,   -193.8936,   -13.99355, -193.9936 }; // 23rd
    // Nov 2024
    // Double_t variable_default[16] = {0.007419462, 0.007327982, 90.972, -62.59230, 571.7103, -83.07132, -112.4498,
    // 538.1335, -132.5532, -129.6732, 1.492008, 686.1278, -13.66355, -193.8936, -13.89355, -193.9936}; // 23rd Double_t
    // variable_mem[16]     = {0.007419462, 0.007327982, 90.972, -62.59230, 571.7103, -83.07132, -112.4498, 538.1335,
    // -132.5532, -129.6732, 1.492008, 686.1278, -13.66355, -193.8936, -13.89355, -193.9936}; // 23rd
    // Nov 2024
    // Double_t variable_default[16] = { 0.007419462, 0.007327982, 90.972,    -62.59230, 571.7103, -83.07132,
    //                                 -112.4498,   538.1335,    -132.5532, -129.6732, 1.492008, 686.1278,
    //                               -13.89355,   -193.8936,   -13.99355, -193.9936 }; // 23rd
    // Double_t variable_mem[16] = { 0.007419462, 0.007327982, 90.972,    -62.59230, 571.7103, -83.07132,
    //                             -112.4498,   538.1335,    -132.5532, -129.6732, 1.492008, 686.1278,
    //                           -13.89355,   -193.8936,   -13.99355, -193.9936 }; // 23rd

    Double_t variable_sigma[16] = { 0.2, 0.2, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1., 0.5, 0.1, 0.1, 0.1, 0.1 };

    Double_t variable[16];

    Int_t nroundsmax = 1; // 4;//000;
    Int_t nstepchi2 = 0;
    Int_t status = 0;

    for (Int_t i = 0; i < nroundsmax; i++)
    {

        // Idea - first calculate chi2 for start values; if i is too large, skip these values and try new set of start
        // values

        // ********************************************************************************************
        gRandom->SetSeed(0);

        //   variable_mem[3] = gRandom->Gaus(variable_default[3], variable_sigma[3]);    // xfi30
        //   variable_mem[4] = gRandom->Gaus(variable_default[4], variable_sigma[4]);    // zfi30
        //  variable_mem[5] = gRandom->Gaus(variable_default[5], variable_sigma[5]);    // xfi32
        //  variable_mem[12] = gRandom->Gaus(variable_default[12], variable_sigma[12]); // thetafi30
        //  variable_mem[14] = gRandom->Gaus(variable_default[14], variable_sigma[14]); // thetafi32

        /*

                          variable_mem[0] = gRandom->Gaus(variable_default[0], variable_sigma[0]) ;		// xfi23
                          variable_mem[1] = gRandom->Gaus(variable_default[1], variable_sigma[1]);		// yfi23
                          variable_mem[2] = gRandom->Gaus(variable_default[2], variable_sigma[2]);		// zfi23
                          variable_mem[3] = gRandom->Gaus(variable_default[3], variable_sigma[3]);		// xfi30
                          variable_mem[4] = gRandom->Gaus(variable_default[4], variable_sigma[4]);		// zfi30
                          variable_mem[5] = gRandom->Gaus(variable_default[5], variable_sigma[5]);		// xfi32
                          variable_mem[6] = gRandom->Gaus(variable_default[6], variable_sigma[6]);		// xfi31
                          variable_mem[7] = gRandom->Gaus(variable_default[7],  variable_sigma[7]);		// zfi31
                          variable_mem[8] = gRandom->Gaus(variable_default[8], variable_sigma[8]);		// xfi33
                          variable_mem[9] = gRandom->Gaus(variable_default[9], variable_sigma[9]);		// xtofd
                          variable_mem[10] = gRandom->Gaus(variable_default[10], variable_sigma[10]);		// ytofd
                          variable_mem[11] = gRandom->Gaus(variable_default[11],variable_sigma[11]);		// ztofd
                          variable_mem[12] = gRandom->Gaus(variable_default[12],variable_sigma[12]);   	// thetafi30
                          variable_mem[13] = gRandom->Gaus(variable_default[13],variable_sigma[13]);  	// thetafi31
                          if(nvariables == 16){
                              variable_mem[14] = gRandom->Gaus(variable_default[14],variable_sigma[14]);	// thetafi32
                              variable_mem[15] = gRandom->Gaus(variable_default[15],variable_sigma[15]);	// thetafi33
                          }

           */

        for (Int_t ivar = 0; ivar < nvariables; ivar++)
        {
            variable[ivar] = variable_mem[ivar];
            // cout<<"TESTING: varmem: "<<variable_mem[ivar]<<", var: "<<variable[ivar]<<endl;
        }
        cout << " " << endl;
        if (nvariables == 14)
        {
            cout << "Start values = {" << variable[0] << ", " << variable[1] << ", " << variable[2] << ", "
                 << variable[3] << ", " << variable[4] << ", " << variable[5] << ", " << variable[6] << ", "
                 << variable[7] << ", " << variable[8] << ", " << variable[9] << ", " << variable[10] << ", "
                 << variable[11] << ", " << variable[12] << ", " << variable[13] << "}; " << endl;
        }
        if (nvariables == 16)
        {
            cout << "Start values = {" << variable[0] << ", " << variable[1] << ", " << variable[2] << ", "
                 << variable[3] << ", " << variable[4] << ", " << variable[5] << ", " << variable[6] << ", "
                 << variable[7] << ", " << variable[8] << ", " << variable[9] << ", " << variable[10] << ", "
                 << variable[11] << ", " << variable[12] << ", " << variable[13] << ", " << variable[14] << ", "
                 << variable[15] << "}; " << endl;
        }

        mini->SetLimitedVariable(0, "xfi23", variable[0], step[0], variable[0] - 2.0, variable[0] + 2.0);
        mini->SetLimitedVariable(1, "yfi23", variable[1], step[1], variable[1] - 2.0, variable[1] + 2.0);
        mini->SetLimitedVariable(2, "zfi23", variable[2], step[2], variable[2] - 2.0, variable[2] + 2.0);
        mini->SetLimitedVariable(3, "x30", variable[3], step[3], variable[3] - 2.5, variable[3] + 2.5);
        mini->SetLimitedVariable(4, "z30", variable[4], step[4], variable[4] - 2.5, variable[4] + 2.5);
        mini->SetLimitedVariable(5, "x32", variable[5], step[5], variable[5] - 2.5, variable[5] + 2.5);
        mini->SetLimitedVariable(6, "x31", variable[6], step[6], variable[6] - 2.5, variable[6] + 2.5);
        mini->SetLimitedVariable(7, "z31", variable[7], step[7], variable[7] - 2.5, variable[7] + 2.5);
        mini->SetLimitedVariable(8, "x33", variable[8], step[8], variable[8] - 2.5, variable[8] + 2.5);
        mini->SetLimitedVariable(9, "xtofd", variable[9], step[9], variable[9] - 2.5, variable[9] + 2.5);
        mini->SetLimitedVariable(10, "ytofd", variable[10], step[10], variable[10] - 2.5, variable[10] + 2.5);
        mini->SetLimitedVariable(11, "ztofd", variable[11], step[11], variable[11] - 2.5, variable[11] + 2.5);
        mini->SetLimitedVariable(12, "a30", variable[12], step[12], variable[12] - 2.5, variable[12] + 2.5);
        mini->SetLimitedVariable(13, "a31", variable[13], step[13], variable[13] - 2.5, variable[13] + 2.5);
        if (nvariables == 16)
        {
            mini->SetLimitedVariable(14, "a32", variable[14], step[14], variable[14] - 2.5, variable[14] + 2.5);
            mini->SetLimitedVariable(15, "a33", variable[15], step[15], variable[15] - 2.5, variable[15] + 2.5);
        }

        mini->FixVariable(0);
        mini->FixVariable(1);
        mini->FixVariable(2);
        mini->FixVariable(3);
        mini->FixVariable(4);
        mini->FixVariable(5);
        mini->FixVariable(6);
        mini->FixVariable(7);
        mini->FixVariable(8);
        mini->FixVariable(9);
        mini->FixVariable(10);
        mini->FixVariable(11);
        mini->FixVariable(12);
        mini->FixVariable(13);
        if (nvariables == 16)
        {
            mini->FixVariable(14);
            mini->FixVariable(15);
        }

        // do the minimization
        mini->Minimize();

        status = mini->Status();
        cout << "Attempt: " << i << endl;
        cout << "Status: " << status << endl;
        cout << "Start chi2: " << mini->MinValue() << endl;
        Double_t chi2_start = mini->MinValue();

        // *************************************************************************************************

        if (chi2_start > 1000.)
        {
            cout << "Start chi2 to large: " << chi2_start << ", new start values will be chosen" << endl;
        }
        else
        {
            cout << "Start chi2 ok: " << chi2_start << ", optimization will start" << endl;

            mini->ReleaseVariable(3);  // xfi30
            mini->ReleaseVariable(4);  // zfi30
            mini->ReleaseVariable(5);  // xfi32
            mini->ReleaseVariable(12); // thetafi30
            mini->ReleaseVariable(14); // thetafi32

            // do the minimization
            mini->Minimize();

            status = mini->Status();
            cout << "Attempt: " << i << endl;
            cout << "Status: " << status << endl;
            cout << "optimized values: " << endl;
            cout.precision(7);

            if (nvariables == 16)
            {
                cout << "{" << mini->X()[0] << ", " << mini->X()[1] << ", " << mini->X()[2] << ", " << mini->X()[3]
                     << ", " << mini->X()[4] << ", " << mini->X()[5] << ", " << mini->X()[6] << ", " << mini->X()[7]
                     << ", " << mini->X()[8] << ", " << mini->X()[9] << ", " << mini->X()[10] << ", " << mini->X()[11]
                     << ", " << mini->X()[12] << ", " << mini->X()[13] << ", " << mini->X()[14] << ", " << mini->X()[15]
                     << "}; // chi2= " << mini->MinValue() << endl;
            }
            else if (nvariables == 14)
            {
                cout << "{" << mini->X()[0] << ", " << mini->X()[1] << ", " << mini->X()[2] << ", " << mini->X()[3]
                     << ", " << mini->X()[4] << ", " << mini->X()[5] << ", " << mini->X()[6] << ", " << mini->X()[7]
                     << ", " << mini->X()[8] << ", " << mini->X()[9] << ", " << mini->X()[10] << ", " << mini->X()[11]
                     << ", " << mini->X()[12] << ", " << mini->X()[13] << "}; // chi2= " << mini->MinValue() << endl;
            }

            Double_t dzeven, dzodd;
            if (!fSimu)
            {
                dzeven = 62.77813349;
                dzodd = 62.77813349;
            }
            else
            {
                // simu
                dzeven = 632.9688 - 570.59330;
                dzodd = 597.7741 - 537.9038;
            }
            if (nvariables == 16)
            {
                cout << "fi23x = " << mini->X()[0] << ", fi23y = " << mini->X()[1] << ", fi23z = " << mini->X()[2]
                     << ", fi30x = " << mini->X()[3] << ", fi30z = " << mini->X()[4] << ", fi31x = " << mini->X()[6]
                     << ", fi31z = " << mini->X()[7] << ", fi32x = " << mini->X()[5]
                     << ", fi32z = " << mini->X()[4] + dzeven << ", fi33x = " << mini->X()[8]
                     << ", fi33z = " << mini->X()[7] + dzodd << ", tofdx = " << mini->X()[9]
                     << ", tofdy = " << mini->X()[10] << ", tofdz = " << mini->X()[11] << ", fi30a = " << mini->X()[12]
                     << ", fi31a = " << mini->X()[13] << ", fi32a = " << mini->X()[14] << ", fi33a = " << mini->X()[15]
                     << endl;
            }
            else if (nvariables == 14)
            {
                cout << "fi23x = " << mini->X()[0] << ", fi23y = " << mini->X()[1] << ", fi23z = " << mini->X()[2]
                     << ", fi30x = " << mini->X()[3] << ", fi30z = " << mini->X()[4] << ", fi31x = " << mini->X()[6]
                     << ", fi31z = " << mini->X()[7] << ", fi32x = " << mini->X()[5]
                     << ", fi32z = " << mini->X()[4] + dzeven << ", fi33x = " << mini->X()[8]
                     << ", fi33z = " << mini->X()[7] + dzodd << ", tofdx = " << mini->X()[9]
                     << ", tofdy = " << mini->X()[10] << ", tofdz = " << mini->X()[11] << ", fi30a = " << mini->X()[12]
                     << ", fi31a = " << mini->X()[13] << ", fi32a = " << mini->X()[12] << ", fi33a = " << mini->X()[13]
                     << endl;
            }
        }

        mini->Clear();
    }

    //   fh_Erel->Write();
    //   fh_psum->Write();
    //   fh_chi2->Write();
}

Bool_t R3BOptimizeGeometryS494::InitPropagator()
{
    FairField* fairField = FairRunAna::Instance()->GetField();
    R3BGladFieldMap* gladField = NULL;
    fFieldPar->printParams();
    if (2 == fFieldPar->GetType())
    {
        gladField = (R3BGladFieldMap*)fairField;

        if (fPropagator)
        {
            delete fPropagator;
        }
        fPropagator = new R3BTPropagator(gladField, fVis);
    }
    else
    {
        LOG(error) << "Unsupported type of field.";
        return kFALSE;
    }
    return kTRUE;
}

ClassImp(R3BOptimizeGeometryS494)
