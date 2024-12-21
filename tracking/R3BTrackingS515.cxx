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

// Created on 8/11/2022 by V.Panin

#include "R3BTrackingS515.h"
#include "R3BBunchedFiberHitData.h"
#include "R3BEventHeader.h"
#include "R3BFrsData.h"
#include "R3BLosHitData.h"
#include "R3BMusicHitData.h"
#include "R3BMwpcHitData.h"
#include "R3BTofdHitData.h"

#include "R3BMCTrack.h"
#include "R3BMDFWrapper.h"
#include "R3BTrack.h"

#include "FairLogger.h"
#include "FairRootManager.h"
#include "FairRunAna.h"
#include "FairRuntimeDb.h"

#include "TCanvas.h"
#include "TClonesArray.h"
#include "TCutG.h"
#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TMath.h"
#include "TVector3.h"
#include <TRandom3.h>
#include <TRandomGen.h>

#include "Math/Factory.h"
#include "Math/Functor.h"
#include "Math/GSLMinimizer.h"
#include "Math/Minimizer.h"
#include "Minuit2/Minuit2Minimizer.h"

#include <array>
#include <cstdlib>
#include <ctime>
#include <fstream>
#include <iostream>
#include <sstream>

R3BTrackingS515* gMDFTracker;

R3BTrackingS515::R3BTrackingS515()
    : R3BTrackingS515("TrackingS515", 1)
{
}

R3BTrackingS515::R3BTrackingS515(const char* name, Int_t iVerbose)
    : FairTask(name, iVerbose)
    , fTrigger(-1)
    , fTpat(-1)
    , fNEvents(0)
    , maxevent(0)
    , DoAlignment(false)
    , fTrackItems(new TClonesArray("R3BTrack"))
    , reference_PoQ(0.)
    , GladCurrent(-1)
    , GladReferenceCurrent(-1)
    , FrsCutMeanZ(-1)
    , FrsCutRadiusZ(-1)
    , FrsCutMeanAoZ(-1)
    , FrsCutRadiusAoZ(-1)
    , FrsBrhoMin(-1)
    , FrsBrhoMax(-1)
    , MusicZMin(-1)
    , MusicZMax(-1)
    , FiberEnergyMin(-1)
    , FiberEnergyMax(-1)
{
}

R3BTrackingS515::~R3BTrackingS515()
{
    if (fTrackItems)
        delete fTrackItems;
}

InitStatus R3BTrackingS515::Init()
{
    LOG(info) << "R3BTrackingS515::Init()";
    FairRootManager* mgr = FairRootManager::Instance();
    if (NULL == mgr)
    {
        LOG(fatal) << "FairRootManager not found";
    }

    fHeader = dynamic_cast<R3BEventHeader*>(mgr->GetObject("EventHeader."));
    if (!fHeader)
    {
        LOG(warn) << "R3BTrackingS515::Init() EventHeader. not found";
    }
    // Reading all detector branches
    assert(DET_MAX + 1 == sizeof(fDetectorNames) / sizeof(fDetectorNames[0]));
    LOG(info) << "Reading " << NOF_FIB_DET << " fiber detectors";
    for (int det = 0; det < DET_MAX; det++)
    {
        fDataItems.push_back(dynamic_cast<TClonesArray*>(mgr->GetObject(Form("%s", fDetectorNames[det]))));
        if (NULL == fDataItems.at(det))
        {
            R3BLOG(fatal, Form("\n\n Cannot find tree branch %s \n\n", fDetectorNames[det]));
        }
    }
    // check if all cuts are properly set
    if (GladCurrent < 0 || GladReferenceCurrent < 0 || FrsCutMeanZ < 0 || FrsCutMeanAoZ < 0 || FrsCutRadiusZ < 0 ||
        FrsCutRadiusAoZ < 0 || FrsBrhoMin < 0 || FrsBrhoMax < 0 || MusicZMin < 0 || MusicZMax < 0 ||
        FiberEnergyMin < 0 || FiberEnergyMax < 0)
    {
        R3BLOG(fatal, Form(" Some cuts are not set or negative values are used\n\n"));
    }
    // Initializing all MDF functions
    LOG(info) << "Reading MDF function for FlightPath";
    MDF_FlightPath = new R3BMDFWrapper(MDF_FlightPath_filename.Data());

    LOG(info) << "Reading MDF function for PoQ";
    MDF_PoQ = new R3BMDFWrapper(MDF_PoQ_filename.Data());

    LOG(info) << "Reading MDF function for TX0";
    MDF_TX0 = new R3BMDFWrapper(MDF_TX0_filename.Data());

    LOG(info) << "Reading MDF function for TY0";
    MDF_TY0 = new R3BMDFWrapper(MDF_TY0_filename.Data());

    LOG(info) << "Reading MDF function for TX1";
    MDF_TX1 = new R3BMDFWrapper(MDF_TX1_filename.Data());

    LOG(info) << "Reading MDF function for TY1";
    MDF_TY1 = new R3BMDFWrapper(MDF_TY1_filename.Data());

    // linking to global pointer (needed by alignment)
    gMDFTracker = this;

    // Request storage of R3BTrack data in the output tree
    mgr->Register("MDFTracks", "MDFTracks data", fTrackItems, kTRUE);

    return kSUCCESS;
}

void R3BTrackingS515::Exec(Option_t* option)
{
    if (fNEvents / 1000. == (int)fNEvents / 1000)
        std::cout << "\rEvents: " << fNEvents << " / " << maxevent << " (" << (int)(fNEvents * 100. / maxevent)
                  << " %) " << std::flush;
    fNEvents += 1;
    FairRootManager* mgr = FairRootManager::Instance();
    if (NULL == mgr)
    {
        LOG(fatal) << "FairRootManager not found";
    }

    // Accessing TClonesArrays from every detector
    auto tofd_DataItems = fDataItems.at(DET_TOFD);
    auto mwpc_DataItems = fDataItems.at(MWPC0_HITDATA);
    auto frs_DataItems = fDataItems.at(FRS_DATA);
    auto music_DataItems = fDataItems.at(MUSIC_HITDATA);
    auto los_DataItems = fDataItems.at(LOS_HITDATA);

    // Minimum multiplicity condition
    if (mwpc_DataItems->GetEntriesFast() != 1 || frs_DataItems->GetEntriesFast() < 1 ||
        music_DataItems->GetEntriesFast() < 1 || los_DataItems->GetEntriesFast() < 1 ||
        fDataItems.at(DET_FI10)->GetEntriesFast() == 0 || fDataItems.at(DET_FI11)->GetEntriesFast() == 0 ||
        fDataItems.at(DET_FI12)->GetEntriesFast() == 0 || tofd_DataItems->GetEntriesFast() < 1)
        return;

    // Get hit data from every upstream detector
    // (for now using only one/first hit)
    auto los_data = dynamic_cast<R3BLosHitData*>(los_DataItems->At(0));

    auto frs_data = dynamic_cast<R3BFrsData*>(frs_DataItems->At(0));
    if (!IsInsideFrsEllipticPIDCut(frs_data->GetZ(), frs_data->GetAq()))
        return;

    auto music_hit = dynamic_cast<R3BMusicHitData*>(music_DataItems->At(0));
    if (music_hit->GetZcharge() < MusicZMin || music_hit->GetZcharge() > MusicZMax)
        return;

    auto mwpc_hit = dynamic_cast<R3BMwpcHitData*>(mwpc_DataItems->At(0));

    // Reading fiber detectors (using only highest energy hit)
    R3BBunchedFiberHitData *f10_hit{}, *f11_hit{}, *f12_hit{}, *this_hit{};
    int f{};
    double energy = 0;
    for (f = 0; f < NOF_FIB_DET; ++f)
    {
        energy = 0;
        for (auto i = 0; i < fDataItems.at(f)->GetEntriesFast(); ++i)
        {
            this_hit = dynamic_cast<R3BBunchedFiberHitData*>(fDataItems.at(f)->At(i));
            if (this_hit->GetEloss() > FiberEnergyMin && this_hit->GetEloss() < FiberEnergyMax &&
                this_hit->GetEloss() > energy)
            {
                if (f == DET_FI10)
                {
                    f10_hit = dynamic_cast<R3BBunchedFiberHitData*>(fDataItems.at(f)->At(i));
                }
                else if (f == DET_FI11)
                {
                    f11_hit = dynamic_cast<R3BBunchedFiberHitData*>(fDataItems.at(f)->At(i));
                }
                else if (f == DET_FI12)
                {
                    f12_hit = dynamic_cast<R3BBunchedFiberHitData*>(fDataItems.at(f)->At(i));
                }
                energy = this_hit->GetEloss();
            }
        }
        if (energy == 0)
            break; // no good energy in this fiber
    }
    if (f != 3)
        return; // good hits should be found in all 3 fibers
    assert(f10_hit && f11_hit && f12_hit && "Not all fibers found.");
    // Reading TOFD data
    R3BTofdHitData* tofd_hit;
    bool is_good_tofd = false;
    for (auto i = 0; i < tofd_DataItems->GetEntriesFast(); ++i)
    {
        tofd_hit = dynamic_cast<R3BTofdHitData*>(tofd_DataItems->At(i));
        if (tofd_hit->GetDetId() == 1) // only hits from first plane
            is_good_tofd = true;
    }
    if (!is_good_tofd)
        return;

    // Set TVector3 for every detector hit
    mwpc_point.SetXYZ(mwpc_hit->GetX() * 0.1, mwpc_hit->GetY() * 0.1, 0.); // converting mm to cm
    f10_point.SetXYZ(f10_hit->GetX(), 0, 0);
    f11_point.SetXYZ(f11_hit->GetX(), 0, 0);
    f12_point.SetXYZ(0, f12_hit->GetY(), 0);
    tofd_point.SetXYZ(tofd_hit->GetX(), tofd_hit->GetY(), 0);

    //-------- Collect detector data for the alignment
    if (DoAlignment)
    {
        det_points align_data;
        align_data.f10.SetXYZ(f10_point.X(), f10_point.Y(), f10_point.Z());
        align_data.f11.SetXYZ(f11_point.X(), f11_point.Y(), f11_point.Z());
        align_data.f12.SetXYZ(f12_point.X(), f12_point.Y(), f12_point.Z());
        align_data.mwpc.SetXYZ(mwpc_point.X(), mwpc_point.Y(), mwpc_point.Z());
        det_points_vec.push_back(align_data);
    }
    //------- Transform detector hits to the lab coordinates
    // This will modify det_point vectors
    TransformPoint(mwpc_point, mwpc_angles, mwpc_position);
    TransformPoint(f10_point, f10_angles, f10_position);
    TransformPoint(f11_point, f11_angles, f11_position);
    TransformPoint(f12_point, f12_angles, f12_position);
    TransformPoint(tofd_point, tofd_angles, tofd_position);

    //------ Additionally extrapolate X and Z coordinates of f12
    // Define two (X,Z) points on the f12 plane:
    TVector3 f12_edge[2];
    f12_edge[0].SetXYZ(-1, 0, 0);
    f12_edge[1].SetXYZ(1, 0, 0);
    TransformPoint(f12_edge[0], f12_angles, f12_position);
    TransformPoint(f12_edge[1], f12_angles, f12_position);
    // Parameterize f12 plane
    double f12_slope = (f12_edge[1].X() - f12_edge[0].X()) / (f12_edge[1].Z() - f12_edge[0].Z());
    double f12_offset = f12_edge[0].X() - f12_slope * f12_edge[0].Z();
    // Get track info from f10 and f11 points
    double track_slope = (f11_point.X() - f10_point.X()) / (f11_point.Z() - f10_point.Z());
    double track_offset = (f11_point.X() - track_slope * f11_point.Z());
    // Extrapolate final X and Z position in f12
    f12_point.SetZ((track_offset - f12_offset) / (f12_slope - track_slope)); // extrapolated
    f12_point.SetX(track_slope * f12_point.Z() + track_offset);              // extrapolated

    // Calculate all necessary MDF values
    // preserve the order, it is expected by the MDF function!
    mdf_data[0] = mwpc_point.X();
    mdf_data[1] = mwpc_point.Y();
    mdf_data[2] = mwpc_point.Z();
    mdf_data[3] = f10_point.X();
    mdf_data[4] = f10_point.Z();
    mdf_data[5] = f11_point.X();
    mdf_data[6] = f11_point.Z();
    mdf_data[7] = (f12_point.Y() - mwpc_point.Y()) / (f12_point.Z() - mwpc_point.Z());

    PoQ = MDF_PoQ->MDF(mdf_data) * GladCurrent / GladReferenceCurrent;
    TX0 = MDF_TX0->MDF(mdf_data);
    TX1 = MDF_TX1->MDF(mdf_data);
    TY0 = MDF_TY0->MDF(mdf_data);
    TY1 = MDF_TY1->MDF(mdf_data);

    //----- Calculate additional FlightPath from MWPC to the target
    double targ_dx = mwpc_position.Z() * TX0;
    double targ_dy = mwpc_position.Z() * TY0;
    double fp_mwpc_to_target = sqrt(pow(targ_dx, 2) + pow(targ_dy, 2) + pow(mwpc_position.Z(), 2));

    //----- Calculate additional FlightPath from F11 to TOFD
    TVector3 vec_f11_to_tofd;
    vec_f11_to_tofd = tofd_position - f11_position;
    double fp_f11_to_tofd = vec_f11_to_tofd.Mag();

    // Final flight path from target to tofd
    FlightPath = MDF_FlightPath->MDF(mdf_data) + fp_mwpc_to_target + fp_f11_to_tofd;

    //----- Calculate Beta
    // ToF = tofd_hit->GetTof() + tof_offset;
    ToF = FlightPath / frs_data->GetBeta() / SPEED_OF_LIGHT + tof_offset;

    Beta = FlightPath / ToF / SPEED_OF_LIGHT;
    Gamma = 1. / sqrt(1 - pow(Beta, 2));
    mdf_AoZ = PoQ / Beta / Gamma / AMU;

    // Final momenutm (PoQ) vector
    TVector3 vec_PoQ(TX0, TY0, 1);
    vec_PoQ.SetMag(PoQ);

    AddTrackData(mwpc_point, vec_PoQ, music_hit->GetZcharge(), mdf_AoZ); // chix, chiy, quality
    return;
}

void R3BTrackingS515::FinishEvent()
{
    for (auto& DataItem : fDataItems)
    {
        DataItem->Clear();
    }
    fTrackItems->Clear();
}

void R3BTrackingS515::FinishTask()
{
    LOG(info) << "Processed " << fNEvents << " events\n\n";
    if (DoAlignment)
        Alignment();
}

R3BTrack* R3BTrackingS515::AddTrackData(TVector3 mw, TVector3 poq, Double_t charge, Double_t aoz)
{
    // Filling output track info
    TClonesArray& clref = *fTrackItems;
    Int_t size = clref.GetEntriesFast();
    return new (clref[size]) R3BTrack(mw.X(), mw.Y(), mw.Z(), poq.X(), poq.Y(), poq.Z(), charge, aoz, 0., 0., 0);
}

void R3BTrackingS515::TransformPoint(TVector3& point, TVector3 rot, TVector3 trans)
{
    TRotation r;
    // First Euler rotation around Y axis
    r.RotateY(rot.Y());
    // get local X axis after first rotation
    TVector3 v3_localX(1, 0, 0);
    v3_localX.SetTheta(r.ThetaX());
    v3_localX.SetPhi(r.PhiX());
    // Second Euler rotation around local X axis
    r.Rotate(rot.X(), v3_localX);
    // get local Z axis after second rotation
    TVector3 v3_localZ(0, 0, 1);
    v3_localZ.SetTheta(r.ThetaZ());
    v3_localZ.SetPhi(r.PhiZ());
    // final rotation around local Z axis
    r.Rotate(rot.Z(), v3_localZ);
    // Applying rotation
    point.Transform(r);
    // Applying translation
    point += trans;
    return;
}

void R3BTrackingS515::SetFrsEllipticPIDCut(double meanZ, double rZ, double meanAoZ, double rAoZ)
{
    FrsCutMeanZ = meanZ;
    FrsCutRadiusZ = rZ;
    FrsCutMeanAoZ = meanAoZ;
    FrsCutRadiusAoZ = rAoZ;
}

bool R3BTrackingS515::IsInsideFrsEllipticPIDCut(double z, double aoz)
{
    Double_t ell = pow(((z - FrsCutMeanZ) / (FrsCutRadiusZ)), 2) + pow(((aoz - FrsCutMeanAoZ) / (FrsCutRadiusAoZ)), 2);
    return (ell <= 1.);
}

void R3BTrackingS515::Alignment()
{
    int run_flag = 0;
    std::cout << "\n\n----------------------------------------------------";
    std::cout << "\n Ready for detector alignment using the following data set:";
    std::cout << "\n\tNumber of reference tracks: " << det_points_vec.size();
    std::cout << "\n\tFRS Brho: min = " << FrsBrhoMin << ", max = " << FrsBrhoMax;
    std::cout << "\n\tFRS PoQ: min = " << FrsBrhoMin / 3.3356 << ", max = " << FrsBrhoMax / 3.3356;
    std::cout << "\n\tP/Z reference: " << reference_PoQ << "GeV/c";
    std::cout << "\n\n Do you want to continue? (0 = no, 1 = yes):  ";
    std::cin >> run_flag;
    if (run_flag == 0)
        return;

    // Now define minimizer
    ROOT::Math::Minimizer* minimizer = ROOT::Math::Factory::CreateMinimizer("Minuit2", "Migrad");
    const UInt_t NVarsFunctor = 20; // number of the alignment offsets: (2 angles + 3 offsets) x 4 detectors
    const double* xs;
    double step[NVarsFunctor];
    double offset[NVarsFunctor];
    Double_t min_offset[NVarsFunctor];
    Double_t max_offset[NVarsFunctor];
    TH1F* hist[NVarsFunctor];
    char hname[100];

    // For every detector: 3 Euler angles and 3 shifts
    for (auto d = 0; d < 4; ++d)
    {
        for (auto a = 0; a < 2; ++a) // angle shifts in rad
        {
            min_offset[d * 5 + a] = -0.04;
            max_offset[d * 5 + a] = 0.04;
            step[d * 5 + a] = 0.001;
        }
        for (auto o = 0; o < 3; ++o) // position shifts in cm
        {
            min_offset[d * 5 + 2 + o] = -4.;
            max_offset[d * 5 + 2 + o] = 4.;
            step[d * 5 + 2 + o] = 0.01;
        }
    }
    for (UInt_t i = 0; i < NVarsFunctor; i++)
    {
        sprintf(hname, "par%d", i);
        hist[i] = new TH1F(hname, hname, 100, min_offset[i], max_offset[i]);
    }
    std::cout << "\n\n-- Perfroming minimization for the detector alignment. Please wait... \n";

    // Setting up Minimizer function parameters
    Double_t precision = 1e-10; // 0 - default precision will be automaticalle determined
    Double_t tolerance = 0.02;
    minimizer->SetMaxFunctionCalls(1000000000); // for Minuit/Minuit2
    minimizer->SetMaxIterations(100);           // for GSL
    minimizer->SetTolerance(tolerance);
    minimizer->SetPrecision(precision);
    minimizer->SetPrintLevel(2);
    minimizer->SetStrategy(0); // to run faster
    ROOT::Math::Functor f(&R3BTrackingS515::AlignmentErrorS515, NVarsFunctor);
    minimizer->SetFunction(f);
    minimizer->SetPrintLevel(0);

    // Getting experimental data and doing minimization
    UInt_t i;
    Int_t bs = 0;
    bool good_minimum = false;
    while (bs < 100) // running several minimizations
    {
        minimizer->Clear();
        for (i = 0; i < NVarsFunctor; i++) // sampling +=50% from limits
        {
            offset[i] = gRandom->Uniform(min_offset[i] * 0.5, max_offset[i] * 0.5);
            minimizer->SetVariable(i, Form("par%d", i), offset[i], step[i]);
            minimizer->SetVariableLimits(i, min_offset[i], max_offset[i]);
        }
        minimizer->Minimize();
        xs = minimizer->X();
        // if(minimizer->Status() !=0) continue; //valid minimum
        // Check if all paramters are "far" from limits
        for (i = 0; i < NVarsFunctor; i++)
        {
            if (fabs((xs[i] - min_offset[i]) / min_offset[i]) < 0.1 ||
                fabs((max_offset[i] - xs[i]) / max_offset[i]) < 0.1)
                break;
        }
        if (i != NVarsFunctor)
        {
            // std::cout << "\n\n -- Parameter is too close to the limit!! ";
            continue;
        }
        else
        {
            std::cout << "\n\n-- Good parameters are found!!";
            good_minimum = true;
        }
        std::cout << "\n-- Minimization No. " << bs << "\n-- Minimum: f(";
        for (i = 0; i < NVarsFunctor; i++)
        {
            std::cout << xs[i];
            hist[i]->Fill(xs[i]);
            if (i != NVarsFunctor - 1)
                std::cout << ",";
        }
        std::cout << ") = " << minimizer->MinValue();
        std::cout << "\n-- Minimizer status: " << minimizer->Status() << std::endl;
        bs++;
    }
    // Outputing final info and histograms
    if (minimizer->MinValue() < tolerance && f(xs) < tolerance)
        std::cout << "-- Minimizer "
                  << "   converged to the right minimum" << std::endl;
    else
    {
        std::cout << "-- Minimizer "
                  << "   failed to converge !!!" << std::endl;
    }
    TCanvas* c1 = new TCanvas("c1", "c1", 1000, 1000);
    c1->Divide(5, 4);
    for (UInt_t v = 0; v < NVarsFunctor; v++)
    {
        c1->cd(v + 1);
        hist[v]->Draw();
    }
    return;
}

double R3BTrackingS515::AlignmentErrorS515(const double* par)
{
    gMDFTracker->mwpc_ang_offset.SetXYZ(0, par[0], par[1]);
    gMDFTracker->mwpc_pos_offset.SetXYZ(par[2], par[3], par[4]);

    gMDFTracker->f10_ang_offset.SetXYZ(0, par[5], par[6]);
    gMDFTracker->f10_pos_offset.SetXYZ(par[7], par[8], par[9]);

    gMDFTracker->f11_ang_offset.SetXYZ(0, par[10], par[11]);
    gMDFTracker->f11_pos_offset.SetXYZ(par[12], par[13], par[14]);

    gMDFTracker->f12_ang_offset.SetXYZ(0, par[15], par[16]);
    gMDFTracker->f12_pos_offset.SetXYZ(par[17], par[18], par[19]);

    Double_t mdf_input[8]; // data container for the MDF function

    double v2 = 0;
    double v = 0;
    int counter = 0;
    for (auto& d : (gMDFTracker->det_points_vec))
    {
        gMDFTracker->mwpc_point = d.mwpc;
        gMDFTracker->f10_point = d.f10;
        gMDFTracker->f11_point = d.f11;
        gMDFTracker->f12_point = d.f12;

        // This will transform "det_point" vectors into lab frame
        gMDFTracker->TransformPoint(gMDFTracker->mwpc_point,
                                    gMDFTracker->GetEulerAnglesMWPC() + gMDFTracker->mwpc_ang_offset,
                                    gMDFTracker->GetPositionMWPC() + gMDFTracker->mwpc_pos_offset);

        gMDFTracker->TransformPoint(gMDFTracker->f10_point,
                                    gMDFTracker->GetEulerAnglesF10() + gMDFTracker->f10_ang_offset,
                                    gMDFTracker->GetPositionF10() + gMDFTracker->f10_pos_offset);

        gMDFTracker->TransformPoint(gMDFTracker->f11_point,
                                    gMDFTracker->GetEulerAnglesF11() + gMDFTracker->f11_ang_offset,
                                    gMDFTracker->GetPositionF11() + gMDFTracker->f11_pos_offset);

        gMDFTracker->TransformPoint(gMDFTracker->f12_point,
                                    gMDFTracker->GetEulerAnglesF12() + gMDFTracker->f12_ang_offset,
                                    gMDFTracker->GetPositionF12() + gMDFTracker->f12_pos_offset);

        mdf_input[0] = gMDFTracker->mwpc_point.X();
        mdf_input[1] = gMDFTracker->mwpc_point.Y();
        mdf_input[2] = gMDFTracker->mwpc_point.Z();
        mdf_input[3] = gMDFTracker->f10_point.X();
        mdf_input[4] = gMDFTracker->f10_point.Z();
        mdf_input[5] = gMDFTracker->f11_point.X();
        mdf_input[6] = gMDFTracker->f11_point.Z();
        mdf_input[7] = (gMDFTracker->f12_point.Y() - mdf_input[1]) / (gMDFTracker->f12_point.Z() - mdf_input[2]);

        v2 += pow((gMDFTracker->Get_MDF_PoQ()->MDF(mdf_input) - gMDFTracker->GetReferencePoQ()), 2);
        counter++;
    }
    v2 /= counter;
    v = sqrt(v2);
    // std::cout << "\nReturning error: " << v;
    return v;
}
ClassImp(R3BTrackingS515);
