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

/******************************************************************************
 *****                   Task for fitting a fragment                      *****
 *****     Created Nov 2016 by D. Kresan                                  *****
 *****     Comments to the algorithm implementation from R. Plag          *****
 *****     Implemented ideal fitter, July 2017, D. Kresan                 *****
 ******************************************************************************/

#include "R3BFragmentTracker.h"
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

using namespace std;

#define SPEED_OF_LIGHT 29.9792458 // cm/ns
#define Amu 0.938272
#define Fair_Amu 0.931494028

R3BFragmentTracker::R3BFragmentTracker(const char* name, Bool_t vis, Int_t verbose)
    : FairTask(name, verbose)
    , fFieldPar(NULL)
    , fPropagator(NULL)
    , fArrayMCTracks(NULL)
    , fDetectors(new R3BTrackingSetup())
    , fArrayFragments(new TClonesArray("R3BTrackingParticle"))
    , fNEvents(0)
    , fVis(vis)
    , fFitter(nullptr)
    , fEnergyLoss(kTRUE)
{
    // this is the list of detectors (active areas) we use for tracking
    fDetectors->AddDetector("target", kTarget, "TargetGeoPar");
    fDetectors->AddDetector("psp", kTargetGlad, "PspGeoPar", "PspDigi");
    fDetectors->AddDetector("fi4", kAfterGlad, "fi4GeoPar", "Fi4Hit");
    fDetectors->AddDetector("fi5", kAfterGlad, "fi5GeoPar", "Fi5Hit");
    fDetectors->AddDetector("fi6", kAfterGlad, "fi6GeoPar", "Fi6Hit");
    fDetectors->AddDetector("tofd", kTof, "tofdGeoPar", "TofdHit");
}

R3BFragmentTracker::~R3BFragmentTracker() {}

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
InitStatus R3BFragmentTracker::Init()
{
    FairRootManager* man = FairRootManager::Instance();

    fArrayMCTracks = dynamic_cast<TClonesArray*>(man->GetObject("MCTrack"));
    if (NULL == fArrayMCTracks)
    {
        LOG(error) << "No MC Track array found in input file.";
        return kERROR;
    }

    man->Register("TrackingParticle", "Tracking", fArrayFragments, kTRUE);

    if (!InitPropagator())
    {
        return kERROR;
    }

    fDetectors->Init();

    fh_mult_psp = new TH1F("h_mult_psp", "Multiplicity PSP", 20, -0.5, 19.5);
    fh_mult_fi4 = new TH1F("h_mult_fi4", "Multiplicity Fi4", 20, -0.5, 19.5);
    fh_mult_fi5 = new TH1F("h_mult_fi5", "Multiplicity Fi5", 20, -0.5, 19.5);
    fh_mult_tofd = new TH1F("h_mult_tofd", "Multiplicity TOFd", 20, -0.5, 19.5);
    fh_eloss_psp_mc = new TH1F("h_eloss_psp_mc", "Energy loss PSP (MC truth)", 2000, 0., 200.);
    fh_eloss_psp = new TH1F("h_eloss_psp", "Energy loss PSP", 2000, 0., 200.);
    fh_eloss_fi4_mc = new TH1F("h_eloss_fi4_mc", "Energy loss Fi4 (MC truth)", 2000, 0., 500.);
    fh_eloss_fi4 = new TH1F("h_eloss_fi4", "Energy loss Fi4", 2000, 0., 500.);
    fh_eloss_fi5_mc = new TH1F("h_eloss_fi5_mc", "Energy loss Fi5 (MC truth)", 2000, 0., 500.);
    fh_eloss_fi5 = new TH1F("h_eloss_fi5", "Energy loss Fi5", 2000, 0., 500.);
    fh_eloss_fi6_mc = new TH1F("h_eloss_fi6_mc", "Energy loss Fi6 (MC truth)", 2000, 0., 1000.);
    fh_eloss_fi6 = new TH1F("h_eloss_fi6", "Energy loss Fi6", 2000, 0., 1000.);
    fh_ncand = new TH1F("h_ncand", "Number of candidates", 100, -0.5, 99.5);
    fh_A_reco1 = new TH1F("h_A_reco1", "Reconstructed mass, step 1", 200., 120., 140.);
    fh_A_reco2 = new TH1F("h_A_reco2", "Reconstructed mass, step 2", 200., 120., 140.);
    fh_mom_res = new TH1F("h_mom_res", "Momentum resolution", 500, -0.05, 0.05);
    fh_mass_res = new TH1F("h_mass_res", "Mass resolution", 10, -5.5, 4.5);
    fh_chi2 = new TH1F("h_chi2", "Chi2", 100, 0., 10.);
    fh_vz_res = new TH1F("h_vz_res", "vz - vz_mc", 200, -1., 1.);
    fh_beta_res = new TH1F("h_beta_res", "beta - beta_mc", 200, -0.1, 0.1);

    Double_t ranges[] = { 1., 1., 1., 1., 1., 10. };
    for (Int_t i = 0; i < 6; i++)
    {
        fh_x_res[i] = new TH1F(Form("h_x_res%d", i), Form("x residual %d", i), 200, -ranges[i], ranges[i]);
        fh_x_pull[i] = new TH1F(Form("h_x_pull%d", i), Form("x pull %d", i), 40, -10., 10.);
    }

    fFitter->Init(fPropagator, fEnergyLoss);

    return kSUCCESS;
}

InitStatus R3BFragmentTracker::ReInit()
{
    SetParContainers();

    if (!InitPropagator())
    {
        return kERROR;
    }

    return kSUCCESS;
}

void R3BFragmentTracker::SetParContainers()
{
    fFieldPar = dynamic_cast<R3BFieldPar*>(FairRuntimeDb::instance()->getContainer("R3BFieldPar"));

    fDetectors->SetParContainers();
}

void R3BFragmentTracker::Exec(const Option_t*)
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

    /* this part needs to be adopted to each experiment / setup
     *
     * First thing is to tell each detector which hit we want to use
     * for the track fitting. This means: manually copy the hit we want to
     * use for tracking from the TClonesArray into the single hit structure.
     *
     * Then call fit_fragment() to optimize the track and calculate
     * particle properties.
     */
    fDetectors->CopyHits();

    R3BTrackingDetector* target = fDetectors->GetByName("target");
    R3BTrackingDetector* psp = fDetectors->GetByName("psp");
    R3BTrackingDetector* fi4 = fDetectors->GetByName("fi4");
    R3BTrackingDetector* fi5 = fDetectors->GetByName("fi5");
    R3BTrackingDetector* fi6 = fDetectors->GetByName("fi6");
    R3BTrackingDetector* tof = fDetectors->GetByName("tofd");

    target->hits.push_back(new R3BHit(0, 0., 0., 0., 0., 0));

    // cout << "Hits: " << psp->hits.size() << "  " << fi4->hits.size() << "  "
    // << fi5->hits.size() << "  " << fi6->hits.size() << "  " << tof->hits.size() << endl;

    // fetch start pos, default momentum and charge from the simulation
    // (just for this test!)
    R3BMCTrack* ion = dynamic_cast<R3BMCTrack*>(fArrayMCTracks->At(0));
    // Check if primary
    if (ion->GetMotherId() != -1)
    {
        // error(); !!!
        // exit(0);
        return;
    }

    // Important: Set charge and initial position and momentum of the particle
    Double_t beta = 1. / TMath::Sqrt(1 + TMath::Power(ion->GetMass() / ion->GetP(), 2));
    R3BTrackingParticle* particle = new R3BTrackingParticle(50.,
                                                            ion->GetStartX(),
                                                            ion->GetStartY(),
                                                            ion->GetStartZ(),
                                                            ion->GetPx(),
                                                            ion->GetPy(),
                                                            ion->GetPz(),
                                                            beta,
                                                            ion->GetMass());
    // LOG(info) << "MC mass " << ion->GetMass();
    // LOG(info) << "MC beta " << beta;

    fh_mult_psp->Fill(psp->hits.size());
    fh_mult_fi4->Fill(fi4->hits.size());
    fh_mult_fi5->Fill(fi5->hits.size());
    fh_mult_tofd->Fill(tof->hits.size());

    /* Note:
     * We DO need to know the beta of the particle in order to calculate
     * its energy loss in matter. For that reason only, we need to separate
     * the momentum into mass and beta and we need to vary these two
     * independently. Pity.
     */

    // try to fit all possible combination of hits. (For now just one)
    // for (all psp hits)
    //    for (all tof hits)
    fPropagator->SetVis(kFALSE);

    Int_t nCand = 0;

    Int_t ipsp = 0;
    Int_t ifi4 = 0;
    Int_t ifi5 = 0;
    Int_t ifi6 = 0;
    Int_t itof = 0;
    if (0 == psp->hits.size())
        ipsp = -1;
    if (0 == fi4->hits.size())
        ifi4 = -1;
    if (0 == fi5->hits.size())
        ifi5 = -1;
    if (0 == fi6->hits.size())
        ifi6 = -1;
    if (0 == tof->hits.size())
        itof = -1;

    {
        do
        {
            if (ipsp >= 0)
                fh_eloss_psp_mc->Fill(psp->hits.at(ipsp)->GetEloss()); // MeV

            do
            {
                if (ifi4 >= 0)
                    fh_eloss_fi4_mc->Fill(fi4->hits.at(ifi4)->GetEloss()); // MeV

                do
                {
                    do
                    {
                        do
                        {
                            tof->res_t = 0.03;

                            Double_t velocity0 = 0.8328 + 0.0003;

                            // Create object for particle which will be fitted
                            R3BTrackingParticle* candidate = new R3BTrackingParticle(
                                particle->GetCharge(), 0., 0., 0., 0., 0., 0., velocity0, 132. * Amu);

                            candidate->AddHit("target", 0);
                            if (ipsp >= 0)
                            {
                                if (psp->hits.at(ipsp)->GetEloss() > 30.)
                                    candidate->AddHit("psp", ipsp);
                                else
                                    candidate->AddHit("psp", -1);
                            }
                            else
                            {
                                candidate->AddHit("psp", -1);
                            }
                            candidate->AddHit("fi4", ifi4);
                            candidate->AddHit("fi5", ifi5);
                            candidate->AddHit("fi6", ifi6);
                            candidate->AddHit("tofd", itof);

                            // find momentum
                            // momin is only a first guess
                            Int_t status = fFitter->FitTrackBackward2D(candidate, fDetectors);

                            nCand += 1;

                            if (TMath::IsNaN(candidate->GetMomentum().Z()))
                            {
                                delete candidate;
                                continue;
                            }

                            if (0 == status)
                            {
                                candidate->SetStartPosition(candidate->GetPosition());
                                candidate->SetStartMomentum(-1. * candidate->GetMomentum());
                                // candidate->SetStartBeta(0.8328);
                                candidate->SetStartBeta(velocity0);
                                candidate->UpdateMomentum();
                                candidate->Reset();

                                // candidate->GetStartPosition().Print();
                                // candidate->GetStartMomentum().Print();
                                // cout << "chi2: " << candidate->GetChi2() << endl;
                                // status = FitFragment(candidate);

                                if (0 == status)
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
                            itof += 1;
                        } while (itof < tof->hits.size());
                        ifi6 += 1;
                    } while (ifi6 < fi6->hits.size());
                    ifi5 += 1;
                } while (ifi5 < fi5->hits.size());
                ifi4 += 1;
            } while (ifi4 < fi4->hits.size());
            ipsp += 1;
        } while (ipsp < psp->hits.size());
    }

    fh_ncand->Fill(nCand);

    R3BTrackingParticle* candidate = new R3BTrackingParticle();
    Double_t minChi2 = 1e10;

    if (fFragments.size() > 0)
    {
        for (auto const& x : fFragments)
        {
            if (x->GetChi2() < minChi2)
            {
                candidate = x;
                minChi2 = x->GetChi2();
            }
        }

        fh_A_reco2->Fill(candidate->GetMass() / Fair_Amu);
        candidate->SetMass(TMath::Nint(candidate->GetMass() / Fair_Amu) * Fair_Amu);
        candidate->UpdateMomentum();
        fFitter->FitTrackBeta(candidate, fDetectors);

        Double_t momentum0 = candidate->GetStartMomentum().Mag();
        LOG(debug1);
        LOG(debug1) << "RESULT : " << momentum0;
        LOG(debug1) << "TRUTH  : " << particle->GetMomentum().Mag();
        LOG(debug1) << "Resolution: " << (momentum0 - particle->GetMomentum().Mag()) / particle->GetMomentum().Mag();
        LOG(debug1) << "Mass   : " << candidate->GetMass() / amu;
        LOG(debug1) << "Truth  : " << particle->GetMass() / amu;
        LOG(debug1) << "Mass resolution : " << (candidate->GetMass() - particle->GetMass()) / particle->GetMass();
        LOG(debug1) << "Beta   : " << candidate->GetStartBeta();
        LOG(debug1) << "Truth  : " << particle->GetStartBeta();
        LOG(debug1) << "Beta resolution : "
                    << (candidate->GetStartBeta() - particle->GetStartBeta()) / particle->GetStartBeta();

        fh_mom_res->Fill((momentum0 - particle->GetStartMomentum().Mag()) / particle->GetStartMomentum().Mag());
        fh_chi2->Fill(candidate->GetChi2());
        int imass = TMath::Nint(candidate->GetMass() / amu);
        int imass_mc = TMath::Nint(particle->GetMass() / amu);
        fh_mass_res->Fill(imass - imass_mc);
        fh_vz_res->Fill(candidate->GetStartPosition().X() - particle->GetStartPosition().X());
        fh_beta_res->Fill((candidate->GetStartBeta() - particle->GetStartBeta()) / particle->GetStartBeta());
        // candidate->GetStartPosition().Print();
        // LOG(info) << (candidate->GetMass() / amu);

        Double_t x_l = 0.;
        Double_t y_l = 0.;
        Int_t iDet = 0;
        fPropagator->SetVis();
        for (auto const& det : fDetectors->GetArray())
        {
            if (kTarget != det->section)
            {
                LOG(debug2) << "Propagate to " << det->GetDetectorName();
                fPropagator->PropagateToDetector(candidate, det);
            }

            if (det->GetDetectorName().EqualTo("psp"))
            { // PSP
                Double_t eloss = det->GetEnergyLoss(candidate);
                fh_eloss_psp->Fill(eloss);
            }

            if (det->GetDetectorName().EqualTo("fi4"))
            { // Fi4
                Double_t eloss = det->GetEnergyLoss(candidate);
                fh_eloss_fi4->Fill(eloss);
            }

            if (det->GetDetectorName().EqualTo("fi4"))
            { // Fi4
                Double_t eloss = det->GetEnergyLoss(candidate);
                fh_eloss_fi4->Fill(eloss);
                cout << "Eloss Fi4: " << eloss << endl;
            }

            if (det->GetDetectorName().EqualTo("fi5"))
            { // Fi5
                Double_t eloss = det->GetEnergyLoss(candidate);
                fh_eloss_fi5->Fill(eloss);
                cout << "Eloss Fi5: " << eloss << endl;
            }

            if (det->GetDetectorName().EqualTo("fi6"))
            { // Fi6
                Double_t eloss = det->GetEnergyLoss(candidate);
                fh_eloss_fi6->Fill(eloss);
                cout << "Eloss Fi6: " << eloss << endl;
            }

            if (fEnergyLoss)
            {
                if (kTof != det->section)
                {
                    Double_t weight = 1.;
                    if (kTarget == det->section)
                    {
                        weight = 0.5;
                    }
                    candidate->PassThroughDetector(det, weight);
                }
            }

            // Convert global track coordinates into local on the det plane
            det->GlobalToLocal(candidate->GetPosition(), x_l, y_l);
            Double_t det_hit_x =
                fDetectors
                    ->GetHit(det->GetDetectorName().Data(), candidate->GetHitIndexByName(det->GetDetectorName().Data()))
                    ->GetX();
            fh_x_res[iDet]->Fill(x_l - det_hit_x);
            fh_x_pull[iDet]->Fill((x_l - det_hit_x) / det->res_x);
            iDet++;
        }
        candidate->Reset();

        //        for (auto const& x : fFragments)
        //        {
        //            delete x;
        //        }
        //        fFragments.clear();
    }

    if (0 == (fNEvents % 10))
    {
        LOG(info) << "Event: " << fNEvents;
    }

    fNEvents += 1;

    delete particle;
}

void R3BFragmentTracker::Finish()
{
    fh_mult_psp->Write();
    fh_mult_fi4->Write();
    fh_mult_fi5->Write();
    fh_mult_tofd->Write();
    fh_eloss_psp_mc->Write();
    fh_eloss_psp->Write();
    fh_eloss_fi4_mc->Write();
    fh_eloss_fi4->Write();
    fh_eloss_fi5_mc->Write();
    fh_eloss_fi5->Write();
    fh_eloss_fi6_mc->Write();
    fh_eloss_fi6->Write();
    fh_ncand->Write();
    fh_A_reco1->Write();
    fh_A_reco2->Write();
    fh_mom_res->Write();
    fh_mass_res->Write();
    fh_chi2->Write();
    for (Int_t i = 0; i < 6; i++)
    {
        fh_x_res[i]->Write();
        fh_x_pull[i]->Write();
    }
    fh_vz_res->Write();
    fh_beta_res->Write();

    if (fVis)
    {
        for (auto const& det : fDetectors->GetArray())
        {
            det->Draw();
        }

        TCanvas* c11 = new TCanvas("c11", "", 10, 10, 1000, 1000);
        c11->Divide(2, 2);
        c11->cd(1);
        fh_mult_fi4->Draw();
        c11->cd(2);
        fh_mult_fi5->Draw();
        c11->cd(3);
        fh_mult_psp->Draw();
        c11->cd(4);
        fh_mult_tofd->Draw();
        c11->cd(0);

        new TCanvas("c2", "", 10, 10, 500, 500);
        fh_eloss_psp_mc->SetLineStyle(2);
        fh_eloss_psp_mc->SetLineColor(2);
        fh_eloss_psp_mc->Draw();
        fh_eloss_psp->Draw("same");

        //        TCanvas* c3 = new TCanvas("c3", "", 200, 10, 1000, 500);
        //        c3->Divide(2, 1);
        //        c3->cd(1);
        //        c3->cd(2);
        //        c3->cd(0);

        TCanvas* c5 = new TCanvas("c5", "", 400, 10, 1000, 500);
        c5->Divide(2, 1);
        c5->cd(1);
        fh_mom_res->Draw();
        c5->cd(2);
        fh_mass_res->Draw();
        c5->cd(0);

        new TCanvas("c6", "", 500, 300, 500, 500);
        fh_A_reco2->Draw();
        fh_A_reco2->GetXaxis()->SetRangeUser(126., 132.);
        //        fh_A_reco2->SetLineColor(2);
        //        fh_A_reco2->Draw("same");

        Int_t bin1 = fh_A_reco2->FindFirstBinAbove(fh_A_reco2->GetMaximum() / 2.);
        Int_t bin2 = fh_A_reco2->FindLastBinAbove(fh_A_reco2->GetMaximum() / 2.);
        Double_t fwhm = fh_A_reco2->GetBinCenter(bin2) - fh_A_reco2->GetBinCenter(bin1);

        LOG(info) << fwhm / fh_A_reco2->GetMean();

        new TCanvas("c7", "", 500, 600, 500, 500);
        fh_chi2->Draw();

        TCanvas* c8 = new TCanvas("c8", "", 10, 10, 1500, 500);
        c8->Divide(3, 2);
        for (Int_t i = 0; i < 6; i++)
        {
            c8->cd(i + 1);
            fh_x_res[i]->Draw();
        }
        c8->cd(0);

        TCanvas* c9 = new TCanvas("c9", "", 10, 10, 1500, 500);
        c9->Divide(3, 2);
        for (Int_t i = 0; i < 6; i++)
        {
            c9->cd(i + 1);
            fh_x_pull[i]->Draw();
        }
        c9->cd(0);

        new TCanvas("c10", "", 10, 10, 500, 500);
        fh_vz_res->Draw();

        new TCanvas("c12", "", 600, 10, 500, 500);
        fh_ncand->Draw();

        new TCanvas("c12", "", 10, 10, 500, 500);
        fh_beta_res->Draw();
    }
}

Bool_t R3BFragmentTracker::InitPropagator()
{
    FairField* fairField = FairRunAna::Instance()->GetField();
    R3BGladFieldMap* gladField = NULL;
    fFieldPar->printParams();
    if (2 == fFieldPar->GetType())
    {
        gladField = dynamic_cast<R3BGladFieldMap*>(fairField);

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

ClassImp(R3BFragmentTracker)
