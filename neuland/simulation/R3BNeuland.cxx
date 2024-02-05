/******************************************************************************
 *   Copyright (C) 2019 GSI Helmholtzzentrum für Schwerionenforschung GmbH    *
 *   Copyright (C) 2019-2024 Members of R3B Collaboration                     *
 *                                                                            *
 *             This software is distributed under the terms of the            *
 *                 GNU General Public Licence (GPL) version 3,                *
 *                    copied verbatim in the file "LICENSE".                  *
 *                                                                            *
 * In applying this license GSI does not waive the privileges and immunities  *
 * granted to it by virtue of its status as an Intergovernmental Organization *
 * or submit itself to any jurisdiction.                                      *
 ******************************************************************************/

#include "R3BNeuland.h"
#include "FairRun.h"
#include "FairRuntimeDb.h"
#include "R3BMCStack.h"
#include "R3BNeulandGeoPar.h"
#include "R3BNeulandPoint.h"
#include "TClonesArray.h"
#include "TGeoBBox.h"
#include "TGeoManager.h"
#include "TParticle.h"
#include "TVirtualMC.h"
#include <FairRootManager.h>
#include <fmt/format.h>

// Initialize variables from Birk' s Law
constexpr auto seconds_to_nanoseconds = 1e9;
constexpr auto BirkdP = 1.032;
constexpr auto BirkC1 = 0.013 / BirkdP;
constexpr auto BirkC2 = 9.6e-6 / (BirkdP * BirkdP);

inline auto GetLightYield(const int charge, const double length, const double edep) -> double
{
    // Apply Birk's law ( Adapted from G3BIRK/Geant3)
    if (charge != 0 && length > 0)
    {
        auto birkC1Mod = BirkC1;

        // Apply correction for higher charge states
        if (TMath::Abs(charge) >= 2)
        {
            birkC1Mod *= 7.2 / 12.6; // NOLINT
        }

        const double dedxcm = 1000. * edep / length;
        const double lightYield = edep / (1. + birkC1Mod * dedxcm + BirkC2 * dedxcm * dedxcm);
        return lightYield;
    }
    return edep; // Rarely very small energy depositions have no length?
}

R3BNeuland::R3BNeuland()
    : R3BNeuland("")
{
}

R3BNeuland::R3BNeuland(const TString& geoFile, const TGeoTranslation& trans, const TGeoRotation& rot)
    : R3BNeuland(geoFile, { trans, rot })
{
}

R3BNeuland::R3BNeuland(const TString& geoFile, const TGeoCombiTrans& combi)
    : R3BDetector("R3BNeuland", kNEULAND, geoFile, combi)
{
}

R3BNeuland::R3BNeuland(int nDP, const TGeoTranslation& trans, const TGeoRotation& rot)
    : R3BNeuland(nDP, { trans, rot })
{
}

R3BNeuland::R3BNeuland(const int nDP, const TGeoCombiTrans& combi)
    : R3BNeuland(fmt::format("neuland_v3_{}dp.geo.root", nDP), combi)
{
}

auto R3BNeuland::ProcessHits(FairVolume* /*v*/) -> bool
{
    // New hit in detector
    if (gMC->IsTrackEntering())
    {
        if (!is_last_hit_done_)
        {
            LOG(warn) << "R3BNeuland: Incomplete hit discarded";
            ResetValues();
        }

        is_last_hit_done_ = kFALSE;
        energy_loss_ = 0.;
        light_yield_ = 0.;
        time_ = gMC->TrackTime() * seconds_to_nanoseconds;
        length_ = gMC->TrackLength();
        gMC->TrackPosition(pos_in_);
        gMC->TrackMomentum(mom_in_);
        gMC->CurrentVolOffID(1, paddle_id_);

        particle_id_ = gMC->TrackPid();
        trackid_pid_map_.emplace(gMC->GetStack()->GetCurrentTrackNumber(), gMC->TrackPid());
        if (auto search = trackid_pid_map_.find(gMC->GetStack()->GetCurrentParentTrackNumber());
            search != trackid_pid_map_.end())
        {
            partent_particle_id_ = search->first;
        }
    }

    // Sum energy loss for all steps in the active volume
    energy_loss_ += gMC->Edep();
    light_yield_ += GetLightYield(gMC->TrackCharge(), gMC->TrackStep(), gMC->Edep());

    // Set additional parameters at exit of active volume. Create R3BNeulandPoint.
    if (gMC->IsTrackExiting() || gMC->IsTrackStop() || gMC->IsTrackDisappeared())
    {
        // Do not save a hit if no energy deposited
        constexpr auto minimum_energy_cutoff = 1e-20;
        if (energy_loss_ < minimum_energy_cutoff || light_yield_ < minimum_energy_cutoff)
        {
            ResetValues();
            return kTRUE;
        }

        track_id_ = gMC->GetStack()->GetCurrentTrackNumber();
        gMC->TrackPosition(pos_out_);
        gMC->TrackMomentum(mom_out_);

        // Add Point
        LOG(debug) << "R3BNeuland: Adding Point at (" << pos_in_.X() << ", " << pos_in_.Y() << ", " << pos_in_.Z()
                   << ") cm,  paddle " << paddle_id_ << ", track " << track_id_ << ", energy loss " << energy_loss_
                   << " GeV " << gMC->GetStack()->GetCurrentParentTrackNumber();

        neuland_points_.get().emplace_back(track_id_,
                                           paddle_id_,
                                           pos_in_.Vect(),
                                           mom_in_.Vect(),
                                           time_,
                                           length_,
                                           energy_loss_,
                                           gMC->CurrentEvent(),
                                           light_yield_,
                                           particle_id_,
                                           partent_particle_id_);

        // Increment number of LandPoints for this track
        auto* stack = dynamic_cast<R3BStack*>(gMC->GetStack());
        stack->AddPoint(kNEULAND);
        ResetValues();
    }

    return kTRUE;
}

auto R3BNeuland::CheckIfSensitive(std::string name) -> bool { return name == "volBC408"; }

void R3BNeuland::EndOfEvent()
{
    if (fVerboseLevel != 0)
    {
        Print();
    }
    Reset();
}

void R3BNeuland::Print(Option_t* /*unused*/) const
{
    LOG(info) << "R3BNeuland: " << neuland_points_.get_constref().size() << " Neuland Points registered in this event";
}

void R3BNeuland::Reset()
{
    neuland_points_.clear();
    ResetValues();
    trackid_pid_map_.clear();
}

void R3BNeuland::ResetValues()
{
    is_last_hit_done_ = kTRUE;
    track_id_ = 0;
    paddle_id_ = -1;
    pos_in_.Clear();
    pos_out_.Clear();
    mom_in_.Clear();
    mom_out_.Clear();
    time_ = length_ = energy_loss_ = light_yield_ = 0;
}

void R3BNeuland::WriteParameterFile()
{
    FairRuntimeDb* rtdb = FairRun::Instance()->GetRuntimeDb();
    neuland_geo_par_ = dynamic_cast<R3BNeulandGeoPar*>(rtdb->getContainer("R3BNeulandGeoPar"));

    // Really bad way to find the Neuland *node* (not the volume!)
    TGeoNode* geoNodeNeuland = nullptr;
    for (int i = 0; i < gGeoManager->GetTopNode()->GetNdaughters(); i++)
    {
        if (TString(gGeoManager->GetTopNode()->GetDaughter(i)->GetVolume()->GetName()) == "volNeuland")
        {
            geoNodeNeuland = gGeoManager->GetTopNode()->GetDaughter(i);
            break;
        }
    }

    if (geoNodeNeuland == nullptr)
    {
        LOG(fatal) << "volNeuland not found";
    }

    neuland_geo_par_->SetNeulandGeoNode(geoNodeNeuland);
    neuland_geo_par_->setChanged();
}

void R3BNeuland::Register()
{
    LOG(info) << "R3BNeuland initialization ...";

    neuland_points_.init();

    WriteParameterFile();
    ResetValues();
}

ClassImp(R3BNeuland);
