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

// Initialize variables from Birk' s Law
static constexpr Double_t BirkdP = 1.032;
static constexpr Double_t BirkC1 = 0.013 / BirkdP;
static constexpr Double_t BirkC2 = 9.6e-6 / (BirkdP * BirkdP);

inline Double_t GetLightYield(const Int_t charge, const Double_t length, const Double_t edep)
{
    // Apply Birk's law ( Adapted from G3BIRK/Geant3)
    if (charge != 0 && length > 0)
    {
        Double_t birkC1Mod = BirkC1;

        // Apply correction for higher charge states
        if (TMath::Abs(charge) >= 2)
        {
            birkC1Mod *= 7.2 / 12.6;
        }

        Double_t dedxcm = 1000. * edep / length;
        Double_t lightYield = edep / (1. + birkC1Mod * dedxcm + BirkC2 * dedxcm * dedxcm);
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
    , fNeulandPoints(new TClonesArray("R3BNeulandPoint"))
{
}

R3BNeuland::R3BNeuland(Int_t nDP, const TGeoTranslation& trans, const TGeoRotation& rot)
    : R3BNeuland(nDP, { trans, rot })
{
}

R3BNeuland::R3BNeuland(const Int_t nDP, const TGeoCombiTrans& combi)
    : R3BNeuland(TString::Format("neuland_v3_%ddp.geo.root", nDP), combi)
{
}

R3BNeuland::~R3BNeuland()
{
    if (fNeulandPoints)
    {
        fNeulandPoints->Delete();
        delete fNeulandPoints;
    }
}

void R3BNeuland::Initialize()
{
    LOG(info) << "R3BNeuland initialization ...";

    FairDetector::Initialize();

    WriteParameterFile();
    ResetValues();
}

Bool_t R3BNeuland::ProcessHits(FairVolume*)
{
    // New hit in detector
    if (gMC->IsTrackEntering())
    {
        if (!fLastHitDone)
        {
            LOG(warn) << "R3BNeuland: Incomplete hit discarded";
            ResetValues();
        }

        fLastHitDone = kFALSE;
        fELoss = 0.;
        fLightYield = 0.;
        fTime = gMC->TrackTime() * 1.0e09;
        fLength = gMC->TrackLength();
        gMC->TrackPosition(fPosIn);
        gMC->TrackMomentum(fMomIn);
        gMC->CurrentVolOffID(1, fPaddleID);
    }

    // Sum energy loss for all steps in the active volume
    fELoss += gMC->Edep();
    fLightYield += GetLightYield(gMC->TrackCharge(), gMC->TrackStep(), gMC->Edep());

    // Set additional parameters at exit of active volume. Create R3BNeulandPoint.
    if (gMC->IsTrackExiting() || gMC->IsTrackStop() || gMC->IsTrackDisappeared())
    {
        // Do not save a hit if no energy deposited
        if (fELoss < 1e-20 || fLightYield < 1e-20)
        {
            ResetValues();
            return kTRUE;
        }

        fTrackID = gMC->GetStack()->GetCurrentTrackNumber();
        gMC->TrackPosition(fPosOut);
        gMC->TrackMomentum(fMomOut);

        // Add Point
        LOG(debug) << "R3BNeuland: Adding Point at (" << fPosIn.X() << ", " << fPosIn.Y() << ", " << fPosIn.Z()
                   << ") cm,  paddle " << fPaddleID << ", track " << fTrackID << ", energy loss " << fELoss << " GeV "
                   << gMC->GetStack()->GetCurrentParentTrackNumber();

        Int_t size = fNeulandPoints->GetEntriesFast();
        new ((*fNeulandPoints)[size]) R3BNeulandPoint(fTrackID,
                                                      fPaddleID,
                                                      fPosIn.Vect(),
                                                      fMomIn.Vect(),
                                                      fTime,
                                                      fLength,
                                                      fELoss,
                                                      gMC->CurrentEvent(),
                                                      fLightYield);

        // Increment number of LandPoints for this track
        auto stack = dynamic_cast<R3BStack*>(gMC->GetStack());
        stack->AddPoint(kNEULAND);
        ResetValues();
    }

    return kTRUE;
}

Bool_t R3BNeuland::CheckIfSensitive(std::string name) { return name == "volBC408"; }

void R3BNeuland::EndOfEvent()
{
    if (fVerboseLevel)
    {
        Print();
    }
    Reset();
}

TClonesArray* R3BNeuland::GetCollection(Int_t iColl) const
{
    if (iColl == 0)
    {
        return fNeulandPoints;
    }
    return nullptr;
}

void R3BNeuland::Register()
{
    FairRootManager::Instance()->Register("NeulandPoints", GetName(), fNeulandPoints, kTRUE);
}

void R3BNeuland::Print(Option_t*) const
{
    LOG(info) << "R3BNeuland: " << fNeulandPoints->GetEntries() << " Neuland Points registered in this event";
}

void R3BNeuland::Reset()
{
    fNeulandPoints->Clear();
    ResetValues();
}

void R3BNeuland::ResetValues()
{
    fLastHitDone = kTRUE;
    fTrackID = 0;
    fPaddleID = -1;
    fPosIn.Clear();
    fPosOut.Clear();
    fMomIn.Clear();
    fMomOut.Clear();
    fTime = fLength = fELoss = fLightYield = 0;
}

void R3BNeuland::WriteParameterFile()
{
    FairRuntimeDb* rtdb = FairRun::Instance()->GetRuntimeDb();
    fNeulandGeoPar = dynamic_cast<R3BNeulandGeoPar*>(rtdb->getContainer("R3BNeulandGeoPar"));

    // Really bad way to find the Neuland *node* (not the volume!)
    TGeoNode* geoNodeNeuland = nullptr;
    for (Int_t i = 0; i < gGeoManager->GetTopNode()->GetNdaughters(); i++)
    {
        if (TString(gGeoManager->GetTopNode()->GetDaughter(i)->GetVolume()->GetName()) == "volNeuland")
        {
            geoNodeNeuland = gGeoManager->GetTopNode()->GetDaughter(i);
            break;
        }
    }

    if (!geoNodeNeuland)
    {
        LOG(fatal) << "volNeuland not found";
    }

    fNeulandGeoPar->SetNeulandGeoNode(geoNodeNeuland);
    fNeulandGeoPar->setChanged();
}

ClassImp(R3BNeuland);
