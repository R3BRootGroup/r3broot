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

#include <FairLogger.h>
#include <FairRootManager.h>
#include <FairVolume.h>

#include <R3BCalifa.h>
#include <R3BCalifaGeometry.h>
#include <R3BCalifaPoint.h>
#include <R3BLogger.h>
#include <R3BMCStack.h>

#include <TClonesArray.h>
#include <TGeoManager.h>
#include <TGeoNode.h>
#include <TParticle.h>
#include <TVirtualMC.h>

#include <iostream>
#include <stdlib.h>

R3BCalifa::R3BCalifa()
    : R3BCalifa("")
{
}

R3BCalifa::R3BCalifa(const std::string& geoFile, const TGeoTranslation& trans, const TGeoRotation& rot)
    : R3BCalifa(geoFile, { trans, rot })
{
}

R3BCalifa::R3BCalifa(const std::string& geoFile, const TGeoCombiTrans& combi)
    : R3BDetector("R3BCalifa", kCALIFA, geoFile, combi)
{
    ResetParameters();
    fCalifaCollection = new TClonesArray("R3BCalifaPoint");
}

R3BCalifa::~R3BCalifa()
{
    if (fCalifaCollection)
    {
        fCalifaCollection->Delete();
        delete fCalifaCollection;
    }
}

void R3BCalifa::Initialize()
{
    FairDetector::Initialize();

    R3BLOG(info, " ");

    TGeoVolume* vol = static_cast<TGeoVolume*>(gGeoManager->GetVolume("CalifaWorld"));
    vol->SetVisibility(kFALSE);

    if (!R3BCalifaGeometry::Instance()->Init(fGeometryVersion))
    {
        R3BLOG(error, "Califa geometry not found");
    }
    return;
}

Bool_t R3BCalifa::ProcessHits(FairVolume* vol)
{
    int crystalId = R3BCalifaGeometry::Instance()->GetCrystalId(TVirtualMC::GetMC()->CurrentVolPath());

    if (TVirtualMC::GetMC()->IsTrackEntering())
    {
        fELoss = 0.;
        fTime = TVirtualMC::GetMC()->TrackTime() * 1.0e09;
        fLengthzero = TVirtualMC::GetMC()->TrackLength();
        TVirtualMC::GetMC()->TrackPosition(fPosIn);
        TVirtualMC::GetMC()->TrackMomentum(fMomIn);
    }

    // Sum energy loss for all steps in the active volume
    fELoss += TVirtualMC::GetMC()->Edep() * 1000.; // in MeV

    // Set additional parameters at exit of active volume. Create R3BCalifaPoint.
    if (TVirtualMC::GetMC()->IsTrackExiting() || TVirtualMC::GetMC()->IsTrackStop() ||
        TVirtualMC::GetMC()->IsTrackDisappeared())
    {
        if (fELoss == 0.)
        {
            return kFALSE;
        }

        fTrackID = TVirtualMC::GetMC()->GetStack()->GetCurrentTrackNumber();
        fVolumeID = vol->getMCid();
        fLength = TVirtualMC::GetMC()->TrackLength();
        fTrackPID = TVirtualMC::GetMC()->TrackPid();

        AddPoint(fTrackID,
                 fVolumeID,
                 fTrackPID,
                 crystalId,
                 TVector3(fPosIn.X(), fPosIn.Y(), fPosIn.Z()),
                 TVector3(fMomIn.Px(), fMomIn.Py(), fMomIn.Pz()),
                 fTime,
                 (fLength - fLengthzero),
                 fELoss);

        // Increment number of CalifaPoints for this track
        auto* stack = dynamic_cast<R3BStack*>(TVirtualMC::GetMC()->GetStack());
        stack->AddPoint(kCALIFA);
        ResetParameters();
    }
    return kTRUE;
}

void R3BCalifa::EndOfEvent()
{
    if (fVerboseLevel > 1)
    {
        Print();
    }
    fCalifaCollection->Clear();
    ResetParameters();
}

void R3BCalifa::Register()
{
    FairRootManager::Instance()->Register("CalifaPoint", GetName(), fCalifaCollection, kTRUE);
}

TClonesArray* R3BCalifa::GetCollection(Int_t iColl) const
{
    if (iColl == 0)
    {
        return fCalifaCollection;
    }
    else
    {
        return nullptr;
    }
}

void R3BCalifa::Print(Option_t* /*opt*/) const
{
    Int_t nPoints = fCalifaCollection->GetEntriesFast();
    LOG(info) << "R3BCalifa: " << nPoints << " points registered in this event";
}

void R3BCalifa::Reset()
{
    fCalifaCollection->Clear();
    ResetParameters();
}

void R3BCalifa::CopyClones(TClonesArray* cl1, TClonesArray* cl2, Int_t offset)
{
    Int_t nEntries = cl1->GetEntriesFast();
    LOG(info) << "R3BCalifa: " << nEntries << " entries to add";
    TClonesArray& clref = *cl2;
    R3BCalifaPoint* oldpoint = NULL;
    for (Int_t i = 0; i < nEntries; i++)
    {
        oldpoint = dynamic_cast<R3BCalifaPoint*>(cl1->At(i));
        Int_t index = oldpoint->GetTrackID() + offset;
        oldpoint->SetTrackID(index);
        new (clref[fPosIndex]) R3BCalifaPoint(*oldpoint);
        fPosIndex++;
    }
    LOG(info) << "R3BCalifa: " << cl2->GetEntriesFast() << " merged entries";
}

R3BCalifaPoint* R3BCalifa::AddPoint(int trackID,
                                    int detID,
                                    int trackPID,
                                    int cryID,
                                    TVector3 posIn,
                                    TVector3 momIn,
                                    double time,
                                    double length,
                                    double eLoss)
{
    TClonesArray& clref = *fCalifaCollection;
    Int_t size = clref.GetEntriesFast();
    if (fVerboseLevel > 1)
    {
        LOG(info) << "R3BCalifa: Adding Point at (" << posIn.X() << ", " << posIn.Y() << ", " << posIn.Z()
                  << ") cm,  crystal " << cryID << ", track " << trackPID << ", energy loss " << eLoss << " MeV";
    }
    return new (clref[size]) R3BCalifaPoint(trackID, detID, trackPID, cryID, posIn, momIn, time, length, eLoss);
}

void R3BCalifa::SelectGeometryVersion(Int_t version) { fGeometryVersion = version; }

Bool_t R3BCalifa::CheckIfSensitive(std::string name)
{
    if (TString(name).Contains("Crystal_"))
    {
        return kTRUE;
    }
    return kFALSE;
}

ClassImp(R3BCalifa);
