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

#include "FairLogger.h"
#include "FairRootManager.h"
#include "FairVolume.h"

#include "R3BLogger.h"
#include "R3BMCStack.h"
#include "R3BRpc.h"
#include "R3BRpcPoint.h"

#include "TClonesArray.h"
#include "TGeoManager.h"
#include "TVirtualMC.h"

R3BRpc::R3BRpc()
    : R3BRpc("")
{
}

R3BRpc::R3BRpc(const TString& geoFile, const TGeoTranslation& trans, const TGeoRotation& rot)
    : R3BRpc(geoFile, { trans, rot })
{
}

R3BRpc::R3BRpc(const TString& geoFile, const TGeoCombiTrans& combi)
    : R3BDetector("R3BRpc", kRPC, geoFile, combi)
    , fRpcCollection(new TClonesArray("R3BRpcPoint"))
    , fPosIndex(0)
{
    ResetParameters();
}

R3BRpc::~R3BRpc()
{
    if (fRpcCollection)
    {
        fRpcCollection->Delete();
        delete fRpcCollection;
    }
}

void R3BRpc::Initialize()
{
    FairDetector::Initialize();

    R3BLOG(info, " ");
    R3BLOG(debug, "Vol. (McId) " << gMC->VolId("strip"));
}

Bool_t R3BRpc::ProcessHits(FairVolume* vol)
{
    // HERE YOU GET the CHANNEL FROM THE VOLUME PATH (BY NAME)
    // Int_t channelId = AFUNCTION OF gMC->CurrentVolPath();
    // OR ALTERNATIVELY FROM THE POSITION OF THE GEANT4 HIT

    if (gMC->IsTrackEntering())
    {
        fELoss = 0.;
        fTime = gMC->TrackTime() * 1.0e09;
        fLength = gMC->TrackLength();
        gMC->TrackPosition(fPosIn);
        gMC->TrackMomentum(fMomIn);
    }

    // Sum energy loss for all steps in the active volume
    fELoss += gMC->Edep();

    // Set additional parameters at exit of active volume. Create R3BRpcPoint.
    if (gMC->IsTrackExiting() || gMC->IsTrackStop() || gMC->IsTrackDisappeared())
    {
        fTrackID = gMC->GetStack()->GetCurrentTrackNumber();
        // fParentTrackID = gMC->GetStack()->GetCurrentParentTrackNumber();
        fVolumeID = vol->getMCid();
        fStripID = vol->getCopyNo();
        // fTrackPID = gMC->TrackPid();
        // fUniqueID = gMC->GetStack()->GetCurrentTrack()->GetUniqueID();
        // updating the value of the track length when exiting or stopping

        if (fELoss == 0.)
            return kFALSE;

        AddPoint(fTrackID,
                 fVolumeID,
                 fStripID,
                 TVector3(fPosIn.X(), fPosIn.Y(), fPosIn.Z()),
                 TVector3(fMomIn.Px(), fMomIn.Py(), fMomIn.Pz()),
                 fTime,
                 fLength,
                 fELoss,
                 gMC->CurrentEvent());

        // Increment number of RpcPoints for this track
        R3BStack* stack = dynamic_cast<R3BStack*>(gMC->GetStack());
        stack->AddPoint(kRPC);
        ResetParameters();
    }
    return kTRUE;
}

void R3BRpc::EndOfEvent()
{
    if (fVerboseLevel)
    {
        Print();
    }
    Reset();
}

void R3BRpc::Register() { FairRootManager::Instance()->Register("RpcPoint", GetName(), fRpcCollection, kTRUE); }

TClonesArray* R3BRpc::GetCollection(Int_t iColl) const
{
    if (iColl == 0)
    {
        return fRpcCollection;
    }
    else
    {
        return nullptr;
    }
}

void R3BRpc::Print(Option_t* option) const
{
    Int_t nPoints = fRpcCollection->GetEntriesFast();
    R3BLOG(info, nPoints << " points registered in this event");
}

void R3BRpc::Reset()
{
    fRpcCollection->Clear();
    ResetParameters();
}

void R3BRpc::CopyClones(TClonesArray* cl1, TClonesArray* cl2, Int_t offset)
{
    Int_t nEntries = cl1->GetEntriesFast();
    LOG(info) << "R3BRpc: " << nEntries << " entries to add";
    TClonesArray& clref = *cl2;
    R3BRpcPoint* oldpoint = NULL;
    for (Int_t i = 0; i < nEntries; i++)
    {
        oldpoint = dynamic_cast<R3BRpcPoint*>(cl1->At(i));
        Int_t index = oldpoint->GetTrackID() + offset;
        oldpoint->SetTrackID(index);
        new (clref[fPosIndex]) R3BRpcPoint(*oldpoint);
        fPosIndex++;
    }
    LOG(info) << "R3BRpc: " << cl2->GetEntriesFast() << " merged entries";
}

R3BRpcPoint* R3BRpc::AddPoint(Int_t trackID,
                              Int_t detID,
                              Int_t ident,
                              TVector3 posIn,
                              TVector3 momIn,
                              Double_t time,
                              Double_t length,
                              Double_t eLoss,
                              UInt_t EventId)
{
    TClonesArray& clref = *fRpcCollection;
    Int_t size = clref.GetEntriesFast();
    if (fVerboseLevel > 1)
    {
        LOG(info) << "R3BRpc: Adding Point at (" << posIn.X() << ", " << posIn.Y() << ", " << posIn.Z()
                  << ") cm,  detector " << detID << ", track " << trackID << ", energy loss " << eLoss * 1e06 << " keV";
    }
    return new (clref[size]) R3BRpcPoint(trackID, detID, ident, posIn, momIn, time, length, eLoss, EventId);
}

Bool_t R3BRpc::CheckIfSensitive(std::string name)
{
    if (TString(name).Contains("strip"))
    {
        // LOG(info) << "Found geometry from ROOT file: " << name;
        return kTRUE;
    }
    return kFALSE;
}

ClassImp(R3BRpc);
