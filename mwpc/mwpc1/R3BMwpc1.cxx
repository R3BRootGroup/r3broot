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

// ------------------------------------------------------------------------
// -----                        R3BMwpc1 source file                  -----
// -----                  Created 06/10/19 by JL Rodriguez            -----
// ------------------------------------------------------------------------

#include "FairGeoInterface.h"
#include "FairGeoLoader.h"
#include "FairVolume.h"
#include <FairRootManager.h>

#include "R3BMCStack.h"
#include "R3BMwpc1.h"
#include "R3BMwpcPoint.h"
#include "TClonesArray.h"
#include "TGeoManager.h"
#include "TObjArray.h"
#include "TParticle.h"
#include "TVirtualMC.h"
#include <stdlib.h>

R3BMwpc1::R3BMwpc1()
    : R3BMwpc1("")
{
}

R3BMwpc1::R3BMwpc1(const TString& geoFile, const TGeoTranslation& trans, const TGeoRotation& rot)
    : R3BMwpc1(geoFile, { trans, rot })
{
}

R3BMwpc1::R3BMwpc1(const TString& geoFile, const TGeoCombiTrans& combi)
    : R3BDetector("R3BMwpc1", kSOFMWPC1, geoFile, combi)
    , fSofMWPCCollection(new TClonesArray("R3BMwpcPoint"))
    , fPosIndex(0)
{
    ResetParameters();
}

R3BMwpc1::~R3BMwpc1()
{
    if (fSofMWPCCollection)
    {
        fSofMWPCCollection->Delete();
        delete fSofMWPCCollection;
    }
}

void R3BMwpc1::Initialize()
{
    FairDetector::Initialize();

    LOG(info) << "R3BMwpc1: initialisation";
    LOG(debug) << "R3BMwpc1: Sens. Vol. (McId) " << gMC->VolId("MWPC1");
}

void R3BMwpc1::SetSpecialPhysicsCuts() { LOG(info) << "R3BMwpc1: Adding customized Physics cut ... "; }

// -----   Public method ProcessHits  --------------------------------------
Bool_t R3BMwpc1::ProcessHits(FairVolume* vol)
{
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

    // Set additional parameters at exit of active volume. Create R3BMwpcPoint.
    if (gMC->IsTrackExiting() || gMC->IsTrackStop() || gMC->IsTrackDisappeared())
    {
        fTrackID = gMC->GetStack()->GetCurrentTrackNumber();
        fVolumeID = vol->getMCid();
        fDetCopyID = 1; // vol->getCopyNo();
        gMC->TrackPosition(fPosOut);
        gMC->TrackMomentum(fMomOut);
        if (fELoss == 0.)
            return kFALSE;

        if (gMC->IsTrackExiting())
        {
            const Double_t* oldpos;
            const Double_t* olddirection;
            Double_t newpos[3];
            Double_t newdirection[3];
            Double_t safety;

            gGeoManager->FindNode(fPosOut.X(), fPosOut.Y(), fPosOut.Z());
            oldpos = gGeoManager->GetCurrentPoint();
            olddirection = gGeoManager->GetCurrentDirection();

            for (Int_t i = 0; i < 3; i++)
            {
                newdirection[i] = -1 * olddirection[i];
            }

            gGeoManager->SetCurrentDirection(newdirection);
            //   TGeoNode *bla = gGeoManager->FindNextBoundary(2);
            safety = gGeoManager->GetSafeDistance();

            gGeoManager->SetCurrentDirection(-newdirection[0], -newdirection[1], -newdirection[2]);

            for (Int_t i = 0; i < 3; i++)
            {
                newpos[i] = oldpos[i] - (3 * safety * olddirection[i]);
            }

            fPosOut.SetX(newpos[0]);
            fPosOut.SetY(newpos[1]);
            fPosOut.SetZ(newpos[2]);
        }

        AddPoint(fTrackID,
                 fVolumeID,
                 fDetCopyID,
                 TVector3(fPosIn.X(), fPosIn.Y(), fPosIn.Z()),
                 TVector3(fPosOut.X(), fPosOut.Y(), fPosOut.Z()),
                 TVector3(fMomIn.Px(), fMomIn.Py(), fMomIn.Pz()),
                 TVector3(fMomOut.Px(), fMomOut.Py(), fMomOut.Pz()),
                 fTime,
                 fLength,
                 fELoss);

        // Increment number of TraPoints for this track
        R3BStack* stack = dynamic_cast<R3BStack*>(gMC->GetStack());
        stack->AddPoint(kSOFMWPC1);

        ResetParameters();
    }

    return kTRUE;
}

// -----   Public method EndOfEvent   -----------------------------------------
void R3BMwpc1::EndOfEvent()
{
    if (fVerboseLevel)
        Print();
    fSofMWPCCollection->Clear();

    ResetParameters();
}
// ----------------------------------------------------------------------------

// -----   Public method Register   -------------------------------------------
void R3BMwpc1::Register() { FairRootManager::Instance()->Register("Mwpc1Point", GetName(), fSofMWPCCollection, kTRUE); }
// ----------------------------------------------------------------------------

// -----   Public method GetCollection   --------------------------------------
TClonesArray* R3BMwpc1::GetCollection(Int_t iColl) const
{
    if (iColl == 0)
        return fSofMWPCCollection;
    else
        return NULL;
}
// ----------------------------------------------------------------------------

// -----   Public method Print   ----------------------------------------------
void R3BMwpc1::Print(Option_t* option) const
{
    Int_t nHits = fSofMWPCCollection->GetEntriesFast();
    LOG(info) << "R3BMwpc1: " << nHits << " points registered in this event";
}
// ----------------------------------------------------------------------------

// -----   Public method Reset   ----------------------------------------------
void R3BMwpc1::Reset()
{
    fSofMWPCCollection->Clear();
    ResetParameters();
}
// ----------------------------------------------------------------------------

// -----   Public method CopyClones   -----------------------------------------
void R3BMwpc1::CopyClones(TClonesArray* cl1, TClonesArray* cl2, Int_t offset)
{
    Int_t nEntries = cl1->GetEntriesFast();
    LOG(info) << "R3BMwpc1: " << nEntries << " entries to add";
    TClonesArray& clref = *cl2;
    R3BMwpcPoint* oldpoint = NULL;
    for (Int_t i = 0; i < nEntries; i++)
    {
        oldpoint = dynamic_cast<R3BMwpcPoint*>(cl1->At(i));
        Int_t index = oldpoint->GetTrackID() + offset;
        oldpoint->SetTrackID(index);
        new (clref[fPosIndex]) R3BMwpcPoint(*oldpoint);
        fPosIndex++;
    }
    LOG(info) << "R3BMwpc1: " << cl2->GetEntriesFast() << " merged entries";
}

// -----   Private method AddPoint   --------------------------------------------
R3BMwpcPoint* R3BMwpc1::AddPoint(Int_t trackID,
                                 Int_t detID,
                                 Int_t detCopyID,
                                 TVector3 posIn,
                                 TVector3 posOut,
                                 TVector3 momIn,
                                 TVector3 momOut,
                                 Double_t time,
                                 Double_t length,
                                 Double_t eLoss)
{
    TClonesArray& clref = *fSofMWPCCollection;
    Int_t size = clref.GetEntriesFast();
    if (fVerboseLevel > 1)
        LOG(info) << "R3BMwpc1: Adding Point at (" << posIn.X() << ", " << posIn.Y() << ", " << posIn.Z()
                  << ") cm,  detector " << detID << ", track " << trackID << ", energy loss " << eLoss * 1e06 << " keV";
    return new (clref[size]) R3BMwpcPoint(trackID, detID, detCopyID, posIn, posOut, momIn, momOut, time, length, eLoss);
}

Bool_t R3BMwpc1::CheckIfSensitive(std::string name)
{
    if (TString(name).Contains("MWPC1"))
    {
        LOG(debug) << "Found MWPC1 geometry from ROOT file: " << name;
        return kTRUE;
    }
    return kFALSE;
}

ClassImp(R3BMwpc1);
