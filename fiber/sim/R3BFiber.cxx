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

// ---------------------------------------------------------------
// -----                       R3BFiber                      -----
// -----          Created 10/05/2022 by J.L. Rodriguez       -----
// ---------------------------------------------------------------

#include "R3BFiber.h"
#include "R3BFibPoint.h"
#include "R3BLogger.h"
#include "R3BMCStack.h"

#include "FairRootManager.h"
#include "FairRuntimeDb.h"
#include "FairVolume.h"

#include "TClonesArray.h"
#include "TGeoManager.h"
#include "TParticle.h"
#include "TVirtualMC.h"

R3BFiber::R3BFiber()
    : R3BFiber("", "", kREF)
{
}

R3BFiber::R3BFiber(const TString& name,
                   const TString& geoFile,
                   DetectorId detid,
                   const TGeoTranslation& trans,
                   const TGeoRotation& rot)
    : R3BFiber(name, geoFile, detid, { trans, rot })
{
}

R3BFiber::R3BFiber(const TString& name, const TString& geoFile, DetectorId detid, const TGeoCombiTrans& combi)
    : R3BDetector("R3B" + name, detid, geoFile, combi)
    , fName(name)
    , fDetId(detid)
    , fFiCollection(new TClonesArray("R3BFibPoint"))
    , fPosIndex(0)
    , kGeoSaved(kFALSE)
    , flGeoPar(new TList())
{
    flGeoPar->SetName(GetName());
    ResetParameters();
}

R3BFiber::~R3BFiber()
{
    if (flGeoPar)
    {
        delete flGeoPar;
    }
    if (fFiCollection)
    {
        fFiCollection->Delete();
        delete fFiCollection;
    }
}

void R3BFiber::Initialize()
{
    FairDetector::Initialize();

    R3BLOG(info, "for fiber " << fName);
    R3BLOG(debug, "Vol (McId) def " << gMC->VolId("Fiber"));
}

void R3BFiber::SetSpecialPhysicsCuts() { R3BLOG(info, ""); }

// -----   Public method ProcessHits  --------------------------------------
Bool_t R3BFiber::ProcessHits(FairVolume* vol)
{
    // 2 Simple Det PLane
    // get Info from DCH planes
    Int_t copyNo = -1;
    Int_t planeNr = -1;
    // Get the Geo info from MC Point
    gMC->CurrentVolID(copyNo);
    gMC->CurrentVolOffID(1, planeNr);

    if (gMC->IsTrackEntering())
    {
        fELoss = 0.;
        // fTime   = gMC->TrackTime() * 1.0e09;
        // fLength = gMC->TrackLength();
        fTime_in = gMC->TrackTime() * 1.0e09;
        fLength_in = gMC->TrackLength();
        gMC->TrackPosition(fPosIn);
        gMC->TrackMomentum(fMomIn);
    }

    // Sum energy loss for all steps in the active volume
    fELoss += gMC->Edep();

    // Set additional parameters at exit of active volume. Create R3BFiberPoint.
    if (gMC->IsTrackExiting() || gMC->IsTrackStop() || gMC->IsTrackDisappeared())
    {
        fTrackID = gMC->GetStack()->GetCurrentTrackNumber();
        fVolumeID = vol->getMotherCopyNo();
        gMC->TrackPosition(fPosOut);
        gMC->TrackMomentum(fMomOut);
        if (fELoss == 0.)
            return kFALSE;

        fTime_out = gMC->TrackTime() * 1.0e09; // also in case particle is stopped in detector, or decays...
        fLength_out = gMC->TrackLength();
        fTime = (fTime_out + fTime_in) / 2.;
        fLength = (fLength_out + fLength_in) / 2.;

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
            // TGeoNode *bla = gGeoManager->FindNextBoundary(2);
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

        AddHit(fTrackID,
               /*fVolumeID*/ /*copyNo*/ planeNr,
               planeNr,
               TVector3(fPosIn.X(), fPosIn.Y(), fPosIn.Z()),
               TVector3(fPosOut.X(), fPosOut.Y(), fPosOut.Z()),
               TVector3(fMomIn.Px(), fMomIn.Py(), fMomIn.Pz()),
               TVector3(fMomOut.Px(), fMomOut.Py(), fMomOut.Pz()),
               fTime,
               fLength,
               fELoss);

        // Increment number of Fi30Points for this track
        R3BStack* stack = dynamic_cast<R3BStack*>(gMC->GetStack());
        stack->AddPoint(fDetId);

        ResetParameters();
    }

    return kTRUE;
}

// -----   Public method EndOfEvent   -----------------------------------------
void R3BFiber::EndOfEvent()
{
    if (fVerboseLevel)
        Print();
    fFiCollection->Clear();

    ResetParameters();
}
// ----------------------------------------------------------------------------

// -----   Public method Register   -------------------------------------------
void R3BFiber::Register()
{
    R3BLOG(debug, "");
    FairRootManager::Instance()->Register(fName + "Point", GetName(), fFiCollection, kTRUE);
}
// ----------------------------------------------------------------------------

// -----   Public method GetCollection   --------------------------------------
TClonesArray* R3BFiber::GetCollection(Int_t iColl) const
{
    if (iColl == 0)
        return fFiCollection;
    else
        return NULL;
}
// ----------------------------------------------------------------------------

// -----   Public method Print   ----------------------------------------------
void R3BFiber::Print(Option_t* option) const
{
    Int_t nHits = fFiCollection->GetEntriesFast();
    LOG(info) << "R3B" << fName << ": " << nHits << " points registered in this event";
}
// ----------------------------------------------------------------------------

// -----   Public method Reset   ----------------------------------------------
void R3BFiber::Reset()
{
    fFiCollection->Clear();
    ResetParameters();
}
// ----------------------------------------------------------------------------

// -----   Public method CopyClones   -----------------------------------------
void R3BFiber::CopyClones(TClonesArray* cl1, TClonesArray* cl2, Int_t offset)
{
    Int_t nEntries = cl1->GetEntriesFast();
    R3BLOG(info, nEntries << " entries to add");
    TClonesArray& clref = *cl2;
    R3BFibPoint* oldpoint = NULL;
    for (Int_t i = 0; i < nEntries; i++)
    {
        oldpoint = dynamic_cast<R3BFibPoint*>(cl1->At(i));
        Int_t index = oldpoint->GetTrackID() + offset;
        oldpoint->SetTrackID(index);
        new (clref[fPosIndex]) R3BFibPoint(*oldpoint);
        fPosIndex++;
    }
    R3BLOG(info, cl2->GetEntriesFast() << " merged entries");
}

// -----   Private method AddHit   --------------------------------------------
R3BFibPoint* R3BFiber::AddHit(Int_t trackID,
                              Int_t detID,
                              Int_t plane,
                              TVector3 posIn,
                              TVector3 posOut,
                              TVector3 momIn,
                              TVector3 momOut,
                              Double_t time,
                              Double_t length,
                              Double_t eLoss)
{
    TClonesArray& clref = *fFiCollection;
    Int_t size = clref.GetEntriesFast();
    if (fVerboseLevel > 1)
    {
        R3BLOG(info,
               "Adding Point at (" << posIn.X() << ", " << posIn.Y() << ", " << posIn.Z() << ") cm,  detector " << detID
                                   << ", track " << trackID << ", energy loss " << eLoss * 1e06 << " keV");
    }
    return new (clref[size]) R3BFibPoint(trackID, detID, plane, posIn, posOut, momIn, momOut, time, length, eLoss);
}

Bool_t R3BFiber::CheckIfSensitive(std::string name)
{
    if (TString(name).Contains(fName))
    {
        LOG(debug) << "Found" << fName << " geometry from ROOT file: " << name;
        return kTRUE;
    }
    return kFALSE;
}

ClassImp(R3BFiber);
