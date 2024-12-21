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

#include "R3BGfi.h"
#include "FairGeoInterface.h"
#include "FairGeoLoader.h"
#include "FairGeoNode.h"
#include "FairGeoRootBuilder.h"
#include "FairRootManager.h"
#include "FairRun.h"
#include "FairRuntimeDb.h"
#include "FairVolume.h"
#include "R3BGeoGfi.h"
#include "R3BGeoGfiPar.h"
#include "R3BGfiPoint.h"
#include "R3BMCStack.h"
#include "TClonesArray.h"
#include "TGeoMCGeometry.h"
#include "TGeoManager.h"
#include "TObjArray.h"
#include "TParticle.h"
#include "TVirtualMC.h"
#include <stdlib.h>

R3BGfi::R3BGfi(const TString& geoFile)
    : R3BDetector("R3BGfi", kGFI)
    , fPos1(-73.274339, 0.069976, 513.649524)
    , fPos2(-147.135037, 0.069976, 729.680342)
    , fRot1(new TGeoRotation())
    , fRot2(new TGeoRotation())
{
    ResetParameters();
    SetGeometryFileName(geoFile);
    fRot1->RotateY(-16.7);
    fRot2->RotateY(-16.7);
    fGfiCollection = new TClonesArray("R3BGfiPoint");
    fPosIndex = 0;
    kGeoSaved = kFALSE;
    flGeoPar = new TList();
    flGeoPar->SetName(GetName());
    fVerboseLevel = 1;
}

R3BGfi::R3BGfi(const TString& geoFile,
               const Double_t x1,
               const Double_t y1,
               const Double_t z1,
               const Double_t rot_y1,
               const Double_t x2,
               const Double_t y2,
               const Double_t z2,
               const Double_t rot_y2)
    : R3BDetector("R3BGfi", kGFI)
    , fPos1(x1, y1, z1)
    , fPos2(x2, y2, z2)
    , fRot1(new TGeoRotation())
    , fRot2(new TGeoRotation())
{
    ResetParameters();
    SetGeometryFileName(geoFile);
    fRot1->RotateY(rot_y1);
    fRot2->RotateY(rot_y2);
    fGfiCollection = new TClonesArray("R3BGfiPoint");
    fPosIndex = 0;
    kGeoSaved = kFALSE;
    flGeoPar = new TList();
    flGeoPar->SetName(GetName());
    fVerboseLevel = 1;
}

R3BGfi::~R3BGfi()
{
    if (flGeoPar)
    {
        delete flGeoPar;
    }
    if (fGfiCollection)
    {
        fGfiCollection->Delete();
        delete fGfiCollection;
    }
}

void R3BGfi::Initialize()
{
    FairDetector::Initialize();

    LOG(info) << "R3BGfi: initialisation";
    LOG(debug) << "R3BGfi: Sci. Vol. (McId) " << gMC->VolId("GFILog");
}

void R3BGfi::SetSpecialPhysicsCuts()
{
    LOG(info) << "-I- R3BGfi: Adding customized Physics cut ... ";

    if (gGeoManager)
    {
        TGeoMedium* pSi = gGeoManager->GetMedium("plasticForGFI");
        if (pSi)
        {
            // Setting processes for Si only
            gMC->Gstpar(pSi->GetId(), "LOSS", 3);
            gMC->Gstpar(pSi->GetId(), "STRA", 1.0);
            gMC->Gstpar(pSi->GetId(), "PAIR", 1.0);
            gMC->Gstpar(pSi->GetId(), "COMP", 1.0);
            gMC->Gstpar(pSi->GetId(), "PHOT", 1.0);
            gMC->Gstpar(pSi->GetId(), "ANNI", 1.0);
            gMC->Gstpar(pSi->GetId(), "BREM", 1.0);
            gMC->Gstpar(pSi->GetId(), "HADR", 1.0);
            gMC->Gstpar(pSi->GetId(), "DRAY", 1.0);
            gMC->Gstpar(pSi->GetId(), "DCAY", 1.0);
            gMC->Gstpar(pSi->GetId(), "MULS", 1.0);
            gMC->Gstpar(pSi->GetId(), "RAYL", 1.0);

            // Setting Energy-CutOff for Si Only
            Double_t cutE = fCutE; // GeV-> 1 keV

            LOG(info) << "-I- R3BGfi: plasticForGFI Medium Id " << pSi->GetId() << " Energy Cut-Off : " << cutE
                      << " GeV";
            // Si
            gMC->Gstpar(pSi->GetId(), "CUTGAM", cutE); /** gammas (GeV)*/
            gMC->Gstpar(pSi->GetId(), "CUTELE", cutE); /** electrons (GeV)*/
            gMC->Gstpar(pSi->GetId(), "CUTNEU", cutE); /** neutral hadrons (GeV)*/
            gMC->Gstpar(pSi->GetId(), "CUTHAD", cutE); /** charged hadrons (GeV)*/
            gMC->Gstpar(pSi->GetId(), "CUTMUO", cutE); /** muons (GeV)*/
            gMC->Gstpar(pSi->GetId(), "BCUTE", cutE);  /** electron bremsstrahlung (GeV)*/
            gMC->Gstpar(pSi->GetId(), "BCUTM", cutE);  /** muon and hadron bremsstrahlung(GeV)*/
            gMC->Gstpar(pSi->GetId(), "DCUTE", cutE);  /** delta-rays by electrons (GeV)*/
            gMC->Gstpar(pSi->GetId(), "DCUTM", cutE);  /** delta-rays by muons (GeV)*/
            gMC->Gstpar(pSi->GetId(), "PPCUTM", -1.);  /** direct pair production by muons (GeV)*/
        }

    } //! gGeoManager
}

// -----   Public method ProcessHits  --------------------------------------
Bool_t R3BGfi::ProcessHits(FairVolume* vol)
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

    // Set additional parameters at exit of active volume. Create R3BGfiPoint.
    if (gMC->IsTrackExiting() || gMC->IsTrackStop() || gMC->IsTrackDisappeared())
    {
        fTrackID = gMC->GetStack()->GetCurrentTrackNumber();
        fVolumeID = vol->getMCid();
        gMC->TrackPosition(fPosOut);
        gMC->TrackMomentum(fMomOut);
        //        if (fELoss == 0.)
        //            return kFALSE;

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

            //       cout << "1st direction: " << olddirection[0] << "," << olddirection[1] << "," << olddirection[2] <<
            //       endl;

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

            if (fPosIn.Z() < 30. && newpos[2] > 30.02)
            {
                LOG(error) << "2nd direction: " << olddirection[0] << "," << olddirection[1] << "," << olddirection[2]
                           << " with safety = " << safety;
                LOG(error) << "oldpos = " << oldpos[0] << "," << oldpos[1] << "," << oldpos[2];
                LOG(error) << "newpos = " << newpos[0] << "," << newpos[1] << "," << newpos[2];
            }

            fPosOut.SetX(newpos[0]);
            fPosOut.SetY(newpos[1]);
            fPosOut.SetZ(newpos[2]);
        }

        AddHit(fTrackID,
               fVolumeID,
               planeNr,
               TVector3(fPosIn.X(), fPosIn.Y(), fPosIn.Z()),
               TVector3(fPosOut.X(), fPosOut.Y(), fPosOut.Z()),
               TVector3(fMomIn.Px(), fMomIn.Py(), fMomIn.Pz()),
               TVector3(fMomOut.Px(), fMomOut.Py(), fMomOut.Pz()),
               fTime,
               fLength,
               fELoss);

        // Increment number of GfiPoints for this track
        R3BStack* stack = dynamic_cast<R3BStack*>(gMC->GetStack());
        stack->AddPoint(kGFI);

        ResetParameters();
    }

    return kTRUE;
}

// -----   Public method EndOfEvent   -----------------------------------------
void R3BGfi::BeginEvent()
{

    //  if (! kGeoSaved ) {
    //      SaveGeoParams();
    //  cout << "-I STS geometry parameters saved " << endl;
    //  kGeoSaved = kTRUE;
    //  }
}
// -----   Public method EndOfEvent   -----------------------------------------
void R3BGfi::EndOfEvent()
{

    if (fVerboseLevel)
        Print();
    fGfiCollection->Clear();

    ResetParameters();
}
// ----------------------------------------------------------------------------

// -----   Public method Register   -------------------------------------------
void R3BGfi::Register() { FairRootManager::Instance()->Register("GFIPoint", GetName(), fGfiCollection, kTRUE); }
// ----------------------------------------------------------------------------

// -----   Public method GetCollection   --------------------------------------
TClonesArray* R3BGfi::GetCollection(Int_t iColl) const
{
    if (iColl == 0)
        return fGfiCollection;
    else
        return NULL;
}
// ----------------------------------------------------------------------------

// -----   Public method Print   ----------------------------------------------
void R3BGfi::Print(Option_t* option) const
{
    Int_t nHits = fGfiCollection->GetEntriesFast();
    LOG(info) << "R3BGfi: " << nHits << " points registered in this event";
}
// ----------------------------------------------------------------------------

// -----   Public method Reset   ----------------------------------------------
void R3BGfi::Reset()
{
    fGfiCollection->Clear();
    ResetParameters();
}
// ----------------------------------------------------------------------------

// -----   Public method CopyClones   -----------------------------------------
void R3BGfi::CopyClones(TClonesArray* cl1, TClonesArray* cl2, Int_t offset)
{
    Int_t nEntries = cl1->GetEntriesFast();
    LOG(info) << "R3BGfi: " << nEntries << " entries to add";
    TClonesArray& clref = *cl2;
    R3BGfiPoint* oldpoint = NULL;
    for (Int_t i = 0; i < nEntries; i++)
    {
        oldpoint = dynamic_cast<R3BGfiPoint*>(cl1->At(i));
        Int_t index = oldpoint->GetTrackID() + offset;
        oldpoint->SetTrackID(index);
        new (clref[fPosIndex]) R3BGfiPoint(*oldpoint);
        fPosIndex++;
    }
    LOG(info) << "R3BGfi: " << cl2->GetEntriesFast() << " merged entries";
}

// -----   Private method AddHit   --------------------------------------------
R3BGfiPoint* R3BGfi::AddHit(Int_t trackID,
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
    TClonesArray& clref = *fGfiCollection;
    Int_t size = clref.GetEntriesFast();
    if (fVerboseLevel > 1)
        LOG(info) << "R3BGfi: Adding Point at (" << posIn.X() << ", " << posIn.Y() << ", " << posIn.Z()
                  << ") cm,  detector " << detID << ", track " << trackID << ", energy loss " << eLoss * 1e06 << " keV";
    return new (clref[size]) R3BGfiPoint(trackID, detID, plane, posIn, posOut, momIn, momOut, time, length, eLoss);
}

// -----   Public method ConstructGeometry   ----------------------------------
void R3BGfi::ConstructGeometry()
{
    TString fileName = GetGeometryFileName();
    if (fileName.EndsWith(".root"))
    {
        LOG(info) << "Constructing GFI geometry from ROOT file " << fileName.Data();
        ConstructRootGeometry();

        TGeoNode* gfi_node = gGeoManager->GetTopVolume()->GetNode("GFI_0");

        TGeoNode* node = gfi_node->GetVolume()->GetNode("GFILogWorld_0");
        TGeoCombiTrans* combtrans = dynamic_cast<TGeoCombiTrans*>((dynamic_cast<TGeoNodeMatrix*>(node))->GetMatrix());
        TGeoTranslation* tr1 = new TGeoTranslation(fPos1.X(), fPos1.Y(), fPos1.Z());
        *combtrans = { *tr1, *fRot1 };

        node = gfi_node->GetVolume()->GetNode("GFILogWorld_1");
        combtrans = dynamic_cast<TGeoCombiTrans*>((dynamic_cast<TGeoNodeMatrix*>(node))->GetMatrix());
        TGeoTranslation* tr2 = new TGeoTranslation(fPos2.X(), fPos2.Y(), fPos2.Z());
        *combtrans = { *tr2, *fRot2 };
    }
    else
    {
        LOG(fatal) << "GFI geometry file is not specified";
        exit(1);
    }
}

Bool_t R3BGfi::CheckIfSensitive(std::string name)
{
    if (TString(name).Contains("GFILog"))
    {
        return kTRUE;
    }
    return kFALSE;
}

ClassImp(R3BGfi)
