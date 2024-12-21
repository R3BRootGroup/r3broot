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

#include "R3BPsp.h"
#include "R3BTGeoPar.h"

#include "FairGeoInterface.h"
#include "FairGeoLoader.h"
#include "FairGeoNode.h"
#include "FairGeoRootBuilder.h"
#include "FairRootManager.h"
#include "FairRun.h"
#include "FairRuntimeDb.h"
#include "FairVolume.h"
#include "R3BGeoPsp.h"
#include "R3BMCStack.h"
#include "R3BPspPoint.h"
#include "TClonesArray.h"
#include "TGeoArb8.h"
#include "TGeoBBox.h"
#include "TGeoBoolNode.h"
#include "TGeoCompositeShape.h"
#include "TGeoCone.h"
#include "TGeoMCGeometry.h"
#include "TGeoManager.h"
#include "TGeoMaterial.h"
#include "TGeoMatrix.h"
#include "TGeoMedium.h"
#include "TGeoPara.h"
#include "TGeoPgon.h"
#include "TGeoSphere.h"
#include "TObjArray.h"
#include "TParticle.h"
#include "TVirtualMC.h"
#include <stdlib.h>

R3BPsp::R3BPsp()
    : R3BPsp("")
{
}

R3BPsp::R3BPsp(const TString& geoFile,
               const TGeoTranslation& trans,
               const TGeoRotation& rot,
               const Float_t z1,
               const Float_t z2,
               const Float_t z3)
    : R3BPsp(geoFile, { trans, rot }, z1, z2, z3)
{
}

R3BPsp::R3BPsp(const TString& geoFile,
               const TGeoCombiTrans& combi,
               const Float_t z1,
               const Float_t z2,
               const Float_t z3)
    : R3BDetector("R3BPsp", kPSP, geoFile, combi)
    , fZ1(z1)
    , fZ2(z2)
    , fZ3(z3)
    , fPspCollection(new TClonesArray("R3BPspPoint"))
    , fPosIndex(0)
    , kGeoSaved(kFALSE)
    , flGeoPar(new TList())
{
    flGeoPar->SetName(GetName());
    ResetParameters();
}

R3BPsp::~R3BPsp()
{
    if (flGeoPar)
    {
        delete flGeoPar;
    }
    if (fPspCollection)
    {
        fPspCollection->Delete();
        delete fPspCollection;
    }
}

void R3BPsp::Initialize()
{
    FairDetector::Initialize();

    LOG(info) << "R3BPsp: initialisation";
    LOG(debug) << "R3BPsp: Vol. (McId) " << gMC->VolId("PSP1Log");
}

void R3BPsp::SetSpecialPhysicsCuts()
{
    LOG(info) << "-I- R3BPsp: Adding customized Physics cut ... ";

    if (gGeoManager)
    {
        TGeoMedium* pSi = gGeoManager->GetMedium("silicon");
        if (pSi && 0 == 1)
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

            LOG(info) << "-I- R3BPsp: silicon Medium Id " << pSi->GetId() << " Energy Cut-Off : " << cutE << " GeV";

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
Bool_t R3BPsp::ProcessHits(FairVolume* vol)
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

    // Set additional parameters at exit of active volume. Create R3BPspPoint.
    if (gMC->IsTrackExiting() || gMC->IsTrackStop() || gMC->IsTrackDisappeared())
    {
        fTrackID = gMC->GetStack()->GetCurrentTrackNumber();
        fVolumeID = vol->getMCid();
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
               fVolumeID,
               planeNr,
               TVector3(fPosIn.X(), fPosIn.Y(), fPosIn.Z()),
               TVector3(fPosOut.X(), fPosOut.Y(), fPosOut.Z()),
               TVector3(fMomIn.Px(), fMomIn.Py(), fMomIn.Pz()),
               TVector3(fMomOut.Px(), fMomOut.Py(), fMomOut.Pz()),
               fTime,
               fLength,
               fELoss);

        // Increment number of PspPoints for this track
        R3BStack* stack = dynamic_cast<R3BStack*>(gMC->GetStack());
        stack->AddPoint(kPSP);

        ResetParameters();
    }

    return kTRUE;
}

// -----   Public method EndOfEvent   -----------------------------------------
void R3BPsp::BeginEvent() {}

// -----   Public method EndOfEvent   -----------------------------------------
void R3BPsp::EndOfEvent()
{
    if (fVerboseLevel)
        Print();
    fPspCollection->Clear();

    ResetParameters();
}
// ----------------------------------------------------------------------------

// -----   Public method Register   -------------------------------------------
void R3BPsp::Register() { FairRootManager::Instance()->Register("PSPPoint", GetName(), fPspCollection, kTRUE); }
// ----------------------------------------------------------------------------

// -----   Public method GetCollection   --------------------------------------
TClonesArray* R3BPsp::GetCollection(Int_t iColl) const
{
    if (iColl == 0)
        return fPspCollection;
    else
        return NULL;
}
// ----------------------------------------------------------------------------

// -----   Public method Print   ----------------------------------------------
void R3BPsp::Print(Option_t* option) const
{
    Int_t nHits = fPspCollection->GetEntriesFast();
    LOG(info) << "R3BPsp: " << nHits << " points registered in this event";
}
// ----------------------------------------------------------------------------

// -----   Public method Reset   ----------------------------------------------
void R3BPsp::Reset()
{
    fPspCollection->Clear();
    ResetParameters();
}
// ----------------------------------------------------------------------------

// -----   Public method CopyClones   -----------------------------------------
void R3BPsp::CopyClones(TClonesArray* cl1, TClonesArray* cl2, Int_t offset)
{
    Int_t nEntries = cl1->GetEntriesFast();
    LOG(info) << "R3BPsp: " << nEntries << " entries to add";
    TClonesArray& clref = *cl2;
    R3BPspPoint* oldpoint = NULL;
    for (Int_t i = 0; i < nEntries; i++)
    {
        oldpoint = dynamic_cast<R3BPspPoint*>(cl1->At(i));
        Int_t index = oldpoint->GetTrackID() + offset;
        oldpoint->SetTrackID(index);
        new (clref[fPosIndex]) R3BPspPoint(*oldpoint);
        fPosIndex++;
    }
    LOG(info) << "R3BPsp: " << cl2->GetEntriesFast() << " merged entries";
}

// -----   Private method AddHit   --------------------------------------------
R3BPspPoint* R3BPsp::AddHit(Int_t trackID,
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
    TClonesArray& clref = *fPspCollection;
    Int_t size = clref.GetEntriesFast();
    if (fVerboseLevel > 1)
    {
        LOG(info) << "R3BPsp: Adding Point at (" << posIn.X() << ", " << posIn.Y() << ", " << posIn.Z()
                  << ") cm,  detector " << detID << ", track " << trackID << ", energy loss " << eLoss * 1e06 << " keV";
    }
    return new (clref[size]) R3BPspPoint(trackID, detID, plane, posIn, posOut, momIn, momOut, time, length, eLoss);
}

// -----   Public method ConstructGeometry   ----------------------------------
void R3BPsp::ConstructGeometry()
{
    TString fileName = GetGeometryFileName();
    if (fileName.EndsWith(".root"))
    {
        LOG(info) << "Constructing PSP geometry from ROOT file " << fileName.Data();
        ConstructRootGeometry();

        // TODO: Check if this works as expected
        if (!fCombiTrans.IsIdentity())
        {
            TGeoNode* n = gGeoManager->GetTopNode()->GetDaughter(gGeoManager->GetTopNode()->GetNdaughters() - 1);
            TGeoCombiTrans* combtrans = dynamic_cast<TGeoCombiTrans*>((dynamic_cast<TGeoNodeMatrix*>(n))->GetMatrix());

            *combtrans = fCombiTrans;
        }

        TGeoNode* psp_node = gGeoManager->GetTopVolume()->GetNode("PSPWorld_0");

        TGeoNode* node = psp_node->GetVolume()->GetNode("PSP1LogWorld_1");
        TGeoCombiTrans* combtrans = dynamic_cast<TGeoCombiTrans*>((dynamic_cast<TGeoNodeMatrix*>(node))->GetMatrix());
        combtrans->SetDz(fZ1);

        node = psp_node->GetVolume()->GetNode("PSP2LogWorld_2");
        combtrans = dynamic_cast<TGeoCombiTrans*>((dynamic_cast<TGeoNodeMatrix*>(node))->GetMatrix());
        combtrans->SetDz(fZ2);

        node = psp_node->GetVolume()->GetNode("PSP3LogWorld_3");
        combtrans = dynamic_cast<TGeoCombiTrans*>((dynamic_cast<TGeoNodeMatrix*>(node))->GetMatrix());
        combtrans->SetDz(fZ3);
    }
    else
    {
        LOG(fatal) << "PSP geometry file is not specified";
        exit(1);
    }
}

Bool_t R3BPsp::CheckIfSensitive(std::string name)
{
    if (TString(name).Contains("PSP1Log") || TString(name).Contains("PSP2Log") || TString(name).Contains("PSP3Log"))
    {
        return kTRUE;
    }
    return kFALSE;
}

ClassImp(R3BPsp)
