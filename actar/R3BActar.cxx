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

#include "R3BActar.h"
#include "FairGeoInterface.h"
#include "FairGeoLoader.h"
#include "FairGeoNode.h"
#include "FairGeoRootBuilder.h"
#include "FairRootManager.h"
#include "FairRun.h"
#include "FairRuntimeDb.h"
#include "FairVolume.h"
#include "GarfieldPhysics.h"
#include "R3BActarPoint.h"
#include "R3BMCStack.h"
#include "TClonesArray.h"
#include "TObjArray.h"
#include "TParticle.h"
#include "TVirtualMC.h"
#include <stdlib.h>

R3BActar::R3BActar()
    : R3BActar("")
{
}

R3BActar::R3BActar(const TString& geoFile, const TGeoTranslation& trans, const TGeoRotation& rot)
    : R3BActar(geoFile, { trans, rot })
{
}

R3BActar::R3BActar(const TString& geoFile, const TGeoCombiTrans& combi)
    : R3BDetector("R3BActar", kACTAR, geoFile, combi)
    , fPspCollection(new TClonesArray("R3BActarPoint"))
    , fPosIndex(0)
    , kGeoSaved(kFALSE)
    , flGeoPar(new TList())
{
    flGeoPar->SetName(GetName());
    ResetParameters();
}

R3BActar::~R3BActar()
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

void R3BActar::Initialize()
{
    FairDetector::Initialize();

    LOG(info) << "R3BActar: initialisation";
    LOG(debug) << "R3BActar: Vol. (McId) " << gMC->VolId("PSP1Log");
}

void R3BActar::SetSpecialPhysicsCuts()
{
    // LOG(info) << "R3BActar: Adding customized Physics cut ... ";

    return;

    if (gGeoManager)
    {
        TGeoMedium* pHe = gGeoManager->GetMedium("helium");
        if (pHe)
        {
            // Setting processes for Si only
            gMC->Gstpar(pHe->GetId(), "LOSS", 3);
            gMC->Gstpar(pHe->GetId(), "STRA", 1.0);
            gMC->Gstpar(pHe->GetId(), "PAIR", 1.0);
            gMC->Gstpar(pHe->GetId(), "COMP", 1.0);
            gMC->Gstpar(pHe->GetId(), "PHOT", 1.0);
            gMC->Gstpar(pHe->GetId(), "ANNI", 1.0);
            gMC->Gstpar(pHe->GetId(), "BREM", 1.0);
            gMC->Gstpar(pHe->GetId(), "HADR", 1.0);
            gMC->Gstpar(pHe->GetId(), "DRAY", 1.0);
            gMC->Gstpar(pHe->GetId(), "DCAY", 1.0);
            gMC->Gstpar(pHe->GetId(), "MULS", 1.0);
            gMC->Gstpar(pHe->GetId(), "RAYL", 1.0);

            // Setting Energy-CutOff for Si Only
            Double_t cutE = fCutE; // GeV-> 1 keV

            LOG(info) << "R3BActar: helium Medium Id " << pHe->GetId() << " Energy Cut-Off : " << cutE << " GeV";

            // Si
            gMC->Gstpar(pHe->GetId(), "CUTGAM", cutE); /** gammas (GeV)*/
            gMC->Gstpar(pHe->GetId(), "CUTELE", cutE); /** electrons (GeV)*/
            gMC->Gstpar(pHe->GetId(), "CUTNEU", cutE); /** neutral hadrons (GeV)*/
            gMC->Gstpar(pHe->GetId(), "CUTHAD", cutE); /** charged hadrons (GeV)*/
            gMC->Gstpar(pHe->GetId(), "CUTMUO", cutE); /** muons (GeV)*/
            gMC->Gstpar(pHe->GetId(), "BCUTE", cutE);  /** electron bremsstrahlung (GeV)*/
            gMC->Gstpar(pHe->GetId(), "BCUTM", cutE);  /** muon and hadron bremsstrahlung(GeV)*/
            gMC->Gstpar(pHe->GetId(), "DCUTE", cutE);  /** delta-rays by electrons (GeV)*/
            gMC->Gstpar(pHe->GetId(), "DCUTM", cutE);  /** delta-rays by muons (GeV)*/
            gMC->Gstpar(pHe->GetId(), "PPCUTM", -1.);  /** direct pair production by muons (GeV)*/
        }
    } //! gGeoManager
}

// -----   Public method ProcessHits  --------------------------------------
Bool_t R3BActar::ProcessHits(FairVolume* vol)
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

    // Set additional parameters at exit of active volume. Create R3BActarPoint.
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
        R3BStack* stack = (R3BStack*)gMC->GetStack();
        stack->AddPoint(kACTAR);

        ResetParameters();
    }

    return kTRUE;
}

// -----   Public method EndOfEvent   -----------------------------------------
void R3BActar::BeginEvent() {}

// -----   Public method EndOfEvent   -----------------------------------------
void R3BActar::EndOfEvent()
{
    /// Update the collected hit information from Garfield interface
    GarfieldPhysics* garfieldPhysics = GarfieldPhysics::GetInstance();

    // get energy deposit from Garfield, convert it in GeV
    Double_t edep = garfieldPhysics->GetEnergyDeposit_MeV() * 1e-03;
    Double_t avalancheSize = garfieldPhysics->GetAvalancheSize();
    Double_t gain = garfieldPhysics->GetGain();

    LOG(info) << "Garfield output:";
    LOG(info) << "  energy deposit: " << edep;
    LOG(info) << "  avalanche size: " << avalancheSize;
    LOG(info) << "            gain: " << gain;

    if (fVerboseLevel)
        Print();
    fPspCollection->Clear();

    ResetParameters();

    garfieldPhysics->Clear();
}
// ----------------------------------------------------------------------------

// -----   Public method Register   -------------------------------------------
void R3BActar::Register() { FairRootManager::Instance()->Register("ActarPoint", GetName(), fPspCollection, kTRUE); }
// ----------------------------------------------------------------------------

// -----   Public method GetCollection   --------------------------------------
TClonesArray* R3BActar::GetCollection(Int_t iColl) const
{
    if (iColl == 0)
        return fPspCollection;
    else
        return NULL;
}
// ----------------------------------------------------------------------------

// -----   Public method Print   ----------------------------------------------
void R3BActar::Print(Option_t* option) const
{
    Int_t nHits = fPspCollection->GetEntriesFast();
    LOG(info) << "R3BActar: " << nHits << " points registered in this event";
}
// ----------------------------------------------------------------------------

// -----   Public method Reset   ----------------------------------------------
void R3BActar::Reset()
{
    fPspCollection->Clear();
    ResetParameters();
}
// ----------------------------------------------------------------------------

// -----   Public method CopyClones   -----------------------------------------
void R3BActar::CopyClones(TClonesArray* cl1, TClonesArray* cl2, Int_t offset)
{
    Int_t nEntries = cl1->GetEntriesFast();
    LOG(info) << "R3BActar: " << nEntries << " entries to add";
    TClonesArray& clref = *cl2;
    R3BActarPoint* oldpoint = NULL;
    for (Int_t i = 0; i < nEntries; i++)
    {
        oldpoint = (R3BActarPoint*)cl1->At(i);
        Int_t index = oldpoint->GetTrackID() + offset;
        oldpoint->SetTrackID(index);
        new (clref[fPosIndex]) R3BActarPoint(*oldpoint);
        fPosIndex++;
    }
    LOG(info) << "R3BActar: " << cl2->GetEntriesFast() << " merged entries";
}

// -----   Private method AddHit   --------------------------------------------
R3BActarPoint* R3BActar::AddHit(Int_t trackID,
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
        LOG(info) << "R3BActar: Adding Point at (" << posIn.X() << ", " << posIn.Y() << ", " << posIn.Z()
                  << ") cm,  detector " << detID << ", track " << trackID << ", energy loss " << eLoss * 1e06 << " keV";
    }
    return new (clref[size]) R3BActarPoint(trackID, detID, plane, posIn, posOut, momIn, momOut, time, length, eLoss);
}

Bool_t R3BActar::CheckIfSensitive(std::string name)
{
    if (TString(name).Contains("Actar_chamber_inner"))
    {
        return kTRUE;
    }
    return kFALSE;
}

ClassImp(R3BActar)
