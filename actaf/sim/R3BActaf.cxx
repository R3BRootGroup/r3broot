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

#include <FairRootManager.h>
#include <FairVolume.h>

#include <TClonesArray.h>
#include <TVirtualMC.h>
#include <stdlib.h>

#include "R3BActaf.h"
#include "R3BActafPoint.h"
#include "R3BLogger.h"
#include "R3BMCStack.h"

R3BActaf::R3BActaf()
    : R3BActaf("")
{
}

R3BActaf::R3BActaf(const TString& geoFile, const TGeoTranslation& trans, const TGeoRotation& rot)
    : R3BActaf(geoFile, { trans, rot })
{
}

R3BActaf::R3BActaf(const TString& geoFile, const TGeoCombiTrans& combi)
    : R3BDetector("R3BActaf", kACTAF, geoFile, combi)
    , fActafCollection(new TClonesArray("R3BActafPoint"))
{
    ResetParameters();
}

R3BActaf::~R3BActaf()
{
    if (fActafCollection)
    {
        fActafCollection->Delete();
        delete fActafCollection;
    }
}

void R3BActaf::Initialize()
{
    FairDetector::Initialize();

    R3BLOG(info, " ");
    // LOG(debug) << "R3BActaf: Vol. (McId) " << gMC->VolId("PSP1Log");
}

// -----   Public method ProcessHits  --------------------------------------
Bool_t R3BActaf::ProcessHits(FairVolume* vol)
{
    // 2 Simple Det PLane
    // get Info from DCH planes
    Int_t copyNo = -1;
    Int_t planeNr = -1;
    // Get the Geo info from MC Point
    TVirtualMC::GetMC()->CurrentVolID(copyNo);
    TVirtualMC::GetMC()->CurrentVolOffID(1, planeNr);

    if (TVirtualMC::GetMC()->IsTrackEntering())
    {
        fELoss = 0.;
        fTime_in = TVirtualMC::GetMC()->TrackTime() * 1.0e09;
        fLength_in = TVirtualMC::GetMC()->TrackLength();
        TVirtualMC::GetMC()->TrackPosition(fPosIn);
        TVirtualMC::GetMC()->TrackMomentum(fMomIn);
    }

    // Sum energy loss for all steps in the active volume
    fELoss += TVirtualMC::GetMC()->Edep();

    // Set additional parameters at exit of active volume. Create R3BActafPoint.
    if (TVirtualMC::GetMC()->IsTrackExiting() || TVirtualMC::GetMC()->IsTrackStop() ||
        TVirtualMC::GetMC()->IsTrackDisappeared())
    {
        fTrackID = gMC->GetStack()->GetCurrentTrackNumber();
        fVolumeID = vol->getMCid();
        gMC->TrackPosition(fPosOut);
        gMC->TrackMomentum(fMomOut);
        if (fELoss == 0.)
        {
            return kFALSE;
        }

        fTime_out =
            TVirtualMC::GetMC()->TrackTime() * 1.0e09; // also in case particle is stopped in detector, or decays...
        fLength_out = TVirtualMC::GetMC()->TrackLength();
        fTime = (fTime_out + fTime_in) / 2.;
        fLength = (fLength_out + fLength_in) / 2.;

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
        auto* stack = static_cast<R3BStack*>(TVirtualMC::GetMC()->GetStack());
        stack->AddPoint(kACTAF);
        ResetParameters();
    }

    return kTRUE;
}

// -----   Public method EndOfEvent   -----------------------------------------
void R3BActaf::EndOfEvent()
{
    if (fVerboseLevel)
    {
        Print();
    }
    Reset();
}

// -----   Public method Register   -------------------------------------------
void R3BActaf::Register() { FairRootManager::Instance()->Register("ActafPoint", GetName(), fActafCollection, kTRUE); }

// -----   Public method GetCollection   --------------------------------------
TClonesArray* R3BActaf::GetCollection(Int_t iColl) const
{
    if (iColl == 0)
    {
        return fActafCollection;
    }
    else
    {
        return nullptr;
    }
}

// -----   Public method Print   ----------------------------------------------
void R3BActaf::Print(Option_t* option) const
{
    Int_t nHits = fActafCollection->GetEntriesFast();
    R3BLOG(info, nHits << " points registered in this event");
}

// -----   Public method Reset   ----------------------------------------------
void R3BActaf::Reset()
{
    fActafCollection->Clear();
    ResetParameters();
}

// -----   Private method AddHit   --------------------------------------------
R3BActafPoint* R3BActaf::AddHit(Int_t trackID,
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
    TClonesArray& clref = *fActafCollection;
    Int_t size = clref.GetEntriesFast();
    return new (clref[size]) R3BActafPoint(trackID, detID, plane, posIn, posOut, momIn, momOut, time, length, eLoss);
}

Bool_t R3BActaf::CheckIfSensitive(std::string name)
{
    if (TString(name).Contains("Actar_chamber_inner"))
    {
        return kTRUE;
    }
    return kFALSE;
}

ClassImp(R3BActaf)
