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

#include "R3BmTof.h"
#include "FairGeoInterface.h"
#include "FairGeoLoader.h"
#include "FairGeoNode.h"
#include "FairGeoRootBuilder.h"
#include "FairRootManager.h"
#include "FairRun.h"
#include "FairRuntimeDb.h"
#include "FairVolume.h"
#include "R3BGeomTof.h"
#include "R3BGeomTofPar.h"
#include "R3BMCStack.h"
#include "R3BTGeoPar.h"
#include "R3BmTofPoint.h"
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

R3BmTof::R3BmTof()
    : R3BmTof("")
{
}

R3BmTof::R3BmTof(const TString& geoFile, const TGeoTranslation& trans, const TGeoRotation& rot)
    : R3BmTof(geoFile, { trans, rot })
{
}

R3BmTof::R3BmTof(const TString& geoFile, const TGeoCombiTrans& combi)
    : R3BDetector("R3BmTof", kMTOF, geoFile, combi)
    , fmTofCollection(new TClonesArray("R3BmTofPoint"))
    , fPosIndex(0)
    , kGeoSaved(kFALSE)
    , flGeoPar(new TList())
{
    flGeoPar->SetName(GetName());
    ResetParameters();
}

R3BmTof::~R3BmTof()
{
    if (flGeoPar)
    {
        delete flGeoPar;
    }
    if (fmTofCollection)
    {
        fmTofCollection->Delete();
        delete fmTofCollection;
    }
}

void R3BmTof::Initialize()
{
    FairDetector::Initialize();

    LOG(info) << "R3BmTof: initialisation";
    LOG(debug) << "R3BmTof: Sci. Vol. (McId) " << gMC->VolId("mTOFLog");

    fTGeoPar = dynamic_cast<R3BTGeoPar*>(FairRuntimeDb::instance()->getContainer("mTofGeoPar"));

    // Position and rotation
    TGeoNode* main_vol = gGeoManager->GetTopVolume()->FindNode("mTOF_0");
    TGeoMatrix* matr = main_vol->GetMatrix();
    fTGeoPar->SetPosXYZ(matr->GetTranslation()[0], matr->GetTranslation()[1], matr->GetTranslation()[2]);
    fTGeoPar->SetRotXYZ(0., -TMath::Abs(TMath::ASin(matr->GetRotationMatrix()[2]) * TMath::RadToDeg()), 0.);

    // Dimensions
    TGeoBBox* box = dynamic_cast<TGeoBBox*>(gGeoManager->GetVolume(gMC->VolId("mTOFLog"))->GetShape());
    fTGeoPar->SetDimXYZ(box->GetDX(), box->GetDY(), box->GetDZ());

    fTGeoPar->setChanged();
}

void R3BmTof::SetSpecialPhysicsCuts()
{
    LOG(info) << "-I- R3BmTof: Adding customized Physics cut ... ";

    if (gGeoManager)
    {
        TGeoMedium* pSi = gGeoManager->GetMedium("plasticFormTOF");
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

            LOG(info) << "-I- R3BmTof: plasticFormTOF Medium Id " << pSi->GetId() << " Energy Cut-Off : " << cutE
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
Bool_t R3BmTof::ProcessHits(FairVolume* vol)
{

    // Simple Det plane

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

    // Set additional parameters at exit of active volume. Create R3BmTofPoint.
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

        AddHit(fTrackID,
               fVolumeID,
               TVector3(fPosIn.X(), fPosIn.Y(), fPosIn.Z()),
               TVector3(fPosOut.X(), fPosOut.Y(), fPosOut.Z()),
               TVector3(fMomIn.Px(), fMomIn.Py(), fMomIn.Pz()),
               TVector3(fMomOut.Px(), fMomOut.Py(), fMomOut.Pz()),
               fTime,
               fLength,
               fELoss);

        // Increment number of mTofPoints for this track
        R3BStack* stack = dynamic_cast<R3BStack*>(gMC->GetStack());
        stack->AddPoint(kMTOF);

        ResetParameters();
    }

    return kTRUE;
}

// ----------------------------------------------------------------------------
// void R3BmTof::SaveGeoParams(){
//
//  cout << " -I Save STS geo params " << endl;
//
//  TFolder *mf = (TFolder*) gDirectory->FindObjectAny("cbmroot");
//  cout << " mf: " << mf << endl;
//  TFolder *stsf = NULL;
//  if (mf ) stsf = (TFolder*) mf->FindObjectAny(GetName());
//  cout << " stsf: " << stsf << endl;
//  if (stsf) stsf->Add( flGeoPar0 ) ;
//  FairRootManager::Instance()->WriteFolder();
//  mf->Write("cbmroot",TObject::kWriteDelete);
//}

// -----   Public method EndOfEvent   -----------------------------------------
void R3BmTof::BeginEvent()
{
    //  if (! kGeoSaved ) {
    //      SaveGeoParams();
    //  cout << "-I STS geometry parameters saved " << endl;
    //  kGeoSaved = kTRUE;
    //  }
}

// -----   Public method EndOfEvent   -----------------------------------------
void R3BmTof::EndOfEvent()
{
    if (fVerboseLevel)
        Print();
    fmTofCollection->Clear();

    ResetParameters();
}
// ----------------------------------------------------------------------------

// -----   Public method Register   -------------------------------------------
void R3BmTof::Register() { FairRootManager::Instance()->Register("mTOFPoint", GetName(), fmTofCollection, kTRUE); }
// ----------------------------------------------------------------------------

// -----   Public method GetCollection   --------------------------------------
TClonesArray* R3BmTof::GetCollection(Int_t iColl) const
{
    if (iColl == 0)
        return fmTofCollection;
    else
        return NULL;
}
// ----------------------------------------------------------------------------

// -----   Public method Print   ----------------------------------------------
void R3BmTof::Print(Option_t* option) const
{
    Int_t nHits = fmTofCollection->GetEntriesFast();
    LOG(info) << "R3BmTof: " << nHits << " points registered in this event";
}
// ----------------------------------------------------------------------------

// -----   Public method Reset   ----------------------------------------------
void R3BmTof::Reset()
{
    fmTofCollection->Clear();
    ResetParameters();
}
// ----------------------------------------------------------------------------

// -----   Public method CopyClones   -----------------------------------------
void R3BmTof::CopyClones(TClonesArray* cl1, TClonesArray* cl2, Int_t offset)
{
    Int_t nEntries = cl1->GetEntriesFast();
    LOG(info) << "R3BmTof: " << nEntries << " entries to add";
    TClonesArray& clref = *cl2;
    R3BmTofPoint* oldpoint = NULL;
    for (Int_t i = 0; i < nEntries; i++)
    {
        oldpoint = dynamic_cast<R3BmTofPoint*>(cl1->At(i));
        Int_t index = oldpoint->GetTrackID() + offset;
        oldpoint->SetTrackID(index);
        new (clref[fPosIndex]) R3BmTofPoint(*oldpoint);
        fPosIndex++;
    }
    LOG(info) << "R3BmTof: " << cl2->GetEntriesFast() << " merged entries";
}

// -----   Private method AddHit   --------------------------------------------
R3BmTofPoint* R3BmTof::AddHit(Int_t trackID,
                              Int_t detID,
                              TVector3 posIn,
                              TVector3 posOut,
                              TVector3 momIn,
                              TVector3 momOut,
                              Double_t time,
                              Double_t length,
                              Double_t eLoss)
{
    TClonesArray& clref = *fmTofCollection;
    Int_t size = clref.GetEntriesFast();
    if (fVerboseLevel > 1)
        LOG(info) << "R3BmTof: Adding Point at (" << posIn.X() << ", " << posIn.Y() << ", " << posIn.Z()
                  << ") cm,  detector " << detID << ", track " << trackID << ", energy loss " << eLoss * 1e06 << " keV";
    return new (clref[size]) R3BmTofPoint(trackID, detID, posIn, posOut, momIn, momOut, time, length, eLoss);
}

Bool_t R3BmTof::CheckIfSensitive(std::string name)
{
    if (TString(name).Contains("mTOFLog"))
    {
        return kTRUE;
    }
    return kFALSE;
}

ClassImp(R3BmTof)
