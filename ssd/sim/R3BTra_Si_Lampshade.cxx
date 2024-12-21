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

// -------------------------------------------------------------------------
// -----                        R3BTra source file                     -----
// -----                  Created 26/03/09  by D.Bertini               -----
// -------------------------------------------------------------------------
#include "R3BTra.h"

#include "R3BGeoTra.h"
#include "R3BGeoTraPar.h"
#include "R3BMCStack.h"
#include "R3BTraPoint.h"

#include "FairGeoInterface.h"
#include "FairGeoLoader.h"
#include "FairGeoNode.h"
#include "FairGeoRootBuilder.h"
#include "FairLogger.h"
#include "FairRootManager.h"
#include "FairRun.h"
#include "FairRuntimeDb.h"
#include "FairVolume.h"

#include "TClonesArray.h"
#include "TGeoMCGeometry.h"
#include "TObjArray.h"
#include "TParticle.h"
#include "TVirtualMC.h"

// includes for modeling
#include "TGeoArb8.h"
#include "TGeoBBox.h"
#include "TGeoBoolNode.h"
#include "TGeoCompositeShape.h"
#include "TGeoCone.h"
#include "TGeoManager.h"
#include "TGeoMaterial.h"
#include "TGeoMatrix.h"
#include "TGeoMedium.h"
#include "TGeoPara.h"
#include "TGeoPgon.h"
#include "TGeoShapeAssembly.h"
#include "TGeoSphere.h"
#include "TGeoTrd1.h"
#include "TGeoTube.h"
#include "TParticle.h"
#include "TVirtualMC.h"

// -----   Default constructor   -------------------------------------------
R3BTra::R3BTra()
    : R3BDetector("R3BTra", kTRUE, kTRA)
{
    ResetParameters();
    fTraCollection = new TClonesArray("R3BTraPoint");
    fPosIndex = 0;
    kGeoSaved = kFALSE;
    flGeoPar = new TList();
    flGeoPar->SetName(GetName());
    fVerboseLevel = 1;
}
// -------------------------------------------------------------------------

// -----   Standard constructor   ------------------------------------------
R3BTra::R3BTra(const char* name, Bool_t active)
    : R3BDetector(name, active, kTRA)
{
    ResetParameters();
    fTraCollection = new TClonesArray("R3BTraPoint");
    fPosIndex = 0;
    kGeoSaved = kFALSE;
    flGeoPar = new TList();
    flGeoPar->SetName(GetName());
    fVerboseLevel = 1;
}
// -------------------------------------------------------------------------

// -----   Destructor   ----------------------------------------------------
R3BTra::~R3BTra()
{

    if (flGeoPar)
        delete flGeoPar;
    if (fTraCollection)
    {
        fTraCollection->Delete();
        delete fTraCollection;
    }
}
// -------------------------------------------------------------------------

// ----   Initialize   -----------------------------------------------------

void R3BTra::Initialize()
{
    FairDetector::Initialize();

    LOG(info);
    LOG(info) << "-I- R3BTra: initialisation";
    LOG(info) << "-I- R3BTra: Sens. Vol. (McId) " << gMC->VolId("TraLog");
}

// -----   Public method ProcessHits  --------------------------------------
Bool_t R3BTra::ProcessHits(FairVolume* vol)
{

    //      LOG(info) << " -I process hit called for:" <<  vol->GetName();
    // Set parameters at entrance of volume. Reset ELoss.

    //    if ( vol ) {
    //        LOG(info) << " Name Id:copy "
    //            << vol->getName() << " : " << vol->getMCid() << " : " << vol->getCopyNo();
    //        Int_t copyNo=0;
    //        LOG(info) << " Geant: " << gMC->CurrentVolID(copyNo) << ":" << copyNo;
    //    }

    if (gMC->IsTrackEntering())
    {
        fELoss = 0.;
        fTime = gMC->TrackTime() * 1.0e09;
        fLength = gMC->TrackLength();
        gMC->TrackPosition(fPosIn);
        gMC->TrackMomentum(fMomIn);
        // LOG(info) << "X,Y,X tracker=" << fPosIn(0) << " " << fPosIn(1) << " " << fPosIn(2);
        // LOG(info) << "track length=" << fLength;
    }

    // Sum energy loss for all steps in the active volume
    fELoss += gMC->Edep();
    // LOG(info) << "Tracker Eloss=" << fELoss << "  " << gMC->Edep();

    // Set additional parameters at exit of active volume. Create R3BTraPoint.
    if (gMC->IsTrackExiting() || gMC->IsTrackStop() || gMC->IsTrackDisappeared())
    {
        fTrackID = gMC->GetStack()->GetCurrentTrackNumber();
        fVolumeID = vol->getMCid();    // getCopyNo(); //getMCid();
        fDetCopyID = vol->getCopyNo(); // added by Marc
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

            //       LOG(info) << "1st direction: " << olddirection[0] << "," << olddirection[1] << "," <<
            //       olddirection[2];

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

            if (fPosIn.Z() < 30. && newpos[2] > 30.02)
            {
                cerr << "2nd direction: " << olddirection[0] << "," << olddirection[1] << "," << olddirection[2]
                     << " with safety = " << safety;
                cerr << "oldpos = " << oldpos[0] << "," << oldpos[1] << "," << oldpos[2];
                cerr << "newpos = " << newpos[0] << "," << newpos[1] << "," << newpos[2];
            }

            fPosOut.SetX(newpos[0]);
            fPosOut.SetY(newpos[1]);
            fPosOut.SetZ(newpos[2]);
        }

        AddHit(fTrackID,
               fVolumeID,
               fDetCopyID, // fdetCopyID, added by Marc
               TVector3(fPosIn.X(), fPosIn.Y(), fPosIn.Z()),
               TVector3(fPosOut.X(), fPosOut.Y(), fPosOut.Z()),
               TVector3(fMomIn.Px(), fMomIn.Py(), fMomIn.Pz()),
               TVector3(fMomOut.Px(), fMomOut.Py(), fMomOut.Pz()),
               fTime,
               fLength,
               fELoss);

        // Increment number of TraPoints for this track
        R3BStack* stack = (R3BStack*)gMC->GetStack();
        stack->AddPoint(kTRA);

        ResetParameters();
    }

    return kTRUE;
}
// ----------------------------------------------------------------------------
// void R3BTra::SaveGeoParams(){
//
//  LOG(info) << " -I Save STS geo params ";
//
//  TFolder *mf = (TFolder*) gDirectory->FindObjectAny("cbmroot");
//  LOG(info) << " mf: " << mf;
//  TFolder *stsf = NULL;
//  if (mf ) stsf = (TFolder*) mf->FindObjectAny(GetName());
//  LOG(info) << " stsf: " << stsf;
//  if (stsf) stsf->Add( flGeoPar0 ) ;
//  FairRootManager::Instance()->WriteFolder();
//  mf->Write("cbmroot",TObject::kWriteDelete);
//}

// -----   Public method EndOfEvent   -----------------------------------------
void R3BTra::BeginEvent()
{

    //  LOG(info) << "-I- begin tracker event called ##########";

    if (gGeoManager)
    {
        TGeoVolume* vol = gGeoManager->FindVolumeFast("TraLog");

        if (vol)
        {
            //    LOG(info) << "id tracker serial number : " << vol->GetNumber();
        }
    }

    //  if (! kGeoSaved ) {
    //      SaveGeoParams();
    //  LOG(info) << "-I STS geometry parameters saved ";
    //  kGeoSaved = kTRUE;
    //  }
}
// -----   Public method EndOfEvent   -----------------------------------------
void R3BTra::EndOfEvent()
{

    if (fVerboseLevel)
        Print();
    fTraCollection->Clear();

    ResetParameters();
}
// ----------------------------------------------------------------------------

// -----   Public method Register   -------------------------------------------
void R3BTra::Register() { FairRootManager::Instance()->Register("TraPoint", GetName(), fTraCollection, kTRUE); }
// ----------------------------------------------------------------------------

// -----   Public method GetCollection   --------------------------------------
TClonesArray* R3BTra::GetCollection(Int_t iColl) const
{
    if (iColl == 0)
        return fTraCollection;
    else
        return NULL;
}
// ----------------------------------------------------------------------------

// -----   Public method Print   ----------------------------------------------
void R3BTra::Print() const
{
    Int_t nHits = fTraCollection->GetEntriesFast();
    LOG(info) << "-I- R3BTra: " << nHits << " points registered in this event.";
}
// ----------------------------------------------------------------------------

// -----   Public method Reset   ----------------------------------------------
void R3BTra::Reset()
{
    fTraCollection->Clear();
    ResetParameters();
}
// ----------------------------------------------------------------------------

// -----   Public method CopyClones   -----------------------------------------
void R3BTra::CopyClones(TClonesArray* cl1, TClonesArray* cl2, Int_t offset)
{
    Int_t nEntries = cl1->GetEntriesFast();
    LOG(info) << "-I- R3BTra: " << nEntries << " entries to add.";
    TClonesArray& clref = *cl2;
    R3BTraPoint* oldpoint = NULL;
    for (Int_t i = 0; i < nEntries; i++)
    {
        oldpoint = (R3BTraPoint*)cl1->At(i);
        Int_t index = oldpoint->GetTrackID() + offset;
        oldpoint->SetTrackID(index);
        new (clref[fPosIndex]) R3BTraPoint(*oldpoint);
        fPosIndex++;
    }
    LOG(info) << " -I- R3BTra: " << cl2->GetEntriesFast() << " merged entries.";
}

// -----   Private method AddHit   --------------------------------------------
R3BTraPoint* R3BTra::AddHit(Int_t trackID,
                            Int_t detID,
                            Int_t detCopyID,
                            TVector3 posIn, // Int_t detCopyID added by Marc
                            TVector3 posOut,
                            TVector3 momIn,
                            TVector3 momOut,
                            Double_t time,
                            Double_t length,
                            Double_t eLoss)
{
    TClonesArray& clref = *fTraCollection;
    Int_t size = clref.GetEntriesFast();
    if (fVerboseLevel > 1)
        LOG(info) << "-I- R3BTra: Adding Point at (" << posIn.X() << ", " << posIn.Y() << ", " << posIn.Z()
                  << ") cm,  detector " << detID << ", track " << trackID << ", energy loss " << eLoss * 1e06 << " keV";
    return new (clref[size]) R3BTraPoint(
        trackID, detID, detCopyID, posIn, posOut, momIn, momOut, time, length, eLoss); // detCopy added by Marc
}

// -----   Public method ConstructGeometry   ----------------------------------
void R3BTra::ConstructGeometry()
{

    // out-of-file geometry definition
    Double_t dx, dy, dz;
    // Double_t par[20];
    Double_t rmin, rmax;
    Double_t a;
    Double_t thx, phx, thy, phy, thz, phz;
    Double_t phi1, phi2, theta1, theta2;
    Double_t z, density, radl, absl, w;
    Double_t tx, ty, tz;
    Int_t nel, numed;

    Double_t PI = 3.141592654;
    Int_t NSide1 = 6;            // Nb of detector in Layer 1
    Double_t WidthMax1 = 7.945;  // Max width of detector (cm)
    Double_t WidthMin1 = 2.25;   // Max width of detector (cm)
    Double_t Thickness1 = 0.005; // Half thickness of detector (cm)
    Double_t Length1 = 19.03;    // length of detector (cm
    Double_t InclAng1 = 14.9;    // angle d'inclinaison with respect to z axis (deg)
    Double_t Rmin1 = 2.;         // beam clearance 2cm radius
    Double_t AngRangeMin1 = 7;   // Min theta angle covered (deg)
    Double_t AngTrap1 = atan((WidthMax1 / 2 - WidthMin1 / 2) / Length1);         // (rad)
    Double_t WidthHalf1 = WidthMax1 - (Length1 / cos(AngTrap1)) * sin(AngTrap1); // width of detector at Length/2

    Int_t NSide2 = 12;           // Nb of detector in Layer 1
    Double_t WidthMax2 = 10.4;   // Max width of detector (cm)
    Double_t WidthMin2 = 1.3;    // Max width of detector (cm)
    Double_t Thickness2 = 0.015; // half thickness of detector (cm)
    Double_t Length2 = 30.6;     // length of detector (cm
    Double_t InclAng2 = 33.7;    // angle d'inclinaison with respect to z axis (deg)
    Double_t Rmin2 = 2.5;        // beam clearance 3cm radius
    Double_t AngRangeMin2 = 7.;  // Min theta angle covered (deg)
    Double_t AngTrap2 = atan((WidthMax2 / 2 - WidthMin2 / 2) / Length2);         // (rad)
    Double_t WidthHalf2 = WidthMax2 - (Length2 / cos(AngTrap2)) * sin(AngTrap2); // width of detector at Length/2

    Int_t NSide3 = 12;           // Nb of detector in Layer 1
    Double_t WidthMax3 = 10.4;   // Max width of detector (cm)
    Double_t WidthMin3 = 1.3;    // Max width of detector (cm)
    Double_t Thickness3 = 0.015; // half thickness of detector (cm)
    Double_t Length3 = 30.6;     // length of detector (cm
    Double_t InclAng3 = 33.7;    // angle d'inclinaison with respect to z axis (deg)
    Double_t Rmin3 = 2.685;      // beam clearance 3cm radius
    Double_t AngRangeMin3 = 7;   // Min theta angle covered (deg)
    Double_t AngTrap3 = atan((WidthMax3 / 2 - WidthMin3 / 2) / Length2);         // (rad)
    Double_t WidthHalf3 = WidthMax3 - (Length3 / cos(AngTrap3)) * sin(AngTrap3); // width of detector at Length/2

    /****************************************************************************/
    // Material definition

    // Mixture: Air
    TGeoMedium* pMed2 = NULL;
    if (gGeoManager->GetMedium("Air"))
    {
        pMed2 = gGeoManager->GetMedium("Air");
    }
    else
    {
        nel = 2;
        density = 0.001290;
        TGeoMixture* pMat2 = new TGeoMixture("Air", nel, density);
        a = 14.006740;
        z = 7.000000;
        w = 0.700000; // N
        pMat2->DefineElement(0, a, z, w);
        a = 15.999400;
        z = 8.000000;
        w = 0.300000; // O
        pMat2->DefineElement(1, a, z, w);
        pMat2->SetIndex(1);
        // Medium: Air
        numed = 1; // medium number
        Double_t par[8];
        par[0] = 0.000000; // isvol
        par[1] = 0.000000; // ifield
        par[2] = 0.000000; // fieldm
        par[3] = 0.000000; // tmaxfd
        par[4] = 0.000000; // stemax
        par[5] = 0.000000; // deemax
        par[6] = 0.000100; // epsil
        par[7] = 0.000000; // stmin
        pMed2 = new TGeoMedium("Air", numed, pMat2, par);
    }

    // Material: Silicon
    TGeoMedium* pMed22 = NULL;
    if (gGeoManager->GetMedium("Silicon"))
    {
        pMed22 = gGeoManager->GetMedium("Silicon");
    }
    else
    {
        a = 28.090000;
        z = 14.000000;
        density = 2.330000;
        radl = 9.351106;
        absl = 456.628489;

        TGeoMaterial* pMat22 = new TGeoMaterial("Silicon", a, z, density);

        pMat22->SetIndex(21);
        // Medium: Silicon
        numed = 21; // medium number
        Double_t par[8];
        par[0] = 0.000000; // isvol
        par[1] = 0.000000; // ifield
        par[2] = 0.000000; // fieldm
        par[3] = 0.000000; // tmaxfd
        par[4] = 0.000000; // stemax
        par[5] = 0.000000; // deemax
        par[6] = 0.000100; // epsil
        par[7] = 0.000000; // stmin
        pMed22 = new TGeoMedium("Silicon", numed, pMat22, par);
    }

    // Material: Copper
    TGeoMedium* pMed25 = NULL;
    if (gGeoManager->GetMedium("Copper"))
    {
        pMed25 = gGeoManager->GetMedium("Copper");
    }
    else
    {
        a = 63.540000;
        z = 29.000000;
        density = 8.960000;
        radl = 1.435029;
        absl = 155.874854;
        TGeoMaterial* pMat25 = new TGeoMaterial("Copper", a, z, density, radl, absl);
        pMat25->SetIndex(24);
        // Medium: Copper
        numed = 24; // medium number
        Double_t par[8];
        par[0] = 0.000000; // isvol
        par[1] = 0.000000; // ifield
        par[2] = 0.000000; // fieldm
        par[3] = 0.000000; // tmaxfd
        par[4] = 0.000000; // stemax
        par[5] = 0.000000; // deemax
        par[6] = 0.000100; // epsil
        par[7] = 0.000000; // stmin
        pMed25 = new TGeoMedium("Copper", numed, pMat25, par);
    }

    // Material: Aluminum
    TGeoMedium* pMed21 = NULL;
    if (gGeoManager->GetMedium("Aluminum"))
    {
        pMed21 = gGeoManager->GetMedium("Aluminum");
    }
    else
    {
        a = 26.980000;
        z = 13.000000;
        density = 2.700000;
        radl = 8.875105;
        absl = 388.793113;
        TGeoMaterial* pMat21 = new TGeoMaterial("Aluminum", a, z, density, radl, absl);
        pMat21->SetIndex(20);
        // Medium: Aluminum
        numed = 20; // medium number
        Double_t par[8];
        par[0] = 0.000000; // isvol
        par[1] = 0.000000; // ifield
        par[2] = 0.000000; // fieldm
        par[3] = 0.000000; // tmaxfd
        par[4] = 0.000000; // stemax
        par[5] = 0.000000; // deemax
        par[6] = 0.000100; // epsil
        par[7] = 0.000000; // stmin
        pMed21 = new TGeoMedium("Aluminum", numed, pMat21, par);
    }

    // Mixture: CarbonFibre
    TGeoMedium* pMed3 = NULL;
    if (gGeoManager->GetMedium("CarbonFibre"))
    {
        pMed3 = gGeoManager->GetMedium("CarbonFibre");
    }
    else
    {
        nel = 3;
        density = 1.690000;
        TGeoMixture* pMat3 = new TGeoMixture("CarbonFibre", nel, density);
        a = 12.010700;
        z = 6.000000;
        w = 0.844907; // C
        pMat3->DefineElement(0, a, z, w);
        a = 1.007940;
        z = 1.000000;
        w = 0.042543; // H
        pMat3->DefineElement(1, a, z, w);
        a = 15.999400;
        z = 8.000000;
        w = 0.112550; // O
        pMat3->DefineElement(2, a, z, w);
        // Medium: CarbonFibre
        numed = 802; // medium number
        pMat3->SetIndex(numed);
        Double_t par[8];
        par[0] = 0.000000; // isvol
        par[1] = 0.000000; // ifield
        par[2] = 0.000000; // fieldm
        par[3] = 0.000000; // tmaxfd
        par[4] = 0.000000; // stemax
        par[5] = 0.000000; // deemax
        par[6] = 0.000100; // epsil
        par[7] = 0.000000; // stmin
        pMed3 = new TGeoMedium("CarbonFibre", numed, pMat3, par);
    }

    // Material: Gold
    TGeoMedium* pMed79 = NULL;
    if (gGeoManager->GetMedium("Gold"))
    {
        LOG(info) << "-I- TGeoManager: Gold Medium already defined ";
        pMed79 = gGeoManager->GetMedium("Gold");
    }
    else
    {
        // Material definition
        // Material: Gold
        a = 196.9665;
        z = 79.000000;
        density = 19.300000;
        TGeoMaterial* pMat79 = new TGeoMaterial("Gold", a, z, density);
        pMat79->SetIndex(600);
        // Medium: Gold
        numed = 79; // medium number
        Double_t par[8];
        par[0] = 0.000000; // isvol
        par[1] = 0.000000; // ifield
        par[2] = 0.000000; // fieldm
        par[3] = 0.000000; // tmaxfd
        par[4] = 0.000000; // stemax
        par[5] = 0.000000; // deemax
        par[6] = 0.000100; // epsil
        par[7] = 0.000000; // stmin
        pMed79 = new TGeoMedium("Gold", numed, pMat79, par);
    }

    // TRANSFORMATION MATRICES

    // Detectors

    // First layer

    // Combi transformation:
    dx = 0.000;
    dy = -((Length1 / 2 + (WidthHalf1 / 2) * tan(AngTrap1)) * sin(InclAng1 * PI / 180.) +
           Rmin1); // considering  real barycentre position
    // dy = -((Length1/2)*sin(InclAng1*PI/180.)+ Rmin1);                                // considering  intersection of
    // 2 medianes
    dz = Length1 * cos(InclAng1 * PI / 180.) / 2 -
         (Length1 * cos(InclAng1 * PI / 180.) - (Rmin1 / tan(AngRangeMin1 * PI / 180.)));
    // Rotation:
    thx = 90.000000;
    phx = 0.000000;
    thy = 90. + InclAng1;
    phy = 90.000000;
    thz = InclAng1;
    phz = 90.000000;
    TGeoRotation* pMatrix3 = new TGeoRotation("", thx, phx, thy, phy, thz, phz); // geant3 angles
    // thx = -InclAng1;
    // thy = 0.;
    // thz = 0;
    // TGeoRotation *pMatrix3 = new TGeoRotation("",thz,thx,thy); // Euler phi (rot/z), theta (rot/x), psi (rot/y)
    // angles
    TGeoCombiTrans* pMatrix2 = new TGeoCombiTrans("", dx, dy, dz, pMatrix3); // bottom inner detector

    // Combi transformation:
    dx = dx - (WidthHalf1 / 2) * (1 / (cos(AngTrap1) * cos(AngTrap1))) *
                  (1 + cos((360. / NSide1) * PI / 180.)); // considering  real barycentre position;
    // dx =  -(WidthHalf1/2)*(1 + cos((360./NSide1)*PI/180));                                     // considering
    // intersection of 2 medianes
    dy = dy + (WidthHalf1 / 2) * (1 / (cos(AngTrap1) * cos(AngTrap1))) *
                  sin((360. / NSide1) * PI / 180.); // considering  real barycentre position;
    // dy = dy +(WidthHalf1/2)*sin((360./NSide1)*PI/180);                                         // considering
    // intersection of 2 medianes
    // Rotation:
    thx = 90.000000;
    phx = 0.000000 - 360. / NSide1;
    thy = 90. + InclAng1;
    phy = 90.000000 - 360. / NSide1;
    thz = InclAng1;
    phz = 90.000000 - 360. / NSide1;
    TGeoRotation* pMatrix5 = new TGeoRotation("", thx, phx, thy, phy, thz, phz); // geant3 angles
    // thx = -InclAng1;
    // thy = 0.;
    // thz = -39.;
    // TGeoRotation *pMatrix5 = new TGeoRotation("",thz,thx,thy); // Euler phi (rot/z), theta (rot/x), psi (rot/y)
    // angles
    TGeoCombiTrans* pMatrix4 = new TGeoCombiTrans("", dx, dy, dz, pMatrix5);

    // Combi transformation:
    dx = dx - (WidthHalf1 / 2) * (1 / (cos(AngTrap1) * cos(AngTrap1))) *
                  (cos((360. / NSide1) * PI / 180.) +
                   cos(2 * (360. / NSide1) * PI / 180.)); // considering  real barycentre position;
    // dx =  dx -(WidthHalf1/2)*(cos((360./NSide1)*PI/180.) + cos(2*(360./NSide1)*3.14159/180));   // considering
    // intersection of 2 medianes
    dy = dy + (WidthHalf1 / 2) * (1 / (cos(AngTrap1) * cos(AngTrap1))) *
                  (sin((360. / NSide1) * PI / 180.) +
                   sin(2 * (360. / NSide1) * PI / 180.)); // considering  real barycentre position;
    // dy =  dy +(WidthHalf1/2)*(sin((360./NSide1)*PI/180.) + sin(2*(360./NSide1)*PI/180));        // considering
    // intersection of 2 medianes
    // Rotation:
    thx = 90.000000;
    phx = 0.000000 - (2 * 360. / NSide1);
    thy = 90. + InclAng1;
    phy = 90.000000 - (2 * 360. / NSide1);
    thz = InclAng1;
    phz = 90.000000 - (2 * 360. / NSide1);
    TGeoRotation* pMatrix7 = new TGeoRotation("", thx, phx, thy, phy, thz, phz);
    TGeoCombiTrans* pMatrix6 = new TGeoCombiTrans("", dx, dy, dz, pMatrix7);

    // Combi transformation:
    dx = dx - (WidthHalf1 / 2) * (1 / (cos(AngTrap1) * cos(AngTrap1))) *
                  (cos(2 * (360. / NSide1) * PI / 180.) +
                   cos(3 * (360. / NSide1) * PI / 180.)); // considering  real barycentre position;
    // dx =  dx -(WidthHalf1/2)*(cos(2*(360./NSide1)*PI/180.) + cos(3*(360./NSide1)*3.14159/180));  // considering
    // intersection of 2 medianes
    dy = dy + (WidthHalf1 / 2) * (1 / (cos(AngTrap1) * cos(AngTrap1))) *
                  (sin(2 * (360. / NSide1) * PI / 180.) +
                   sin(3 * (360. / NSide1) * PI / 180.)); // considering  real barycentre position;
    // dy =  dy +(WidthHalf1/2)*(sin(2*(360./NSide1)*PI/180.) + sin(3*(360./NSide1)*PI/180));       // considering
    // intersection of 2 medianes
    // Rotation:
    thx = 90.000000;
    phx = 0.000000 - (3 * 360. / NSide1);
    thy = 90. + InclAng1;
    phy = 90.000000 - (3 * 360. / NSide1);
    thz = InclAng1;
    phz = 90.000000 - (3 * 360. / NSide1);
    TGeoRotation* pMatrix9 = new TGeoRotation("", thx, phx, thy, phy, thz, phz);
    TGeoCombiTrans* pMatrix8 = new TGeoCombiTrans("", dx, dy, dz, pMatrix9);

    // Combi transformation:
    dx = dx - (WidthHalf1 / 2) * (1 / (cos(AngTrap1) * cos(AngTrap1))) *
                  (cos(3 * (360. / NSide1) * PI / 180.) +
                   cos(4 * (360. / NSide1) * PI / 180.)); // considering  real barycentre position;
    // dx =  dx -(WidthHalf1/2)*(cos(3*(360./NSide1)*PI/180.) + cos(4*(360./NSide1)*3.14159/180));  // considering
    // intersection of 2 medianes
    dy = dy + (WidthHalf1 / 2) * (1 / (cos(AngTrap1) * cos(AngTrap1))) *
                  (sin(3 * (360. / NSide1) * PI / 180.) +
                   sin(4 * (360. / NSide1) * PI / 180.)); // considering  real barycentre position;
    // dy =  dy +(WidthHalf1/2)*(sin(3*(360./NSide1)*PI/180.) + sin(4*(360./NSide1)*PI/180));       // considering
    // intersection of 2 medianes
    // Rotation:
    thx = 90.000000;
    phx = 0.000000 - (4 * 360. / NSide1);
    thy = 90. + InclAng1;
    phy = 90.000000 - (4 * 360. / NSide1);
    thz = InclAng1;
    phz = 90.000000 - (4 * 360. / NSide1);
    TGeoRotation* pMatrix11 = new TGeoRotation("", thx, phx, thy, phy, thz, phz);
    TGeoCombiTrans* pMatrix10 = new TGeoCombiTrans("", dx, dy, dz, pMatrix11);

    // Combi transformation:
    dx = dx - (WidthHalf1 / 2) * (1 / (cos(AngTrap1) * cos(AngTrap1))) *
                  (cos(4 * (360. / NSide1) * PI / 180.) +
                   cos(5 * (360. / NSide1) * PI / 180.)); // considering  real barycentre position;
    // dx =  dx -(WidthHalf1/2)*(cos(4*(360./NSide1)*PI/180.) + cos(5*(360./NSide1)*3.14159/180)); // considering
    // intersection of 2 medianes
    dy = dy + (WidthHalf1 / 2) * (1 / (cos(AngTrap1) * cos(AngTrap1))) *
                  (sin(4 * (360. / NSide1) * PI / 180.) +
                   sin(5 * (360. / NSide1) * PI / 180.)); // considering  real barycentre position;
    // dy =  dy +(WidthHalf1/2)*(sin(4*(360./NSide1)*PI/180.) + sin(5*(360./NSide1)*PI/180));      // considering
    // intersection of 2 medianes
    // Rotation:
    thx = 90.000000;
    phx = 0.000000 - (5 * 360. / NSide1);
    thy = 90. + InclAng1;
    phy = 90.000000 - (5 * 360. / NSide1);
    thz = InclAng1;
    phz = 90.000000 - (5 * 360. / NSide1);
    TGeoRotation* pMatrix13 = new TGeoRotation("", thx, phx, thy, phy, thz, phz);
    TGeoCombiTrans* pMatrix12 = new TGeoCombiTrans("", dx, dy, dz, pMatrix13);

    // Combi transformation:
    dx = dx - (WidthHalf1 / 2) * (1 / (cos(AngTrap1) * cos(AngTrap1))) *
                  (cos(5 * (360. / NSide1) * PI / 180.) +
                   cos(6 * (360. / NSide1) * PI / 180.)); // considering  real barycentre position;
    // dx =  dx -(WidthHalf1/2)*(cos(5*(360./NSide1)*PI/180.) + cos(6*(360./NSide1)*3.14159/180)); // considering
    // intersection of 2 medianes
    dy = dy + (WidthHalf1 / 2) * (1 / (cos(AngTrap1) * cos(AngTrap1))) *
                  (sin(5 * (360. / NSide1) * PI / 180.) +
                   sin(6 * (360. / NSide1) * PI / 180.)); // considering  real barycentre position;
    // dy =  dy +(WidthHalf1/2)*(sin(5*(360./NSide1)*PI/180.) + sin(6*(360./NSide1)*PI/180));      // considering
    // intersection of 2 medianes
    // Rotation:
    thx = 90.000000;
    phx = 0.000000 - (6 * 360. / NSide1);
    thy = 90. + InclAng1;
    phy = 90.000000 - (6 * 360. / NSide1);
    thz = InclAng1;
    phz = 90.000000 - (6 * 360. / NSide1);
    TGeoRotation* pMatrix15 = new TGeoRotation("", thx, phx, thy, phy, thz, phz);
    TGeoCombiTrans* pMatrix14 = new TGeoCombiTrans("", dx, dy, dz, pMatrix15);

    // Combi transformation:
    dx = dx - (WidthHalf1 / 2) * (1 / (cos(AngTrap1) * cos(AngTrap1))) *
                  (cos(6 * (360. / NSide1) * PI / 180.) +
                   cos(7 * (360. / NSide1) * PI / 180.)); // considering  real barycentre position;
    // dx =  dx -(WidthHalf1/2)*(cos(6*(360./NSide1)*PI/180. + cos(7*(360./NSide1)*3.14159/180));  // considering
    // intersection of 2 medianes
    dy = dy + (WidthHalf1 / 2) * (1 / (cos(AngTrap1) * cos(AngTrap1))) *
                  (sin(6 * (360. / NSide1) * PI / 180.) +
                   sin(7 * (360. / NSide1) * PI / 180.)); // considering  real barycentre position;
    // dy =  dy +(WidthHalf1/2)*(sin(6*(360./NSide1)*PI/180.) + sin(7*(360./NSide1)*PI/180));      // considering
    // intersection of 2 medianes
    // Rotation:
    thx = 90.000000;
    phx = 0.000000 - (7 * 360. / NSide1);
    thy = 90. + InclAng1;
    phy = 90.000000 - (7 * 360. / NSide1);
    thz = InclAng1;
    phz = 90.000000 - (7 * 360. / NSide1);
    TGeoRotation* pMatrix17 = new TGeoRotation("", thx, phx, thy, phy, thz, phz);
    TGeoCombiTrans* pMatrix16 = new TGeoCombiTrans("", dx, dy, dz, pMatrix17);

    // Combi transformation:
    dx = dx - (WidthHalf1 / 2) * (1 / (cos(AngTrap1) * cos(AngTrap1))) *
                  (cos(7 * (360. / NSide1) * PI / 180.) +
                   cos(8 * (360. / NSide1) * PI / 180.)); // considering  real barycentre position;
    // dx =  dx -(WidthHalf1/2)*(cos(7*(360./NSide1)*PI/180. + cos(8*(360./NSide1)*3.14159/180)); // considering
    // intersection of 2 medianes
    dy = dy + (WidthHalf1 / 2) * (1 / (cos(AngTrap1) * cos(AngTrap1))) *
                  (sin(7 * (360. / NSide1) * PI / 180.) +
                   sin(8 * (360. / NSide1) * PI / 180.)); // considering  real barycentre position;
    // dy =  dy +(WidthHalf1/2)*(sin(7*(360./NSide1)*PI/180.) + sin(8*(360./NSide1)*PI/180));     // considering
    // intersection of 2 medianes
    // Rotation:
    thx = 90.000000;
    phx = 0.000000 - (8 * 360. / NSide1);
    thy = 90. + InclAng1;
    phy = 90.000000 - (8 * 360. / NSide1);
    thz = InclAng1;
    phz = 90.000000 - (8 * 360. / NSide1);
    TGeoRotation* pMatrix19 = new TGeoRotation("", thx, phx, thy, phy, thz, phz);
    TGeoCombiTrans* pMatrix18 = new TGeoCombiTrans("", dx, dy, dz, pMatrix19);

    // Combi transformation:
    dx = dx - (WidthHalf1 / 2) * (1 / (cos(AngTrap1) * cos(AngTrap1))) *
                  (cos(8 * (360. / NSide1) * PI / 180.) +
                   cos(9 * (360. / NSide1) * PI / 180.)); // considering  real barycentre position;
    // dx =  dx -(WidthHalf1/2)*(cos(8*(360./NSide1)*PI/180. + cos(9*(360./NSide1)*3.14159/180)); // considering
    // intersection of 2 medianes
    dy = dy + (WidthHalf1 / 2) * (1 / (cos(AngTrap1) * cos(AngTrap1))) *
                  (sin(8 * (360. / NSide1) * PI / 180.) +
                   sin(9 * (360. / NSide1) * PI / 180.)); // considering  real barycentre position;
    // dy =  dy +(WidthHalf1/2)*(sin(8*(360./NSide1)*PI/180.) + sin(9*(360./NSide1)*PI/180));     // considering
    // intersection of 2 medianes
    // Rotation:
    thx = 90.000000;
    phx = 0.000000 - (9 * 360. / NSide1);
    thy = 90. + InclAng1;
    phy = 90.000000 - (9 * 360. / NSide1);
    thz = InclAng1;
    phz = 90.000000 - (9 * 360. / NSide1);
    TGeoRotation* pMatrix21 = new TGeoRotation("", thx, phx, thy, phy, thz, phz);
    TGeoCombiTrans* pMatrix20 = new TGeoCombiTrans("", dx, dy, dz, pMatrix21);

    // Second layer

    // Combi transformation:
    dx = 0.000;
    dy = -((Length2 / 2 + (WidthHalf2 / 2) * tan(AngTrap2)) * sin(InclAng2 * PI / 180.) +
           Rmin2); // considering  real barycentre position
    // dy = -((Length2/2)*sin(InclAng2*PI/180.)+ Rmin2);                               // considering  intersection of 2
    // medianes
    dz = Length2 * cos(InclAng2 * PI / 180.) / 2 -
         (Length2 * cos(InclAng2 * PI / 180.) - (Rmin2 / tan(AngRangeMin2 * PI / 180.)));
    // Rotation:
    thx = 90.000000;
    phx = 0.000000;
    thy = 90. + InclAng2;
    phy = 90.000000;
    thz = InclAng2;
    phz = 90.000000;
    TGeoRotation* pMatrix67 = new TGeoRotation("", thx, phx, thy, phy, thz, phz); // geant3 angles
    TGeoCombiTrans* pMatrix66 = new TGeoCombiTrans("", dx, dy, dz, pMatrix67);

    // Combi transformation:
    dx = dx - (WidthHalf2 / 2) * (1 / (cos(AngTrap2) * cos(AngTrap2))) *
                  (1 + cos((360. / NSide2) * PI / 180.)); // considering  real barycentre position;
    // dx =  -(WidthHalf2/2)*(1 + cos((360./NSide2)*PI/180));                                     // considering
    // intersection of 2 medianes
    dy = dy + (WidthHalf2 / 2) * (1 / (cos(AngTrap2) * cos(AngTrap2))) *
                  sin((360. / NSide2) * PI / 180.); // considering  real barycentre position;
    // dy = dy +(WidthHalf2/2)*sin((360./NSide2)*PI/180);                                         // considering
    // intersection of 2 medianes
    // Rotation:
    thx = 90.000000;
    phx = 0.000000 - 360. / NSide2;
    thy = 90. + InclAng2;
    phy = 90.000000 - 360. / NSide2;
    thz = InclAng2;
    phz = 90.000000 - 360. / NSide2;
    TGeoRotation* pMatrix69 = new TGeoRotation("", thx, phx, thy, phy, thz, phz);
    TGeoCombiTrans* pMatrix68 = new TGeoCombiTrans("", dx, dy, dz, pMatrix69);

    // Combi transformation:
    dx = dx - (WidthHalf2 / 2) * (1 / (cos(AngTrap2) * cos(AngTrap2))) *
                  (cos((360. / NSide2) * PI / 180.) +
                   cos(2 * (360. / NSide2) * PI / 180.)); // considering  real barycentre position;
    // dx =  dx -(WidthHalf2/2)*(cos((360./NSide2)*PI/180.) + cos(2*(360./NSide2)*3.14159/180));   // considering
    // intersection of 2 medianes
    dy = dy + (WidthHalf2 / 2) * (1 / (cos(AngTrap2) * cos(AngTrap2))) *
                  (sin((360. / NSide2) * PI / 180.) +
                   sin(2 * (360. / NSide2) * PI / 180.)); // considering  real barycentre position;
    // dy =  dy +(WidthHalf2/2)*(sin((360./NSide2)*PI/180.) + sin(2*(360./NSide2)*PI/180));        // considering
    // intersection of 2 medianes
    // Rotation:
    thx = 90.000000;
    phx = 0.000000 - (2 * 360. / NSide2);
    thy = 90. + InclAng2;
    phy = 90.000000 - (2 * 360. / NSide2);
    thz = InclAng2;
    phz = 90.000000 - (2 * 360. / NSide2);
    TGeoRotation* pMatrix71 = new TGeoRotation("", thx, phx, thy, phy, thz, phz);
    TGeoCombiTrans* pMatrix70 = new TGeoCombiTrans("", dx, dy, dz, pMatrix71);

    // Combi transformation:
    dx = dx - (WidthHalf2 / 2) * (1 / (cos(AngTrap2) * cos(AngTrap2))) *
                  (cos(2 * (360. / NSide2) * PI / 180.) +
                   cos(3 * (360. / NSide2) * PI / 180.)); // considering  real barycentre position;
    // dx =  dx -(WidthHalf2/2)*(cos(2*(360./NSide2)*PI/180.) + cos(3*(360./NSide2)*3.14159/180)); // considering
    // intersection of 2 medianes
    dy = dy + (WidthHalf2 / 2) * (1 / (cos(AngTrap2) * cos(AngTrap2))) *
                  (sin(2 * (360. / NSide2) * PI / 180.) +
                   sin(3 * (360. / NSide2) * PI / 180.)); // considering  real barycentre position;
    // dy =  dy +(WidthHalf2/2)*(sin(2*(360./NSide2)*PI/180.) + sin(3*(360./NSide2)*PI/180));      // considering
    // intersection of 2 medianes
    // Rotation:
    thx = 90.000000;
    phx = 0.000000 - (3 * 360. / NSide2);
    thy = 90. + InclAng2;
    phy = 90.000000 - (3 * 360. / NSide2);
    thz = InclAng2;
    phz = 90.000000 - (3 * 360. / NSide2);
    TGeoRotation* pMatrix73 = new TGeoRotation("", thx, phx, thy, phy, thz, phz);
    TGeoCombiTrans* pMatrix72 = new TGeoCombiTrans("", dx, dy, dz, pMatrix73);

    // Combi transformation:
    dx = dx - (WidthHalf2 / 2) * (1 / (cos(AngTrap2) * cos(AngTrap2))) *
                  (cos(3 * (360. / NSide2) * PI / 180.) +
                   cos(4 * (360. / NSide2) * PI / 180.)); // considering  real barycentre position;
    // dx =  dx -(WidthHalf2/2)*(cos(3*(360./NSide2)*PI/180.) + cos(4*(360./NSide2)*3.14159/180));  // considering
    // intersection of 2 medianes
    dy = dy + (WidthHalf2 / 2) * (1 / (cos(AngTrap2) * cos(AngTrap2))) *
                  (sin(3 * (360. / NSide2) * PI / 180.) +
                   sin(4 * (360. / NSide2) * PI / 180.)); // considering  real barycentre position;
    // dy =  dy +(WidthHalf2/2)*(sin(3*(360./NSide2)*PI/180.) + sin(4*(360./NSide2)*PI/180));       // considering
    // intersection of 2 medianes
    // Rotation:
    thx = 90.000000;
    phx = 0.000000 - (4 * 360. / NSide2);
    thy = 90. + InclAng2;
    phy = 90.000000 - (4 * 360. / NSide2);
    thz = InclAng2;
    phz = 90.000000 - (4 * 360. / NSide2);
    TGeoRotation* pMatrix75 = new TGeoRotation("", thx, phx, thy, phy, thz, phz);
    TGeoCombiTrans* pMatrix74 = new TGeoCombiTrans("", dx, dy, dz, pMatrix75);

    // Combi transformation:
    dx = dx - (WidthHalf2 / 2) * (1 / (cos(AngTrap2) * cos(AngTrap2))) *
                  (cos(4 * (360. / NSide2) * PI / 180.) +
                   cos(5 * (360. / NSide2) * PI / 180.)); // considering  real barycentre position;
    // dx =  dx -(WidthHalf2/2)*(cos(4*(360./NSide2)*PI/180.) + cos(5*(360./NSide2)*3.14159/180));  // considering
    // intersection of 2 medianes
    dy = dy + (WidthHalf2 / 2) * (1 / (cos(AngTrap2) * cos(AngTrap2))) *
                  (sin(4 * (360. / NSide2) * PI / 180.) +
                   sin(5 * (360. / NSide2) * PI / 180.)); // considering  real barycentre position;
    // dy =  dy +(WidthHalf2/2)*(sin(4*(360./NSide2)*PI/180.) + sin(5*(360./NSide2)*PI/180));       // considering
    // intersection of 2 medianes
    // Rotation:
    thx = 90.000000;
    phx = 0.000000 - (5 * 360. / NSide2);
    thy = 90. + InclAng2;
    phy = 90.000000 - (5 * 360. / NSide2);
    thz = InclAng2;
    phz = 90.000000 - (5 * 360. / NSide2);
    TGeoRotation* pMatrix77 = new TGeoRotation("", thx, phx, thy, phy, thz, phz);
    TGeoCombiTrans* pMatrix76 = new TGeoCombiTrans("", dx, dy, dz, pMatrix77);

    // Combi transformation:
    dx = dx - (WidthHalf2 / 2) * (1 / (cos(AngTrap2) * cos(AngTrap2))) *
                  (cos(5 * (360. / NSide2) * PI / 180.) +
                   cos(6 * (360. / NSide2) * PI / 180.)); // considering  real barycentre position;
    // dx =  dx -(WidthHalf2/2)*(cos(5*(360./NSide2)*PI/180.) + cos(6*(360./NSide2)*3.14159/180));  // considering
    // intersection of 2 medianes
    dy = dy + (WidthHalf2 / 2) * (1 / (cos(AngTrap2) * cos(AngTrap2))) *
                  (sin(5 * (360. / NSide2) * PI / 180.) +
                   sin(6 * (360. / NSide2) * PI / 180.)); // considering  real barycentre position;
    // dy =  dy +(WidthHalf2/2)*(sin(5*(360./NSide2)*PI/180.) + sin(6*(360./NSide2)*PI/180));       // considering
    // intersection of 2 medianes
    // Rotation:
    thx = 90.000000;
    phx = 0.000000 - (6 * 360. / NSide2);
    thy = 90. + InclAng2;
    phy = 90.000000 - (6 * 360. / NSide2);
    thz = InclAng2;
    phz = 90.000000 - (6 * 360. / NSide2);
    TGeoRotation* pMatrix79 = new TGeoRotation("", thx, phx, thy, phy, thz, phz);
    TGeoCombiTrans* pMatrix78 = new TGeoCombiTrans("", dx, dy, dz, pMatrix79);

    // Combi transformation:
    dx = dx - (WidthHalf2 / 2) * (1 / (cos(AngTrap2) * cos(AngTrap2))) *
                  (cos(6 * (360. / NSide2) * PI / 180.) +
                   cos(7 * (360. / NSide2) * PI / 180.)); // considering  real barycentre position;
    // dx =  dx -(WidthHalf2/2)*(cos(6*(360./NSide2)*PI/180.) + cos(7*(360./NSide2)*3.14159/180)); // considering
    // intersection of 2 medianes
    dy = dy + (WidthHalf2 / 2) * (1 / (cos(AngTrap2) * cos(AngTrap2))) *
                  (sin(6 * (360. / NSide2) * PI / 180.) +
                   sin(7 * (360. / NSide2) * PI / 180.)); // considering  real barycentre position;
    // dy =  dy +(WidthHalf2/2)*(sin(6*(360./NSide2)*PI/180.) + sin(7*(360./NSide2)*PI/180));      // considering
    // intersection of 2 medianes
    // Rotation:
    thx = 90.000000;
    phx = 0.000000 - (7 * 360. / NSide2);
    thy = 90. + InclAng2;
    phy = 90.000000 - (7 * 360. / NSide2);
    thz = InclAng2;
    phz = 90.000000 - (7 * 360. / NSide2);
    TGeoRotation* pMatrix81 = new TGeoRotation("", thx, phx, thy, phy, thz, phz);
    TGeoCombiTrans* pMatrix80 = new TGeoCombiTrans("", dx, dy, dz, pMatrix81);

    // Combi transformation:
    dx = dx - (WidthHalf2 / 2) * (1 / (cos(AngTrap2) * cos(AngTrap2))) *
                  (cos(7 * (360. / NSide2) * PI / 180.) +
                   cos(8 * (360. / NSide2) * PI / 180.)); // considering  real barycentre position;
    // dx =  dx -(WidthHalf2/2)*(cos(7*(360./NSide2)*PI/180.) + cos(8*(360./NSide2)*3.14159/180)); // considering
    // intersection of 2 medianes
    dy = dy + (WidthHalf2 / 2) * (1 / (cos(AngTrap2) * cos(AngTrap2))) *
                  (sin(7 * (360. / NSide2) * PI / 180.) +
                   sin(8 * (360. / NSide2) * PI / 180.)); // considering  real barycentre position;
    // dy =  dy +(WidthHalf2/2)*(sin(7*(360./NSide2)*PI/180.) + sin(8*(360./NSide2)*PI/180));      // considering
    // intersection of 2 medianes
    // Rotation:
    thx = 90.000000;
    phx = 0.000000 - (8 * 360. / NSide2);
    thy = 90. + InclAng2;
    phy = 90.000000 - (8 * 360. / NSide2);
    thz = InclAng2;
    phz = 90.000000 - (8 * 360. / NSide2);
    TGeoRotation* pMatrix83 = new TGeoRotation("", thx, phx, thy, phy, thz, phz);
    TGeoCombiTrans* pMatrix82 = new TGeoCombiTrans("", dx, dy, dz, pMatrix83);

    // Combi transformation:
    dx = dx - (WidthHalf2 / 2) * (1 / (cos(AngTrap2) * cos(AngTrap2))) *
                  (cos(8 * (360. / NSide2) * PI / 180.) +
                   cos(9 * (360. / NSide2) * PI / 180.)); // considering  real barycentre position;
    // dx =  dx -(WidthHalf2/2)*(cos(8*(360./NSide2)*PI/180.) + cos(9*(360./NSide2)*3.14159/180));  // considering
    // intersection of 2 medianes
    dy = dy + (WidthHalf2 / 2) * (1 / (cos(AngTrap2) * cos(AngTrap2))) *
                  (sin(8 * (360. / NSide2) * PI / 180.) +
                   sin(9 * (360. / NSide2) * PI / 180.)); // considering  real barycentre position;
    // dy =  dy +(WidthHalf2/2)*(sin(8*(360./NSide2)*PI/180.) + sin(9*(360./NSide2)*PI/180));       // considering
    // intersection of 2 medianes
    // Rotation:
    thx = 90.000000;
    phx = 0.000000 - (9 * 360. / NSide2);
    thy = 90. + InclAng2;
    phy = 90.000000 - (9 * 360. / NSide2);
    thz = InclAng2;
    phz = 90.000000 - (9 * 360. / NSide2);
    TGeoRotation* pMatrix85 = new TGeoRotation("", thx, phx, thy, phy, thz, phz);
    TGeoCombiTrans* pMatrix84 = new TGeoCombiTrans("", dx, dy, dz, pMatrix85);

    // Combi transformation:
    dx = dx - (WidthHalf2 / 2) * (1 / (cos(AngTrap2) * cos(AngTrap2))) *
                  (cos(9 * (360. / NSide2) * PI / 180.) +
                   cos(10 * (360. / NSide2) * PI / 180.)); // considering  real barycentre position;
    // dx =  dx -(WidthHalf2/2)*(cos(9*(360./NSide2)*PI/180.) + cos(10*(360./NSide2)*3.14159/180)); // considering
    // intersection of 2 medianes
    dy = dy + (WidthHalf2 / 2) * (1 / (cos(AngTrap2) * cos(AngTrap2))) *
                  (sin(9 * (360. / NSide2) * PI / 180.) +
                   sin(10 * (360. / NSide2) * PI / 180.)); // considering  real barycentre position;
    // dy =  dy +(WidthHalf2/2)*(sin(9*(360./NSide2)*PI/180.) + sin(10*(360./NSide2)*PI/180));     // considering
    // intersection of 2 medianes
    // Rotation:
    thx = 90.000000;
    phx = 0.000000 - (10 * 360. / NSide2);
    thy = 90. + InclAng2;
    phy = 90.000000 - (10 * 360. / NSide2);
    thz = InclAng2;
    phz = 90.000000 - (10 * 360. / NSide2);
    TGeoRotation* pMatrix87 = new TGeoRotation("", thx, phx, thy, phy, thz, phz);
    TGeoCombiTrans* pMatrix86 = new TGeoCombiTrans("", dx, dy, dz, pMatrix87);

    // Combi transformation:
    dx = dx - (WidthHalf2 / 2) * (1 / (cos(AngTrap2) * cos(AngTrap2))) *
                  (cos(10 * (360. / NSide2) * PI / 180.) +
                   cos(11 * (360. / NSide2) * PI / 180.)); // considering  real barycentre position;
    // dx =  dx -(WidthHalf2/2)*(cos(10*(360./NSide2)*PI/180.) + cos(11*(360./NSide2)*3.14159/180)); // considering
    // intersection of 2 medianes
    dy = dy + (WidthHalf2 / 2) * (1 / (cos(AngTrap2) * cos(AngTrap2))) *
                  (sin(10 * (360. / NSide2) * PI / 180.) +
                   sin(11 * (360. / NSide2) * PI / 180.)); // considering  real barycentre position;
    // dy =  dy +(WidthHalf2/2)*(sin(10*(360./NSide2)*PI/180.) + sin(11*(360./NSide2)*PI/180));      // considering
    // intersection of 2 medianes
    // Rotation:
    thx = 90.000000;
    phx = 0.000000 - (11 * 360. / NSide2);
    thy = 90. + InclAng2;
    phy = 90.000000 - (11 * 360. / NSide2);
    thz = InclAng2;
    phz = 90.000000 - (11 * 360. / NSide2);
    TGeoRotation* pMatrix89 = new TGeoRotation("", thx, phx, thy, phy, thz, phz);
    TGeoCombiTrans* pMatrix88 = new TGeoCombiTrans("", dx, dy, dz, pMatrix89);

    // Third layer

    // Combi transformation:
    dx = 0.000;
    dy = -((Length3 / 2 + (WidthHalf3 / 2) * tan(AngTrap3)) * sin(InclAng3 * PI / 180.) +
           Rmin3); // considering  real barycentre position
    // dy = -((Length3/2)*sin(InclAng3*PI/180.)+ Rmin3);                               // considering  intersection of 2
    // medianes
    dz = Length3 * cos(InclAng3 * PI / 180.) / 2 -
         (Length3 * cos(InclAng3 * PI / 180.) - (Rmin3 / tan(AngRangeMin3 * PI / 180.)));
    // Rotation:
    thx = 90.000000;
    phx = 0.000000;
    thy = 90. + InclAng3;
    phy = 90.000000;
    thz = InclAng3;
    phz = 90.000000;
    TGeoRotation* pMatrix35 = new TGeoRotation("", thx, phx, thy, phy, thz, phz); // geant3 angles
    TGeoCombiTrans* pMatrix34 = new TGeoCombiTrans("", dx, dy, dz, pMatrix35);

    // Combi transformation:

    dx = dx - (WidthHalf3 / 2) * (1 / (cos(AngTrap3) * cos(AngTrap3))) *
                  (1 + cos((360. / NSide3) * PI / 180.)); // considering  real barycentre position;
    // dx =  -(WidthHalf3/2)*(1 + cos((360./NSide3)*PI/180));                                     // considering
    // intersection of 2 medianes
    dy = dy + (WidthHalf3 / 2) * (1 / (cos(AngTrap3) * cos(AngTrap3))) *
                  sin((360. / NSide3) * PI / 180.); // considering  real barycentre position;
    // dy = dy +(WidthHalf3/2)*sin((360./NSide3)*PI/180);                                         // considering
    // intersection of 2 medianes
    // Rotation:
    thx = 90.000000;
    phx = 0.000000 - 360. / NSide3;
    thy = 90. + InclAng3;
    phy = 90.000000 - 360. / NSide3;
    thz = InclAng3;
    phz = 90.000000 - 360. / NSide3;
    TGeoRotation* pMatrix37 = new TGeoRotation("", thx, phx, thy, phy, thz, phz);
    TGeoCombiTrans* pMatrix36 = new TGeoCombiTrans("", dx, dy, dz, pMatrix37);

    // Combi transformation:
    dx = dx - (WidthHalf3 / 2) * (1 / (cos(AngTrap3) * cos(AngTrap3))) *
                  (cos((360. / NSide3) * PI / 180.) +
                   cos(2 * (360. / NSide3) * PI / 180.)); // considering  real barycentre position;
    // dx =  dx -(WidthHalf3/2)*(cos((360./NSide3)*PI/180.) + cos(2*(360./NSide3)*3.14159/180));   // considering
    // intersection of 2 medianes
    dy = dy + (WidthHalf3 / 2) * (1 / (cos(AngTrap3) * cos(AngTrap3))) *
                  (sin((360. / NSide3) * PI / 180.) +
                   sin(2 * (360. / NSide3) * PI / 180.)); // considering  real barycentre position;
    // dy =  dy +(WidthHalf3/2)*(sin((360./NSide3)*PI/180.) + sin(2*(360./NSide3)*PI/180));        // considering
    // intersection of 2 medianes
    // Rotation:
    thx = 90.000000;
    phx = 0.000000 - (2 * 360. / NSide3);
    thy = 90. + InclAng3;
    phy = 90.000000 - (2 * 360. / NSide3);
    thz = InclAng3;
    phz = 90.000000 - (2 * 360. / NSide3);
    TGeoRotation* pMatrix39 = new TGeoRotation("", thx, phx, thy, phy, thz, phz);
    TGeoCombiTrans* pMatrix38 = new TGeoCombiTrans("", dx, dy, dz, pMatrix39);

    // Combi transformation:
    dx = dx - (WidthHalf3 / 2) * (1 / (cos(AngTrap3) * cos(AngTrap3))) *
                  (cos(2 * (360. / NSide3) * PI / 180.) +
                   cos(3 * (360. / NSide3) * PI / 180.)); // considering  real barycentre position;
    // dx =  dx -(WidthHalf3/2)*(cos(2*(360./NSide3)*PI/180.) + cos(3*(360./NSide3)*3.14159/180)); // considering
    // intersection of 2 medianes
    dy = dy + (WidthHalf3 / 2) * (1 / (cos(AngTrap3) * cos(AngTrap3))) *
                  (sin(2 * (360. / NSide3) * PI / 180.) +
                   sin(3 * (360. / NSide3) * PI / 180.)); // considering  real barycentre position;
    // dy =  dy +(WidthHalf3/2)*(sin(2*(360./NSide3)*PI/180.) + sin(3*(360./NSide3)*PI/180));      // considering
    // intersection of 2 medianes
    // Rotation:
    thx = 90.000000;
    phx = 0.000000 - (3 * 360. / NSide3);
    thy = 90. + InclAng3;
    phy = 90.000000 - (3 * 360. / NSide3);
    thz = InclAng3;
    phz = 90.000000 - (3 * 360. / NSide3);
    TGeoRotation* pMatrix41 = new TGeoRotation("", thx, phx, thy, phy, thz, phz);
    TGeoCombiTrans* pMatrix40 = new TGeoCombiTrans("", dx, dy, dz, pMatrix41);

    // Combi transformation:
    dx = dx - (WidthHalf3 / 2) * (1 / (cos(AngTrap3) * cos(AngTrap3))) *
                  (cos(3 * (360. / NSide3) * PI / 180.) +
                   cos(4 * (360. / NSide3) * PI / 180.)); // considering  real barycentre position;
    // dx =  dx -(WidthHalf3/2)*(cos(3*(360./NSide3)*PI/180.) + cos(4*(360./NSide3)*3.14159/180)); // considering
    // intersection of 2 medianes
    dy = dy + (WidthHalf3 / 2) * (1 / (cos(AngTrap3) * cos(AngTrap3))) *
                  (sin(3 * (360. / NSide3) * PI / 180.) +
                   sin(4 * (360. / NSide3) * PI / 180.)); // considering  real barycentre position;
    // dy =  dy +(WidthHalf3/2)*(sin(3*(360./NSide3)*PI/180.) + sin(4*(360./NSide3)*PI/180));      // considering
    // intersection of 2 medianes
    // Rotation:
    thx = 90.000000;
    phx = 0.000000 - (4 * 360. / NSide3);
    thy = 90. + InclAng3;
    phy = 90.000000 - (4 * 360. / NSide3);
    thz = InclAng3;
    phz = 90.000000 - (4 * 360. / NSide3);
    TGeoRotation* pMatrix43 = new TGeoRotation("", thx, phx, thy, phy, thz, phz);
    TGeoCombiTrans* pMatrix42 = new TGeoCombiTrans("", dx, dy, dz, pMatrix43);

    // Combi transformation:
    dx = dx - (WidthHalf3 / 2) * (1 / (cos(AngTrap3) * cos(AngTrap3))) *
                  (cos(4 * (360. / NSide3) * PI / 180.) +
                   cos(5 * (360. / NSide3) * PI / 180.)); // considering  real barycentre position;
    // dx =  dx -(WidthHalf3/2)*(cos(4*(360./NSide3)*PI/180.) + cos(5*(360./NSide3)*3.14159/180)); // considering
    // intersection of 2 medianes
    dy = dy + (WidthHalf3 / 2) * (1 / (cos(AngTrap3) * cos(AngTrap3))) *
                  (sin(4 * (360. / NSide3) * PI / 180.) +
                   sin(5 * (360. / NSide3) * PI / 180.)); // considering  real barycentre position;
    // dy =  dy +(WidthHalf3/2)*(sin(4*(360./NSide3)*PI/180.) + sin(5*(360./NSide3)*PI/180));      // considering
    // intersection of 2 medianes
    // Rotation:
    thx = 90.000000;
    phx = 0.000000 - (5 * 360. / NSide3);
    thy = 90. + InclAng3;
    phy = 90.000000 - (5 * 360. / NSide3);
    thz = InclAng3;
    phz = 90.000000 - (5 * 360. / NSide3);
    TGeoRotation* pMatrix45 = new TGeoRotation("", thx, phx, thy, phy, thz, phz);
    TGeoCombiTrans* pMatrix44 = new TGeoCombiTrans("", dx, dy, dz, pMatrix45);

    // Combi transformation:
    dx = dx - (WidthHalf3 / 2) * (1 / (cos(AngTrap3) * cos(AngTrap3))) *
                  (cos(5 * (360. / NSide3) * PI / 180.) +
                   cos(6 * (360. / NSide3) * PI / 180.)); // considering  real barycentre position;
    // dx =  dx -(WidthHalf3/2)*(cos(5*(360./NSide3)*PI/180.) + cos(6*(360./NSide3)*3.14159/180)); // considering
    // intersection of 2 medianes
    dy = dy + (WidthHalf3 / 2) * (1 / (cos(AngTrap3) * cos(AngTrap3))) *
                  (sin(5 * (360. / NSide3) * PI / 180.) +
                   sin(6 * (360. / NSide3) * PI / 180.)); // considering  real barycentre position;
    // dy =  dy +(WidthHalf3/2)*(sin(5*(360./NSide3)*PI/180.) + sin(6*(360./NSide3)*PI/180));      // considering
    // intersection of 2 medianes
    // Rotation:
    thx = 90.000000;
    phx = 0.000000 - (6 * 360. / NSide3);
    thy = 90. + InclAng3;
    phy = 90.000000 - (6 * 360. / NSide3);
    thz = InclAng3;
    phz = 90.000000 - (6 * 360. / NSide3);
    TGeoRotation* pMatrix47 = new TGeoRotation("", thx, phx, thy, phy, thz, phz);
    TGeoCombiTrans* pMatrix46 = new TGeoCombiTrans("", dx, dy, dz, pMatrix47);

    // Combi transformation:
    dx = dx - (WidthHalf3 / 2) * (1 / (cos(AngTrap3) * cos(AngTrap3))) *
                  (cos(6 * (360. / NSide3) * PI / 180.) +
                   cos(7 * (360. / NSide3) * PI / 180.)); // considering  real barycentre position;
    // dx =  dx -(WidthHalf3/2)*(cos(6*(360./NSide3)*PI/180.) + cos(7*(360./NSide3)*3.14159/180)); // considering
    // intersection of 2 medianes
    dy = dy + (WidthHalf3 / 2) * (1 / (cos(AngTrap3) * cos(AngTrap3))) *
                  (sin(6 * (360. / NSide3) * PI / 180.) +
                   sin(7 * (360. / NSide3) * PI / 180.)); // considering  real barycentre position;
    // dy =  dy +(WidthHalf3/2)*(sin(6*(360./NSide3)*PI/180.) + sin(7*(360./NSide3)*PI/180));      // considering
    // intersection of 2 medianes
    // Rotation:
    thx = 90.000000;
    phx = 0.000000 - (7 * 360. / NSide3);
    thy = 90. + InclAng3;
    phy = 90.000000 - (7 * 360. / NSide3);
    thz = InclAng3;
    phz = 90.000000 - (7 * 360. / NSide3);
    TGeoRotation* pMatrix49 = new TGeoRotation("", thx, phx, thy, phy, thz, phz);
    TGeoCombiTrans* pMatrix48 = new TGeoCombiTrans("", dx, dy, dz, pMatrix49);

    // Combi transformation:
    dx = dx - (WidthHalf3 / 2) * (1 / (cos(AngTrap3) * cos(AngTrap3))) *
                  (cos(7 * (360. / NSide3) * PI / 180.) +
                   cos(8 * (360. / NSide3) * PI / 180.)); // considering  real barycentre position;
    // dx =  dx -(WidthHalf3/2)*(cos(7*(360./NSide3)*PI/180.) + cos(8*(360./NSide3)*3.14159/180)); // considering
    // intersection of 2 medianes
    dy = dy + (WidthHalf3 / 2) * (1 / (cos(AngTrap3) * cos(AngTrap3))) *
                  (sin(7 * (360. / NSide3) * PI / 180.) +
                   sin(8 * (360. / NSide3) * PI / 180.)); // considering  real barycentre position;
    // dy =  dy +(WidthHalf3/2)*(sin(7*(360./NSide3)*PI/180.) + sin(8*(360./NSide3)*PI/180));      // considering
    // intersection of 2 medianes
    // Rotation:
    thx = 90.000000;
    phx = 0.000000 - (8 * 360. / NSide3);
    thy = 90. + InclAng3;
    phy = 90.000000 - (8 * 360. / NSide3);
    thz = InclAng3;
    phz = 90.000000 - (8 * 360. / NSide3);
    TGeoRotation* pMatrix51 = new TGeoRotation("", thx, phx, thy, phy, thz, phz);
    TGeoCombiTrans* pMatrix50 = new TGeoCombiTrans("", dx, dy, dz, pMatrix51);

    // Combi transformation:
    dx = dx - (WidthHalf3 / 2) * (1 / (cos(AngTrap3) * cos(AngTrap3))) *
                  (cos(8 * (360. / NSide3) * PI / 180.) +
                   cos(9 * (360. / NSide3) * PI / 180.)); // considering  real barycentre position;
    // dx =  dx -(WidthHalf3/2)*(cos(8*(360./NSide3)*PI/180.) + cos(9*(360./NSide3)*3.14159/180)); // considering
    // intersection of 2 medianes
    dy = dy + (WidthHalf3 / 2) * (1 / (cos(AngTrap3) * cos(AngTrap3))) *
                  (sin(8 * (360. / NSide3) * PI / 180.) +
                   sin(9 * (360. / NSide3) * PI / 180.)); // considering  real barycentre position;
    // dy =  dy +(WidthHalf3/2)*(sin(8*(360./NSide3)*PI/180.) + sin(9*(360./NSide3)*PI/180));      // considering
    // intersection of 2 medianes
    // Rotation:
    thx = 90.000000;
    phx = 0.000000 - (9 * 360. / NSide3);
    thy = 90. + InclAng3;
    phy = 90.000000 - (9 * 360. / NSide3);
    thz = InclAng3;
    phz = 90.000000 - (9 * 360. / NSide3);
    TGeoRotation* pMatrix53 = new TGeoRotation("", thx, phx, thy, phy, thz, phz);
    TGeoCombiTrans* pMatrix52 = new TGeoCombiTrans("", dx, dy, dz, pMatrix53);

    // Combi transformation:
    dx = dx - (WidthHalf3 / 2) * (1 / (cos(AngTrap3) * cos(AngTrap3))) *
                  (cos(9 * (360. / NSide3) * PI / 180.) +
                   cos(10 * (360. / NSide3) * PI / 180.)); // considering  real barycentre position;
    // dx =  dx -(WidthHalf3/2)*(cos(9*(360./NSide3)*PI/180.) + cos(10*(360./NSide3)*3.14159/180)); // considering
    // intersection of 2 medianes
    dy = dy + (WidthHalf3 / 2) * (1 / (cos(AngTrap3) * cos(AngTrap3))) *
                  (sin(9 * (360. / NSide3) * PI / 180.) +
                   sin(10 * (360. / NSide3) * PI / 180.)); // considering  real barycentre position;
    // dy =  dy +(WidthHalf3/2)*(sin(9*(360./NSide3)*PI/180.) + sin(10*(360./NSide3)*PI/180));      // considering
    // intersection of 2 medianes
    // Rotation:
    thx = 90.000000;
    phx = 0.000000 - (10 * 360. / NSide3);
    thy = 90. + InclAng3;
    phy = 90.000000 - (10 * 360. / NSide3);
    thz = InclAng3;
    phz = 90.000000 - (10 * 360. / NSide3);
    TGeoRotation* pMatrix55 = new TGeoRotation("", thx, phx, thy, phy, thz, phz);
    TGeoCombiTrans* pMatrix54 = new TGeoCombiTrans("", dx, dy, dz, pMatrix55);

    // Combi transformation:
    dx = dx - (WidthHalf3 / 2) * (1 / (cos(AngTrap3) * cos(AngTrap3))) *
                  (cos(10 * (360. / NSide3) * PI / 180.) +
                   cos(11 * (360. / NSide3) * PI / 180.)); // considering  real barycentre position;
    // dx =  dx -(WidthHalf3/2)*(cos(10*(360./NSide3)*PI/180.) + cos(11*(360./NSide3)*3.14159/180)); // considering
    // intersection of 2 medianes
    dy = dy + (WidthHalf3 / 2) * (1 / (cos(AngTrap3) * cos(AngTrap3))) *
                  (sin(10 * (360. / NSide3) * PI / 180.) +
                   sin(11 * (360. / NSide3) * PI / 180.)); // considering  real barycentre position;
    // dy =  dy +(WidthHalf3/2)*(sin(10*(360./NSide3)*PI/180.) + sin(11*(360./NSide3)*PI/180));      // considering
    // intersection of 2 medianes
    // Rotation:
    thx = 90.000000;
    phx = 0.000000 - (11 * 360. / NSide3);
    thy = 90. + InclAng3;
    phy = 90.000000 - (11 * 360. / NSide3);
    thz = InclAng3;
    phz = 90.000000 - (11 * 360. / NSide3);
    TGeoRotation* pMatrix57 = new TGeoRotation("", thx, phx, thy, phy, thz, phz);
    TGeoCombiTrans* pMatrix56 = new TGeoCombiTrans("", dx, dy, dz, pMatrix57);

    // Chamber

    // Combi transformation:
    dx = 0.000000;
    dy = 0.000000;
    dz = -15.00000;
    // Rotation:
    thx = 90.000000;
    phx = 0.000000;
    thy = 90.000000;
    phy = 90.000000;
    thz = 0.000000;
    phz = 0.000000;
    TGeoRotation* pMatrix171 = new TGeoRotation("", thx, phx, thy, phy, thz, phz);
    TGeoCombiTrans* pMatrix170 = new TGeoCombiTrans("", dx, dy, dz, pMatrix171);

    // Combi transformation:
    dx = 0.000000;
    dy = 0.000000;
    dz = 20.000000;
    // Rotation:
    thx = 90.000000;
    phx = 0.000000;
    thy = 90.000000;
    phy = 90.000000;
    thz = 0.000000;
    phz = 0.000000;
    TGeoRotation* pMatrix173 = new TGeoRotation("", thx, phx, thy, phy, thz, phz);
    TGeoCombiTrans* pMatrix172 = new TGeoCombiTrans("", dx, dy, dz, pMatrix173);

    // Combi transformation:
    dx = 0.000000;
    dy = 0.000000;
    dz = -30.000000;
    // Rotation:
    thx = 90.000000;
    phx = 0.000000;
    thy = 90.000000;
    phy = 90.000000;
    thz = 0.000000;
    phz = 0.000000;
    TGeoRotation* pMatrix175 = new TGeoRotation("", thx, phx, thy, phy, thz, phz);
    TGeoCombiTrans* pMatrix174 = new TGeoCombiTrans("", dx, dy, dz, pMatrix175);

    // Gold Cone
    // Combi transformation:
    dx = 0.000000;
    dy = 0.000000;
    dz = 7.10000;
    // Rotation:
    thx = 90.000000;
    phx = 0.000000;
    thy = 90.000000;
    phy = 90.000000;
    thz = 0.000000;
    phz = 0.000000;
    TGeoRotation* pMatrix177 = new TGeoRotation("", thx, phx, thy, phy, thz, phz);
    TGeoCombiTrans* pMatrix176 = new TGeoCombiTrans("", dx, dy, dz, pMatrix177);

    // Shape: World type: TGeoBBox
    TGeoVolume* pWorld = gGeoManager->GetTopVolume();
    pWorld->SetVisLeaves(kTRUE);

    /*

    // Create a global Mother Volume
    dx = 100.000000;
    dy = 100.000000;
    dz = 100.000000;
    TGeoShape *pBoxWorld = new TGeoBBox("TraBoxWorld", dx,dy,dz);
    TGeoVolume*
    pWorld  = new TGeoVolume("TraBoxLogWorld",pBoxWorld, pMed2);
    pWorld->SetVisLeaves(kTRUE);
    TGeoCombiTrans *pGlobalc = GetGlobalPosition();

    // add the sphere as Mother Volume
    pAWorld->AddNode(pWorld, 0, pGlobalc);

    */

    // SHAPES, VOLUMES AND GEOMETRICAL HIERARCHY
    // Shape: SiVacuumSphereWorld type: TGeoSphere

    // Si Shape & volume: TraBox type: TGeoBBox
    // dx = 15.00000;
    // dy = 2.000000;
    // dz = 0.005000;
    // Volume: TraLog
    // TGeoVolume *TraLog1 = gGeoManager->MakeBox("TraLog1",pMed22,dx,dy,dz);
    TGeoVolume* TraLog1 =
        gGeoManager->MakeTrd1("TraLog1", pMed22, WidthMax1 / 2., WidthMin1 / 2., Thickness1, Length1 / 2.);

    // Si Shape & volume: TraBox type: TGeoBBox
    dx = 15.00000;
    dy = 4.50000;
    dz = 0.015000;
    // Volume: TraLog
    // TGeoVolume *TraLog2 = gGeoManager->MakeBox("TraLog2",pMed22,dx,dy,dz);
    TGeoVolume* TraLog2 =
        gGeoManager->MakeTrd1("TraLog2", pMed22, WidthMax3 / 2, WidthMin3 / 2, Thickness3, Length3 / 2);

    // Si Shape & volume: TraBox type: TGeoBBox
    dx = 15.00000;
    dy = 5.500000;
    dz = 0.015000;
    // Volume: TraLog
    // TGeoVolume *TraLog3 = gGeoManager->MakeBox("TraLog3",pMed22,dx,dy,dz);
    TGeoVolume* TraLog3 =
        gGeoManager->MakeTrd1("TraLog3", pMed22, WidthMax2 / 2, WidthMin2 / 2, Thickness2, Length2 / 2);

    // Si Shape & volume: TraBox type: TGeoBBox
    dx = 13.500000;
    dy = 7.1500000;
    dz = 0.015000;
    // Volume: TraLog
    TGeoVolume* TraLog4 = gGeoManager->MakeBox("TraLog4", pMed22, dx, dy, dz);

    // Si Shape & volume: TraBox type: TGeoBBox
    dx = 13.500000;
    dy = 8.3500000;
    dz = 0.015000;
    // Volume: TraLog
    TGeoVolume* TraLog5 = gGeoManager->MakeBox("TraLog5", pMed22, dx, dy, dz);

    // Shape: MontagePlatform type: TGeoTubeSeg
    rmin = 5.0000;
    rmax = 11.000000;
    dz = 0.0150000;
    phi1 = 0.000000;
    phi2 = 360.000000;
    TGeoShape* Tra6 = new TGeoTubeSeg("Tra6", rmin, rmax, dz, phi1, phi2);
    // Volume: MontagePlatformLog
    TGeoVolume* TraLog6 = new TGeoVolume("TraLog6", Tra6, pMed22);

    // Shape: Chamber type: TGeoTubeSeg
    rmin = 25.00000;
    rmax = 25.200000;
    dz = 15.0000;
    phi1 = 0.000000;
    phi2 = 360.000000;
    TGeoShape* pChamberBarrel = new TGeoTubeSeg("ChamberBarrel", rmin, rmax, dz, phi1, phi2);
    // Volume: ChamberLog
    TGeoVolume* pChamberBarrelLog = new TGeoVolume("ChamberBarrelLog", pChamberBarrel, pMed3);

    // Shape: Chamber type: TGeoCone
    Double_t rmax1 = 25.20000;
    Double_t rmin1 = 25.40000;
    Double_t rmax2 = 2.650000;
    Double_t rmin2 = 2.450000;
    dz = 20.0000;
    TGeoShape* pChamberEnd1 = new TGeoCone("ChamberEnd1", dz, rmin1, rmax1, rmin2, rmax2);
    // Volume: pGoldCone
    TGeoVolume* pChamberEnd1Log = new TGeoVolume("ChamberEnd1Log", pChamberEnd1, pMed3);

    // Shape: Chamber type: TGeoTubeSeg
    rmin = 2.00000;
    rmax = 25.200000;
    dz = 0.20000;
    phi1 = 0.000000;
    phi2 = 360.000000;
    TGeoShape* pChamberEnd2 = new TGeoTubeSeg("ChamberEnd2", rmin, rmax, dz, phi1, phi2);
    // Volume: ChamberLog
    TGeoVolume* pChamberEnd2Log = new TGeoVolume("ChamberEnd2Log", pChamberEnd2, pMed3);

    // Shape: MontagePlatform type: TGeoTubeSeg
    rmin = 2.750000;
    rmax = 18.000000;
    dz = 0.150000;
    phi1 = 0.000000;
    phi2 = 360.000000;
    TGeoShape* pMontagePlatform = new TGeoTubeSeg("MontagePlatform", rmin, rmax, dz, phi1, phi2);
    // Volume: MontagePlatformLog
    TGeoVolume* pMontagePlatformLog = new TGeoVolume("MontagePlatformLog", pMontagePlatform, pMed25);

    // Shape: BeamPipe type: TGeoTubeSeg
    rmin = 1.800000;
    rmax = 2.000000;
    dz = 20.00000;
    phi1 = 0.000000;
    phi2 = 360.000000;
    TGeoShape* pBeamPipe = new TGeoTubeSeg("BeamPipe", rmin, rmax, dz, phi1, phi2);
    // Volume: BeamPipeLog
    TGeoVolume* pBeamPipeLog = new TGeoVolume("BeamPipeLog", pBeamPipe, pMed21);

    // Shape: MontageRing type: TGeoTubeSeg
    rmin = 12.000000;
    rmax = 15.000000;
    dz = 0.250000;
    phi1 = 0.000000;
    phi2 = 360.000000;
    TGeoShape* pMontageRing = new TGeoTubeSeg("MontageRing", rmin, rmax, dz, phi1, phi2);
    // Volume: MontageRingLog
    TGeoVolume* pMontageRingLog = new TGeoVolume("MontageRingLog", pMontageRing, pMed25);

    // Shape: targetWheel type: TGeoTubeSeg
    rmin = 0.410000;
    rmax = 5.300000;
    dz = 0.100000;
    phi1 = 0.000000;
    phi2 = 360.000000;
    TGeoShape* ptargetWheel = new TGeoTubeSeg("targetWheel", rmin, rmax, dz, phi1, phi2);
    // Volume: targetWheelLog
    TGeoVolume* ptargetWheelLog = new TGeoVolume("targetWheelLog", ptargetWheel, pMed21);

    // Shape: targetWheel2 type: TGeoTubeSeg
    rmin = 9.500000;
    rmax = 10.000000;
    dz = 0.100000;
    phi1 = 0.000000;
    phi2 = 360.000000;
    TGeoShape* ptargetWheel2 = new TGeoTubeSeg("targetWheel2", rmin, rmax, dz, phi1, phi2);
    // Volume: targetWheel2Log
    TGeoVolume* ptargetWheel2Log = new TGeoVolume("targetWheel2Log", ptargetWheel2, pMed21);

    // Shape: Target6 type: TGeoCone
    // Double_t rmax1 = 6.270000;
    // Double_t rmin1 = 6.260000;
    // Double_t rmax2 = 1.930000;
    // Double_t rmin2 = 1.920000;
    // dz = 7.0000;
    // TGeoShape *pGoldCone = new TGeoCone("GoldCone",dz,rmin1,rmax1,rmin2,rmax2);
    // Volume: pGoldCone
    // TGeoVolume*
    // pGoldConelog = new TGeoVolume("GoldConelog",pGoldCone, pMed79);

    // Shape: innerElectronicBox type: TGeoBBox
    dx = 3.800000;
    dy = 3.300000;
    dz = 0.800000;
    TGeoShape* pinnerElectronicBox = new TGeoBBox("innerElectronicBox", dx, dy, dz);
    // Volume: innerElectronicsLog
    TGeoVolume* pinnerElectronicsLog = new TGeoVolume("innerElectronicsLog", pinnerElectronicBox, pMed2);

    //
    // Make elementary assembly of the whole structure.
    //

    TGeoVolume* aTra = new TGeoVolumeAssembly("ATRA");

    // aTra->AddNode(ptargetWheel2Log,1, pMatrix30);

    // aTra->AddNode(ptargetWheelLog ,1, pMatrix28);

    // aTra->AddNode(pMontageRingLog ,1, pMatrix16);

    // aTra->AddNode(pBeamPipeLog,1,pMatrix58);

    aTra->AddNode(pChamberBarrelLog, 1, pMatrix170);
    aTra->AddNode(pChamberEnd1Log, 1, pMatrix172);
    aTra->AddNode(pChamberEnd2Log, 1, pMatrix174);
    // aTra->AddNode(pGoldConelog, 1, pMatrix176);

    AddSensitiveVolume(TraLog1);
    AddSensitiveVolume(TraLog2);
    AddSensitiveVolume(TraLog3);

    fNbOfSensitiveVol += 1;

    // First layer

    aTra->AddNode(TraLog1, 1, pMatrix2);
    aTra->AddNode(TraLog1, 2, pMatrix4);
    aTra->AddNode(TraLog1, 3, pMatrix6);
    aTra->AddNode(TraLog1, 4, pMatrix8);
    aTra->AddNode(TraLog1, 5, pMatrix10);
    aTra->AddNode(TraLog1, 6, pMatrix12);

    // Second layer

    aTra->AddNode(TraLog2, 7, pMatrix66);
    aTra->AddNode(TraLog2, 8, pMatrix68);
    aTra->AddNode(TraLog2, 9, pMatrix70);
    aTra->AddNode(TraLog2, 10, pMatrix72);
    aTra->AddNode(TraLog2, 11, pMatrix74);
    aTra->AddNode(TraLog2, 12, pMatrix76);
    aTra->AddNode(TraLog2, 13, pMatrix78);
    aTra->AddNode(TraLog2, 14, pMatrix80);
    aTra->AddNode(TraLog2, 15, pMatrix82);
    aTra->AddNode(TraLog2, 16, pMatrix84);
    aTra->AddNode(TraLog2, 17, pMatrix86);
    aTra->AddNode(TraLog2, 18, pMatrix88);

    // Third layer

    aTra->AddNode(TraLog3, 19, pMatrix34);
    aTra->AddNode(TraLog3, 20, pMatrix36);
    aTra->AddNode(TraLog3, 21, pMatrix38);
    aTra->AddNode(TraLog3, 22, pMatrix40);
    aTra->AddNode(TraLog3, 23, pMatrix42);
    aTra->AddNode(TraLog3, 24, pMatrix44);
    aTra->AddNode(TraLog3, 25, pMatrix46);
    aTra->AddNode(TraLog3, 26, pMatrix48);
    aTra->AddNode(TraLog3, 27, pMatrix50);
    aTra->AddNode(TraLog3, 28, pMatrix52);
    aTra->AddNode(TraLog3, 29, pMatrix54);
    aTra->AddNode(TraLog3, 30, pMatrix56);

    TGeoRotation* rotg = new TGeoRotation();
    rotg->RotateX(0.);
    rotg->RotateY(0.);
    rotg->RotateZ(0.);
    dx = tx = 0.0;
    dy = ty = 0.0;
    dz = tz = 0.0;

    TGeoCombiTrans* t0 = new TGeoCombiTrans(tx, ty, tz, rotg);

    // aTra->AddNode(TraLog1,2, GetGlobalPosition(t0));
    // aTra->AddNode(TraLog3,2, GetGlobalPosition(t0));
    pWorld->AddNode(aTra, 1, GetGlobalPosition(t0));
}

/*
void R3BTra::ConstructGeometry() {

  FairGeoLoader*    geoLoad = FairGeoLoader::Instance();
  FairGeoInterface* geoFace = geoLoad->getGeoInterface();
  R3BGeoTra*       stsGeo  = new R3BGeoTra();
  stsGeo->setGeomFile(GetGeometryFileName());
  geoFace->addGeoModule(stsGeo);

  Bool_t rc = geoFace->readSet(stsGeo);
  if (rc) stsGeo->create(geoLoad->getGeoBuilder());
  TList* volList = stsGeo->getListOfVolumes();
  // store geo parameter
  FairRun *fRun = FairRun::Instance();
  FairRuntimeDb *rtdb= FairRun::Instance()->GetRuntimeDb();
  R3BGeoTraPar* par=(R3BGeoTraPar*)(rtdb->getContainer("R3BGeoTraPar"));
  TObjArray *fSensNodes = par->GetGeoSensitiveNodes();
  TObjArray *fPassNodes = par->GetGeoPassiveNodes();

  TListIter iter(volList);
  FairGeoNode* node   = NULL;
  FairGeoVolume *aVol=NULL;

  while( (node = (FairGeoNode*)iter.Next()) ) {
      aVol = dynamic_cast<FairGeoVolume*> ( node );
       if ( node->isSensitive()  ) {
           fSensNodes->AddLast( aVol );
       }else{
           fPassNodes->AddLast( aVol );
       }
  }
  par->setChanged();
  par->setInputVersion(fRun->GetRunId(),1);
  ProcessNodes( volList );

}
*/

ClassImp(R3BTra)
