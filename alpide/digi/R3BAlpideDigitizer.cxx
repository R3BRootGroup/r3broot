/******************************************************************************
 *   Copyright (C) 2022 GSI Helmholtzzentrum für Schwerionenforschung GmbH    *
 *   Copyright (C) 2022-2025 Members of R3B Collaboration                     *
 *                                                                            *
 *             This software is distributed under the terms of the            *
 *                 GNU General Public Licence (GPL) version 3,                *
 *                    copied verbatim in the file "LICENSE".                  *
 *                                                                            *
 * In applying this license GSI does not waive the privileges and immunities  *
 * granted to it by virtue of its status as an Intergovernmental Organization *
 * or submit itself to any jurisdiction.                                      *
 ******************************************************************************/

// ------------------------------------------------------------
// -----         R3BAlpideDigitizer source file           -----
// -----    Created 12/10/22 by J.L. Rodriguez Sanchez    -----
// ------------------------------------------------------------

#include "R3BAlpideDigitizer.h"
#include "R3BAlpideGeometry.h"
#include "R3BAlpideMappingPar.h"
#include "R3BAlpidePoint.h"
#include "R3BLogger.h"
#include "R3BMCTrack.h"
#include "R3BTGeoPar.h"

#include <FairRootManager.h>
#include <FairRuntimeDb.h>
#include <TClonesArray.h>

#include <TMath.h>
#include <TRandom.h>
#include <TVector3.h>
#include <string>

// R3BAlpideDigitizer: Default Constructor --------------------------
R3BAlpideDigitizer::R3BAlpideDigitizer()
    : R3BAlpideDigitizer("R3BAlpideDigitizer", 1)
{
}

// R3BAlpideDigitizer: Standard Constructor --------------------------
R3BAlpideDigitizer::R3BAlpideDigitizer(const TString& name, Int_t iVerbose)
    : FairTask(name + "Digitizer", iVerbose)
    , fName(name)
    , fGeoversion(2024)
    , fMCTrack(nullptr)
    , fAlpidePoints(nullptr)
    , fAlpideHits(nullptr)
    , fAlpideGeo(nullptr)
    , fMappingPar(nullptr)
    , fLabframe(false)
    , fsigma(0.0005) // in cm
{
}

// Virtual R3BAlpideDigitizer: Destructor ----------------------------
R3BAlpideDigitizer::~R3BAlpideDigitizer()
{
    R3BLOG(debug1, "");
    if (fAlpideHits)
    {
        fAlpideHits->Delete();
        delete fAlpideHits;
    }
}

void R3BAlpideDigitizer::SetParContainers()
{
    FairRuntimeDb* rtdb = FairRuntimeDb::instance();
    fMappingPar = dynamic_cast<R3BAlpideMappingPar*>(rtdb->getContainer("alpideMappingPar"));
    R3BLOG_IF(warn, !fMappingPar, "Could not get access to alpideMappingPar");
    R3BLOG_IF(info, fMappingPar, "Container alpideMappingPar found.");
}

void R3BAlpideDigitizer::SetParameter()
{
    //--- Parameter Container ---
    if (fMappingPar)
    {
        fGeoversion = fMappingPar->GetGeoVersion();
        R3BLOG(info, "Geometry version: " << fGeoversion);
    }
}

// ----   Public method Init  -----------------------------------------
InitStatus R3BAlpideDigitizer::Init()
{
    R3BLOG(info, "Init digitizer for " << fName);

    FairRootManager* ioman = FairRootManager::Instance();
    R3BLOG_IF(fatal, !ioman, "FairRootManager not found.");

    fMCTrack = static_cast<TClonesArray*>(ioman->GetObject("MCTrack"));
    fAlpidePoints = dynamic_cast<TClonesArray*>(ioman->GetObject(fName + "Point"));
    R3BLOG_IF(fatal, !fAlpidePoints, fName << "Point not found.");

    // Register output array fAlpideHits
    fAlpideHits = new TClonesArray("R3BAlpideHitData");
    ioman->Register(fName + "HitData", "Digital response in " + fName, fAlpideHits, kTRUE);

    SetParameter();

    fAlpideGeo = R3BAlpideGeometry::Instance();
    R3BLOG_IF(error, !fAlpideGeo->Init(fGeoversion), "Alpide geometry " << fGeoversion << " not found");

    return kSUCCESS;
}

// -----   Public method Execution   --------------------------------------------
void R3BAlpideDigitizer::Exec(Option_t*)
{
    Reset();
    // Reading the Input -- Point Data --
    Int_t nHits = fAlpidePoints->GetEntriesFast();
    if (nHits == 0)
    {
        return;
    }
    // Data from Point level
    R3BAlpidePoint** pointData = new R3BAlpidePoint*[nHits];
    Int_t TrackId = 0, PID = 0;
    Double_t x = 0., y = 0., z = 0.;
    TVector3 vpos;
    for (Int_t i = 0; i < nHits; i++)
    {
        fRot.SetToIdentity();
        pointData[i] = dynamic_cast<R3BAlpidePoint*>(fAlpidePoints->At(i));
        TrackId = pointData[i]->GetTrackID();
        auto sid = pointData[i]->GetSensorID();

        auto Track = dynamic_cast<R3BMCTrack*>(fMCTrack->At(TrackId));
        PID = Track->GetPdgCode();

        // if (PID > 1000080160) // Z=8 and A=16
        // {
        Double_t fX_in = pointData[i]->GetXIn();
        Double_t fY_in = pointData[i]->GetYIn();
        Double_t fZ_in = pointData[i]->GetZIn();
        Double_t fX_out = pointData[i]->GetXOut();
        Double_t fY_out = pointData[i]->GetYOut();
        Double_t fZ_out = pointData[i]->GetZOut();

        x = ((fX_in + fX_out) / 2.);
        y = ((fY_in + fY_out) / 2.);
        z = ((fZ_in + fZ_out) / 2.);
        vpos.SetXYZ(x, y, z);

        fTrans = fAlpideGeo->GetTranslation(sid);
        TVector3 vtpos = (vpos - fTrans);

        fRot = fAlpideGeo->GetRotation(sid);
        // std::cout <<"Rot "<< fRot.XX() <<" "<< fRot.XY() <<" "<< fRot.XZ() << std::endl;
        // std::cout <<"Rot "<< fRot.YX() <<" "<< fRot.YY() <<" "<< fRot.YZ() << std::endl;
        // std::cout <<"Rot "<< fRot.ZX() <<" "<< fRot.ZY() <<" "<< fRot.ZZ() << std::endl;

        TVector3 localpos = fRot.Inverse() * vtpos;
        localpos.SetX(gRandom->Gaus(0., fsigma) + localpos.X());
        localpos.SetZ(gRandom->Gaus(0., fsigma) + localpos.Z());

        auto eLoss = static_cast<double>(pointData[i]->GetEnergyLoss());

        auto fClusterSize = static_cast<uint16_t>(eLoss * pixelSizeFactor + pixelSizeOffset);

        if (fLabframe)
        {
            // Lab frame
            TVector3 labpos = fRot * localpos + fTrans;
            AddHitData(sid, fClusterSize, labpos.X() * 10., labpos.Y() * 10., labpos.Z() * 10.); // mm
        }
        else
        {
            // Sensor frame
            AddHitData(sid, fClusterSize, localpos.Z() * 10., localpos.X() * 10.); // mm
        }
        //}
    }
    if (pointData)
    {
        delete[] pointData;
    }
    R3BLOG(info, fAlpideHits->GetEntriesFast() << " points registered in this event");
    return;
}

// -----   Public method ReInit   ----------------------------------------------
InitStatus R3BAlpideDigitizer::ReInit()
{
    SetParContainers();
    SetParameter();
    return kSUCCESS;
}

// -----   Public method Reset   -----------------------------------------------
void R3BAlpideDigitizer::Reset()
{
    R3BLOG(debug, "Clearing R3B" + fName + "digitizer structure");
    if (fAlpideHits)
    {
        fAlpideHits->Clear();
    }
    fRot.SetToIdentity();
    return;
}

// -----   Private method AddHitData  -------------------------------------------
R3BAlpideHitData* R3BAlpideDigitizer::AddHitData(UInt_t sid, uint16_t clustersize, Double_t x, Double_t y, Double_t z)
{
    // It fills the R3BAlpideHitData
    TClonesArray& clref = *fAlpideHits;
    Int_t size = clref.GetEntriesFast();
    return new (clref[size]) R3BAlpideHitData(sid, clustersize, x, y, z);
}

ClassImp(R3BAlpideDigitizer)
