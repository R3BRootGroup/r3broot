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

// ------------------------------------------------------------
// -----         R3BMwpcDigitizer source file             -----
// -----      Created 03/11/19 by JL Rodriguez            -----
// ------------------------------------------------------------

#include "R3BMwpcDigitizer.h"
#include "R3BLogger.h"
#include "R3BMCTrack.h"
#include "R3BMwpcPoint.h"
#include "R3BTGeoPar.h"

#include "FairLogger.h"
#include "FairRootManager.h"
#include "FairRuntimeDb.h"

#include "TClonesArray.h"
#include "TMath.h"
#include "TRandom.h"
#include "TVector3.h"
#include <iostream>
#include <string>

// R3BMwpcDigitizer: Default Constructor --------------------------
R3BMwpcDigitizer::R3BMwpcDigitizer()
    : R3BMwpcDigitizer("R3BMwpcDigitizer", 1)
{
}

// R3BMwpcDigitizer: Standard Constructor --------------------------
R3BMwpcDigitizer::R3BMwpcDigitizer(const TString& name, Int_t iVerbose)
    : FairTask(name + "Digi", iVerbose)
    , fName(name)
    , fMCTrack(NULL)
    , fMwpcPoints(NULL)
    , fMwpcHits(NULL)
    , fsigma_x(0.125)
    , fsigma_y(0.25)
{
}

// Virtual R3BMwpcDigitizer: Destructor ----------------------------
R3BMwpcDigitizer::~R3BMwpcDigitizer()
{
    R3BLOG(debug, " for " << fName);
    if (fMwpcHits)
    {
        delete fMwpcHits;
    }
}

void R3BMwpcDigitizer::SetParContainers()
{
    FairRuntimeDb* rtdb = FairRuntimeDb::instance();

    fMwpcGeoPar = dynamic_cast<R3BTGeoPar*>(rtdb->getContainer(fName + "GeoPar"));
    if (!fMwpcGeoPar)
    {
        R3BLOG(error, "Could not get access to " + fName + "GeoPar container.");
        return;
    }
    else
    {
        R3BLOG(info, "Container " + fName + "GeoPar found.");
    }
}

void R3BMwpcDigitizer::SetParameter()
{
    fsigma_x = fMwpcGeoPar->GetSigmaX();
    fsigma_y = fMwpcGeoPar->GetSigmaY();

    fRot.RotateX(-fMwpcGeoPar->GetRotX() * TMath::DegToRad());
    fRot.RotateY(-fMwpcGeoPar->GetRotY() * TMath::DegToRad());
    fRot.RotateZ(-fMwpcGeoPar->GetRotZ() * TMath::DegToRad());

    fTrans.SetXYZ(fMwpcGeoPar->GetPosX(), fMwpcGeoPar->GetPosY(), fMwpcGeoPar->GetPosZ());
}

// ----   Public method Init  -----------------------------------------
InitStatus R3BMwpcDigitizer::Init()
{
    R3BLOG(info, " for " << fName);

    // Get input array
    FairRootManager* ioman = FairRootManager::Instance();
    R3BLOG_IF(fatal, !ioman, "FairRootManager not found");

    fMCTrack = dynamic_cast<TClonesArray*>(ioman->GetObject("MCTrack"));
    fMwpcPoints = dynamic_cast<TClonesArray*>(ioman->GetObject(fName + "Point"));
    R3BLOG_IF(fatal, !fMwpcPoints, fName + "Point not found");

    // Register output array fMwpcHits
    fMwpcHits = new TClonesArray("R3BMwpcHitData");
    ioman->Register(fName + "HitData", "Digital response in " + fName, fMwpcHits, kTRUE);

    SetParameter();
    return kSUCCESS;
}

// -----   Public method Execution   --------------------------------------------
void R3BMwpcDigitizer::Exec(Option_t* opt)
{
    Reset();
    // Reading the Input -- Point Data --
    Int_t nHits = fMwpcPoints->GetEntriesFast();
    if (nHits == 0)
        return;
    // Data from Point level
    R3BMwpcPoint** pointData;
    pointData = new R3BMwpcPoint*[nHits];
    Int_t TrackId = 0, PID = 0;
    Double_t x = 0., y = 0., z = 0.;
    TVector3 vpos;
    for (Int_t i = 0; i < nHits; i++)
    {
        pointData[i] = dynamic_cast<R3BMwpcPoint*>(fMwpcPoints->At(i));
        TrackId = pointData[i]->GetTrackID();

        R3BMCTrack* Track = dynamic_cast<R3BMCTrack*>(fMCTrack->At(TrackId));
        PID = Track->GetPdgCode();

        if (PID > 1000080160) // Z=8 and A=16
        {
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

            vpos = fRot * (vpos - fTrans);

            AddHitData(vpos.X() * 10. + gRandom->Gaus(0., fsigma_x),
                       vpos.Y() * 10. + gRandom->Gaus(0., fsigma_y)); // mm
        }
    }
    if (pointData)
    {
        delete[] pointData;
    }

    LOG(info) << "R3B" << fName << "Digitizer: " << fMwpcHits->GetEntriesFast() << " points registered in this event";
    return;
}

// -----   Public method ReInit   ----------------------------------------------
InitStatus R3BMwpcDigitizer::ReInit()
{
    SetParContainers();
    SetParameter();
    return kSUCCESS;
}

// -----   Public method Reset   -----------------------------------------------
void R3BMwpcDigitizer::Reset()
{
    R3BLOG(debug, " for " << fName);
    if (fMwpcHits)
    {
        fMwpcHits->Clear();
    }
}

// -----   Private method AddHitData  -------------------------------------------
R3BMwpcHitData* R3BMwpcDigitizer::AddHitData(Double_t x, Double_t y)
{
    // It fills the R3BMwpcHitData
    TClonesArray& clref = *fMwpcHits;
    Int_t size = clref.GetEntriesFast();
    return new (clref[size]) R3BMwpcHitData(x, y);
}

ClassImp(R3BMwpcDigitizer);
