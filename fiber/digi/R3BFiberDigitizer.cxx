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
// -----                    R3BFiberDigitizer                -----
// -----          Created 10/05/2022 by J.L. Rodriguez       -----
// ---------------------------------------------------------------

#include "R3BFiberDigitizer.h"
#include "FairRootManager.h"
#include "FairRuntimeDb.h"
#include "R3BFibPoint.h"
#include "R3BLogger.h"
#include "R3BMCTrack.h"
#include "R3BTGeoPar.h"
#include "TClonesArray.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TMath.h"
#include "TRandom.h"
#include "TVector3.h"
#include <iostream>
#include <string>

R3BFiberDigitizer::R3BFiberDigitizer(const TString& name)
    : R3BFiberDigitizer(name, 0.001, 0.01, 1.0)
{
}

R3BFiberDigitizer::R3BFiberDigitizer(const TString& name, Double_t e, Double_t t, Double_t y)
    : FairTask("R3B" + name + "Digitizer")
    , fName(name)
    , fFiPoints(NULL)
    , fFiHits(NULL)
{
    fMinPID = 1000020030; // Z=2 and A=3
    esigma = e;
    tsigma = t;
    ysigma = y;
    fiber_thickness = 0.10000;
    fiber_nbr = 512;
    air_layer = 0.01; // relative to fiber_thickness
    detector_width = fiber_nbr * fiber_thickness * (1 + air_layer);
}

R3BFiberDigitizer::~R3BFiberDigitizer()
{
    if (fFiPoints)
        delete fFiPoints;
}

void R3BFiberDigitizer::SetEnergyResolution(Double_t e) { esigma = e; }

void R3BFiberDigitizer::SetTimeResolution(Double_t t) { tsigma = t; }

void R3BFiberDigitizer::SetYPositionResolution(Double_t y) { ysigma = y; }

void R3BFiberDigitizer::SetParContainers()
{
    FairRuntimeDb* rtdb = FairRuntimeDb::instance();
    fFiGeoPar = dynamic_cast<R3BTGeoPar*>(rtdb->getContainer(fName + "GeoPar"));
    if (!fFiGeoPar)
    {
        R3BLOG(error, "R3BFiberDigitizer::SetParContainers() : Could not get access to " + fName + "GeoPar container.");
        return;
    }
    else
    {
        R3BLOG(info, "R3BFiberDigitizer::SetParContainers() : Container " + fName + "GeoPar found.");
    }
}

void R3BFiberDigitizer::SetParameter()
{
    if (fFiGeoPar)
    {
        ysigma = fFiGeoPar->GetSigmaY();
        xsigma = fFiGeoPar->GetSigmaX();
        fRot.RotateX(-fFiGeoPar->GetRotX() * TMath::DegToRad());
        fRot.RotateY(-fFiGeoPar->GetRotY() * TMath::DegToRad());
        fRot.RotateZ(-fFiGeoPar->GetRotZ() * TMath::DegToRad());

        fTrans.SetXYZ(fFiGeoPar->GetPosX(), fFiGeoPar->GetPosY(), fFiGeoPar->GetPosZ());
    }
}

InitStatus R3BFiberDigitizer::Init()
{
    R3BLOG(info, "for fiber " << fName);
    // Get input array
    FairRootManager* ioman = FairRootManager::Instance();
    R3BLOG_IF(fatal, !ioman, "FairRootManager not found");

    fFiPoints = dynamic_cast<TClonesArray*>(ioman->GetObject(fName + "Point"));
    R3BLOG_IF(fatal, !fFiPoints, fName + "Point not found");

    fMCTrack = dynamic_cast<TClonesArray*>(ioman->GetObject("MCTrack"));

    // Register output array DchDigi
    fFiHits = new TClonesArray("R3BBunchedFiberHitData");
    ioman->Register(fName + "Hit", "Digital response in" + fName, fFiHits, kTRUE);

    // for sigmas
    rand = new TRandom3();

    SetParameter();
    return kSUCCESS;
}

void R3BFiberDigitizer::Exec(Option_t* opt)
{
    Reset();
    // Reading the Input -- Point Data --
    Int_t nHits = fFiPoints->GetEntriesFast();
    if (nHits == 0)
    {
        return;
    }
    // Data from Point level
    R3BFibPoint** pointData = new R3BFibPoint*[nHits];
    Int_t fiber = 0;
    Int_t TrackId = 0, PID = 0, mother = -1;
    Double_t x = 0., y = 0., z = 0., time = 0.;
    TVector3 vpos;

    for (Int_t i = 0; i < nHits; i++)
    {
        pointData[i] = dynamic_cast<R3BFibPoint*>(fFiPoints->At(i));
        TrackId = pointData[i]->GetTrackID();

        R3BMCTrack* Track = dynamic_cast<R3BMCTrack*>(fMCTrack->At(TrackId));
        PID = Track->GetPdgCode();
        // mother = Track->GetMotherId();

        if (PID >= fMinPID)
        {
            x = (pointData[i]->GetXIn() + pointData[i]->GetXOut()) / 2.;
            y = (pointData[i]->GetYIn() + pointData[i]->GetYOut()) / 2.;
            z = (pointData[i]->GetZIn() + pointData[i]->GetZOut()) / 2.;
            vpos.SetXYZ(x, y, z);

            // vpos = fRot * (vpos - fTrans);
            vpos = (vpos - fTrans);
            // vpos = fRot * (vpos);
            time = pointData[i]->GetTime() + rand->Gaus(0., tsigma);

            fiber = (int)std::round(vpos.X() / fiber_thickness);
            // Add hit data
            AddHitData(1,
                       vpos.X() * 10. + rand->Gaus(0., 0.1),
                       vpos.Y() * 10. + rand->Gaus(0., 0.1),
                       1,
                       time,
                       fiber,
                       1,
                       1,
                       1,
                       1);
        }
    }
    if (pointData)
        delete[] pointData;
    return;
}

// -----   Public method ReInit   ----------------------------------------------
InitStatus R3BFiberDigitizer::ReInit()
{
    SetParContainers();
    SetParameter();
    return kSUCCESS;
}
void R3BFiberDigitizer::Reset()
{
    if (fFiHits)
        fFiHits->Clear();
}

// -----   Private method AddR3BHitData  -------------------------------------------
R3BBunchedFiberHitData* R3BFiberDigitizer::AddHitData(Int_t detId,
                                                      Double_t x,
                                                      Double_t y,
                                                      Double_t eloss,
                                                      Double_t time,
                                                      Int_t a_fiber_id,
                                                      Double_t a_mapmt_time_ns,
                                                      Double_t a_spmt_time_ns,
                                                      Double_t a_mapmt_tot_ns,
                                                      Double_t a_spmt_tot_ns)
{
    // It fills the R3BBunchedFiberHitData
    TClonesArray& clref = *fFiHits;
    Int_t size = clref.GetEntriesFast();
    return new (clref[size]) R3BBunchedFiberHitData(
        detId, x, y, eloss, time, a_fiber_id, a_mapmt_time_ns, a_spmt_time_ns, a_mapmt_tot_ns, a_spmt_tot_ns);
}

ClassImp(R3BFiberDigitizer);
