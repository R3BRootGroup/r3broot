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

// -------------------------------------------------------------------
// -----                     R3BAmsDigitizer                     -----
// -----           Created 09/09/22 by J.L. Rodriguez-Sanchez    -----
// -------------------------------------------------------------------

#include "R3BAmsDigitizer.h"
#include "FairRootManager.h"
#include "R3BAmsStripCalData.h"
#include "R3BLogger.h"
#include "R3BTraPoint.h"
#include "TClonesArray.h"
#include "TMath.h"
#include "TVector3.h"

R3BAmsDigitizer::R3BAmsDigitizer()
    : FairTask("R3BAmsDigitizer", false)
    , fPointData(NULL)
    , fCalData(NULL)
{
}

R3BAmsDigitizer::~R3BAmsDigitizer()
{
    R3BLOG(debug1, "");
    if (fPointData)
    {
        fPointData->Delete();
        delete fPointData;
    }
}

InitStatus R3BAmsDigitizer::Init()
{
    R3BLOG(info, "");

    FairRootManager* rootManager = FairRootManager::Instance();

    R3BLOG_IF(fatal, !rootManager, "FairRootManager not found");

    fPointData = dynamic_cast<TClonesArray*>(rootManager->GetObject("TraPoint"));
    if (!fPointData)
    {
        LOG(fatal) << "Init: No CrystalPoint CA";
        return kFATAL;
    }

    fCalData = new TClonesArray("R3BAmsStripCalData");
    rootManager->Register("AmsStripCalData", "AmsCal", fCalData, kTRUE);

    return kSUCCESS;
}

void R3BAmsDigitizer::Exec(Option_t* option)
{
    // Reset entries in output arrays, local arrays
    Reset();

    // Reading the Input -- Point data --
    Int_t nHits = fPointData->GetEntries();
    if (!nHits)
        return;

    auto pointData = new R3BTraPoint*[nHits];
    for (Int_t i = 0; i < nHits; i++)
    {
        pointData[i] = dynamic_cast<R3BTraPoint*>(fPointData->At(i));

        Int_t detid = pointData[i]->GetDetectorID();
        Int_t stripid = pointData[i]->GetDetCopyID() - 1;
        Int_t sideid = pointData[i]->GetDetCopyID() - 1 < 640 ? 0 : 1;
        if (sideid == 1)
            stripid -= 640;
        Double_t energy = pointData[i]->GetEnergyLoss();
        Int_t nStripCals = fCalData->GetEntriesFast();
        Bool_t existHit = 0;
        if (nStripCals == 0)
        {
            AddCal(detid, sideid, stripid, energy);
        }
        else
        {
            for (Int_t j = 0; j < nStripCals; j++)
            {
                if ((dynamic_cast<R3BAmsStripCalData*>(fCalData->At(j)))->GetStripId() == stripid &&
                    (dynamic_cast<R3BAmsStripCalData*>(fCalData->At(j)))->GetSideId() == sideid &&
                    (dynamic_cast<R3BAmsStripCalData*>(fCalData->At(j)))->GetDetId() == detid)
                {
                    (dynamic_cast<R3BAmsStripCalData*>(fCalData->At(j)))->AddMoreEnergy(energy);
                    existHit = 1; // to avoid the creation of a new StripCal Hit
                    break;
                }
            }
            if (!existHit)
                AddCal(detid, sideid, stripid, energy);
        }
        existHit = 0;
    }
}

void R3BAmsDigitizer::Reset()
{
    R3BLOG(debug1, "Clearing StripCalData Structure");
    if (fCalData)
        fCalData->Clear();
}

R3BAmsStripCalData* R3BAmsDigitizer::AddCal(Int_t detid, Int_t sideid, Int_t stripid, Double_t energy)
{
    TClonesArray& clref = *fCalData;
    Int_t size = clref.GetEntriesFast();

    R3BLOG_IF(info,
              fVerbose,
              "Adding StripCalData "
                  << " with unique identifier " << detid << " strip " << stripid << " side " << sideid
                  << " energy=" << energy);

    return new (clref[size]) R3BAmsStripCalData(detid, sideid, stripid, energy);
}

ClassImp(R3BAmsDigitizer);
