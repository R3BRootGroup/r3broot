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

// ----------------------------------------------------------------
// -----                    R3BPspxCal2Hit                    -----
// -----	    created 09-03-2016 by I. Syndikus          -----
// -----             Modified  Dec 2019  by M. Holl	        -----
// ----------------------------------------------------------------

// Fair headers
#include "FairLogger.h"
#include "FairRootManager.h"
#include "FairRunOnline.h"
#include "FairRuntimeDb.h"

// PSP headers
#include "R3BEventHeader.h"
#include "R3BPspxCal2Hit.h"
#include "R3BPspxCalData.h"
#include "R3BPspxHitData.h"
#include "R3BPspxHitPar.h"

#include "TClonesArray.h"
#include <cmath>
#include <iostream>
#include <limits>

R3BPspxCal2Hit::R3BPspxCal2Hit()
    : FairTask("PspxCal2Hit", 1)
    , fCalItems()
    , fHitItems()
    , fOnline(kFALSE)
{
}

R3BPspxCal2Hit::R3BPspxCal2Hit(const char* name, Int_t iVerbose)
    : FairTask(name, iVerbose)
    , fCalItems()
    , fHitItems()
    , fOnline(kFALSE)
{
}

R3BPspxCal2Hit::~R3BPspxCal2Hit()
{

    LOG(info) << "R3BPspxCal2Hit: Delete instance";
    for (Int_t i = 0; i < fCalItems.size(); i++)
    {
        delete fCalItems[i];
    }
    for (Int_t i = 0; i < fHitItems.size(); i++)
    {
        delete fHitItems[i];
    }
}

// -- SetParameters
// -- Read calibration parameters
void R3BPspxCal2Hit::SetParameters()
{

    LOG(info) << "In R3BPspxCal2Hit::SetParameters()";
    //--- Parameter Container ---
    Int_t nDet = fHitPar->GetNumDetectors(); // Number of Detectors/Faces
    orientation.resize(nDet);
    detSize.resize(nDet);
    posOffset.resize(nDet);
    posSlope.resize(nDet);
    LOG(info) << "Position Parameters";
    Int_t parOffset = 1; // 1 "header" parameter.
    for (Int_t d = 0; d < nDet; d++)
    {
        TArrayF par = fHitPar->GetHitPosPar(); // Array with the parameters
        orientation[d] = par.At(parOffset + 2);
        detSize[d] = par.At(parOffset + 3);
        posOffset[d] = par.At(parOffset + 4);
        posSlope[d] = par.At(parOffset + 5);
        parOffset += 6; // move to next line in parameter file.
        LOG(info) << "Det: " << d << "\torientation: " << orientation[d] << "\tdetSize: " << detSize[d]
                  << "\toffset: " << posOffset[d] << "\tslope: " << posSlope[d];
    }

    eOffset.resize(nDet);
    eGain.resize(nDet);
    eRange.resize(nDet);
    LOG(info) << "Energy Parameters";
    parOffset = 1; // 1 "header" parameter.
    for (Int_t d = 0; d < nDet; d++)
    {
        TArrayF par = fHitPar->GetHitEPar(); // Array with the parameters
        eOffset[d] = par.At(parOffset + 2);
        eGain[d] = par.At(parOffset + 3);
        eRange[d] = par.At(parOffset + 4);
        parOffset += 5; // move to next line in parameter file.
        LOG(info) << "Det: " << d << "\toffset: " << eOffset[d] << "\tslope: " << eGain[d] << "\trange: " << eRange[d];
    }
}

// -- Init
// -- Initialize output data. Read input data.
InitStatus R3BPspxCal2Hit::Init()
{
    FairRootManager* fMan = FairRootManager::Instance();
    if (!fMan)
    {
        LOG(error) << "R3BPspxCal2Hit::Init() Root-manager not found.";
        return kFATAL;
    }

    // R3BEventHeader for trigger information, if needed!
    fHeader = dynamic_cast<R3BEventHeader*>(fMan->GetObject("R3BEventHeader"));

    const char xy[2] = { 'x', 'y' }; // orientation of detector face
    // Figure out how many detectors were registered by the reader
    for (Int_t d = 0;; d++)
    {
        TClonesArray* tmp[2];

        for (Int_t f = 0; f < 2; f++)
        {
            tmp[f] = dynamic_cast<TClonesArray*>(
                fMan->GetObject(Form("Pspx%d_%cCal", d + 1, xy[f]))); // = branch name in TTree
        }
        if (tmp[0] == NULL && tmp[1] == NULL)
        {
            if (d == 0)
            {
                printf("Couldn't get handle on PSPX cal items\n");
                return kFATAL;
            }
            break;
        }
        for (Int_t f = 0; f < 2; f++)
        {
            fCalItems.push_back(tmp[f]);
            fHitItems.push_back(new TClonesArray("R3BPspxHitData"));
            if (!fOnline)
            {
                fMan->Register(
                    Form("Pspx%d_%cHit", d + 1, xy[f]), Form("Pspx%d_%c", d + 1, xy[f]), fHitItems.back(), kTRUE);
            }
            else
            {
                fMan->Register(
                    Form("Pspx%d_%cHit", d + 1, xy[f]), Form("Pspx%d_%c", d + 1, xy[f]), fHitItems.back(), kFALSE);
            }
        }
    }

    SetParameters();
    return kSUCCESS;
}

// -- Initialisation
// -- Initialize/Read parameter file for conversion.
void R3BPspxCal2Hit::SetParContainers()
{
    FairRuntimeDb* rtdb = FairRuntimeDb::instance();
    if (!rtdb)
    {
        LOG(error) << "R3BPspxCal2Hit::FairRuntimeDb not opened!";
        return;
    }
    else
    {
        LOG(info) << "R3BPspxCal2Hit::SetParContainers()";
    }

    fHitPar = dynamic_cast<R3BPspxHitPar*>(rtdb->getContainer("R3BPspxHitPar"));
    if (!fHitPar)
    {
        LOG(error) << "R3BPspxCal2Hit::Could not get access to R3BPspxHitPar-Container.";
        return;
    }

    fHitPar->printParams();
}
// --------------------------------------------------------------------

// -- ReInit
// -- Initialize/Read parameter file for conversion.
InitStatus R3BPspxCal2Hit::ReInit()
{
    LOG(info) << " R3BPspxCal2Hit :: ReInit() ";
    /*
        fHitPar = (R3BPspxHitPar*)FairRuntimeDb::instance()->getContainer("R3BPspxHitPar");

        if (!fHitPar)
        {
            LOG(error) << "Could not get access to R3BPspxHitPar-Container.";
            return kFATAL;
        }*/

    SetParContainers();
    SetParameters();
    return kSUCCESS;
}

void R3BPspxCal2Hit::Exec(Option_t* option)
{
    const char xy[2] = { 'x', 'y' }; // orientation of detector face
    for (Int_t d = 0; d < fCalItems.size(); d++)
    {
        if (!fCalItems[d])
        {
            for (Int_t f = 0; f < 2; f++)
            {
                printf("Cannot access PSPX%d_%c cal items\n", (d / 2) + 1, xy[f]);
            }
            return;
        }
        Int_t nCal = fCalItems[d]->GetEntries();

        for (Int_t i = 0; i < nCal; i++)
        {

            R3BPspxCalData* calData = dynamic_cast<R3BPspxCalData*>(fCalItems[d]->At(i)); // get cal level event

            Float_t energy = calData->GetEnergy() * eGain[d] + eOffset[d]; // convert energy to MeV
            Float_t pos =
                calData->GetPos() * orientation[d] /** detSize[d]*/; // convert position to mm, flip axis if necessary
            pos = pos * posSlope[d] + posOffset[d];                  // correct position for detector offset and tilt
            new ((*fHitItems[d])[fHitItems[d]->GetEntriesFast()])
                R3BPspxHitData(energy, pos); // register hit level event
        }
    }
}

void R3BPspxCal2Hit::FinishEvent()
{
    for (Int_t i = 0; i < fCalItems.size(); i++)
    {
        fCalItems[i]->Clear();
    }
    for (Int_t i = 0; i < fHitItems.size(); i++)
    {
        fHitItems[i]->Clear();
    }
}

void R3BPspxCal2Hit::FinishTask() {}

ClassImp(R3BPspxCal2Hit)
