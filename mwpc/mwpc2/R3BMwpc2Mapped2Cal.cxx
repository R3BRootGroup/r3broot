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
// -----            R3BMwpc2Mapped2Cal source file            -----
// -----     Created 10/10/19 by J.L. Rodriguez-Sanchez       -----
// ----------------------------------------------------------------

// ROOT headers
#include "TClonesArray.h"
#include "TMath.h"
#include "TRandom.h"

// Fair headers
#include "FairLogger.h"
#include "FairRootManager.h"
#include "FairRunAna.h"
#include "FairRuntimeDb.h"

// MWPC headers
#include "R3BMwpc2CalPar.h"
#include "R3BMwpc2Mapped2Cal.h"
#include "R3BMwpcCalData.h"
#include "R3BMwpcMappedData.h"

// R3BMwpc2Mapped2Cal: Default Constructor --------------------------
R3BMwpc2Mapped2Cal::R3BMwpc2Mapped2Cal()
    : R3BMwpc2Mapped2Cal("R3BMwpc2Mapped2Cal", 1)
{
}

// R3BMwpc2Mapped2Cal: Standard Constructor --------------------------
R3BMwpc2Mapped2Cal::R3BMwpc2Mapped2Cal(const char* name, Int_t iVerbose)
    : FairTask(name, iVerbose)
    , NumPadX(0)
    , NumPadY(0)
    , NumParams(0)
    , CalParams(NULL)
    , fCal_Par(NULL)
    , fMwpcMappedDataCA(NULL)
    , fMwpcCalDataCA(NULL)
    , fOnline(kFALSE)
{
}

// Virtual R3BMwpc2Mapped2Cal: Destructor
R3BMwpc2Mapped2Cal::~R3BMwpc2Mapped2Cal()
{
    LOG(info) << "R3BMwpc2Mapped2Cal: Delete instance";
    if (fMwpcCalDataCA)
        delete fMwpcCalDataCA;
}

void R3BMwpc2Mapped2Cal::SetParContainers()
{
    // Parameter Container
    // Reading PadCalPar from FairRuntimeDb
    FairRuntimeDb* rtdb = FairRuntimeDb::instance();
    if (!rtdb)
    {
        LOG(error) << "FairRuntimeDb not opened!";
    }

    fCal_Par = dynamic_cast<R3BMwpc2CalPar*>(rtdb->getContainer("mwpc2CalPar"));
    if (!fCal_Par)
    {
        LOG(error) << "R3BMwpc2Mapped2Cal::Init() Couldn't get handle on mwpc2CalPar container";
    }
    else
    {
        LOG(info) << "R3BMwpc2Mapped2Cal:: mwpc2CalPar container open";
    }
    return;
}

void R3BMwpc2Mapped2Cal::SetParameter()
{
    //--- Parameter Container ---
    NumPadX = fCal_Par->GetNumPadsX();           // Number of Pads in X
    NumPadY = fCal_Par->GetNumPadsY();           // Number of Pads in Y
    NumParams = fCal_Par->GetNumParametersFit(); // Number of parameters in the Fit

    LOG(info) << "R3BMwpc2Mapped2Cal: NumPadX: " << NumPadX;
    LOG(info) << "R3BMwpc2Mapped2Cal: NumPadY: " << NumPadY;
    LOG(info) << "R3BMwpc2Mapped2Cal: Number of fit parameters: " << NumParams;

    CalParams = new TArrayF();
    Int_t array_size = (NumPadX + NumPadY) * NumParams;
    CalParams->Set(array_size);
    CalParams = fCal_Par->GetPadCalParams(); // Array with the Cal parameters
    return;
}

// -----   Public method Init   --------------------------------------------
InitStatus R3BMwpc2Mapped2Cal::Init()
{
    LOG(info) << "R3BMwpc2Mapped2Cal::Init()";

    // INPUT DATA
    FairRootManager* rootManager = FairRootManager::Instance();
    if (!rootManager)
    {
        return kFATAL;
    }

    fMwpcMappedDataCA = dynamic_cast<TClonesArray*>(rootManager->GetObject("Mwpc2MappedData"));
    if (!fMwpcMappedDataCA)
    {
        return kFATAL;
    }

    // OUTPUT DATA
    // Calibrated data
    fMwpcCalDataCA = new TClonesArray("R3BMwpcCalData", 10);
    rootManager->Register("Mwpc2CalData", "MWPC2 Cal", fMwpcCalDataCA, !fOnline);

    SetParameter();
    return kSUCCESS;
}

// -----   Public method ReInit   ----------------------------------------------
InitStatus R3BMwpc2Mapped2Cal::ReInit()
{
    SetParContainers();
    SetParameter();
    return kSUCCESS;
}

// -----   Public method Execution   --------------------------------------------
void R3BMwpc2Mapped2Cal::Exec(Option_t* option)
{
    // Reset entries in output arrays, local arrays
    Reset();

    if (!fCal_Par)
    {
        LOG(error) << "NO Container Parameter!, pedestals will be set to zero";
    }

    // Reading the Input -- Mapped Data --
    Int_t nHits = fMwpcMappedDataCA->GetEntriesFast();
    if (nHits > (NumPadX + NumPadY))
        LOG(warn) << "R3BMwpc2Mapped2Cal: nHits>(NumPadX+NumPadY)";
    if (nHits == 0)
        return;

    R3BMwpcMappedData** mappedData = new R3BMwpcMappedData*[nHits];
    Int_t planeId = 0;
    Int_t padId = 0;
    Float_t charge = 0.0;
    Float_t pedestal = 0.0;
    Int_t nbpad = 0;

    for (Int_t i = 0; i < nHits; i++)
    {
        mappedData[i] = dynamic_cast<R3BMwpcMappedData*>(fMwpcMappedDataCA->At(i));
        planeId = mappedData[i]->GetPlane();
        padId = mappedData[i]->GetPad() - 1;
        if (planeId == 1) // X
            nbpad = padId * NumParams;
        else if (planeId == 2) // X
            nbpad = (padId + NumPadX / 2) * NumParams;
        else if (planeId == 3) // Y
            nbpad = (padId + NumPadX) * NumParams;
        else
            LOG(error) << "Plane " << planeId << " does not exist in MWPC2";

        pedestal = CalParams->GetAt(nbpad);
        charge = mappedData[i]->GetQ() - pedestal;

        // We accept the hit if the charge is larger than zero
        if (charge > 0)
        {
            AddCalData(planeId, padId + 1, charge);
        }
    }
    if (mappedData)
        delete[] mappedData;
    return;
}

// -----   Protected method Finish   --------------------------------------------
void R3BMwpc2Mapped2Cal::Finish() {}

// -----   Public method Reset   ------------------------------------------------
void R3BMwpc2Mapped2Cal::Reset()
{
    LOG(debug) << "Clearing Mwpc2CalData Structure";
    if (fMwpcCalDataCA)
        fMwpcCalDataCA->Clear();
}

// -----   Private method AddCalData  --------------------------------------------
R3BMwpcCalData* R3BMwpc2Mapped2Cal::AddCalData(Int_t plane, Int_t pad, Float_t charge)
{
    // It fills the R3BMwpcCalData
    TClonesArray& clref = *fMwpcCalDataCA;
    Int_t size = clref.GetEntriesFast();
    return new (clref[size]) R3BMwpcCalData(plane, pad, charge);
}

ClassImp(R3BMwpc2Mapped2Cal);
