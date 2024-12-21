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

// ----------------------------------------------------------------------
// -----		  R3BMwpc1Mapped2Cal 			      -----
// -----          Created 15/10/19 by G. García Jiménez             -----
// -----  by modifying J.L. Rodriguez-Sanchez  classes for Mwpc2    -----
// ----------------------------------------------------------------------

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
#include "R3BLogger.h"
#include "R3BMwpc1CalPar.h"
#include "R3BMwpc1Mapped2Cal.h"
#include "R3BMwpcCalData.h"
#include "R3BMwpcMappedData.h"

// R3BMwpc1Mapped2Cal: Default Constructor --------------------------
R3BMwpc1Mapped2Cal::R3BMwpc1Mapped2Cal()
    : R3BMwpc1Mapped2Cal("R3BMwpc1Mapped2Cal", 1)
{
}

// R3BMwpc1Mapped2Cal: Standard Constructor --------------------------
R3BMwpc1Mapped2Cal::R3BMwpc1Mapped2Cal(const char* name, Int_t iVerbose)
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

// Virtual R3BMwpc1Mapped2Cal: Destructor
R3BMwpc1Mapped2Cal::~R3BMwpc1Mapped2Cal()
{
    R3BLOG(debug1, "Destructor");
    if (fMwpcCalDataCA)
        delete fMwpcCalDataCA;
}

void R3BMwpc1Mapped2Cal::SetParContainers()
{
    // Parameter Container
    // Reading PadCalPar from FairRuntimeDb
    FairRuntimeDb* rtdb = FairRuntimeDb::instance();
    R3BLOG_IF(error, !rtdb, "FairRuntimeDb not found");

    fCal_Par = dynamic_cast<R3BMwpc1CalPar*>(rtdb->getContainer("mwpc1CalPar"));
    if (!fCal_Par)
    {
        R3BLOG(error, "Couldn't get handle on mwpc1CalPar container.");
    }
    else
    {
        R3BLOG(info, "mwpc1CalPar container found.");
    }
    return;
}

void R3BMwpc1Mapped2Cal::SetParameter()
{
    //--- Parameter Container ---
    NumPadX = fCal_Par->GetNumPadsX();           // Number of Pads in X
    NumPadY = fCal_Par->GetNumPadsY();           // Number of Pads in Y
    NumParams = fCal_Par->GetNumParametersFit(); // Number of parameters in the Fit

    R3BLOG(info, "NumPadX: " << NumPadX);
    R3BLOG(info, "NumPadY: " << NumPadY);
    R3BLOG(info, "Number of fit parameters: " << NumParams);

    CalParams = new TArrayF();
    Int_t array_size = (NumPadX + NumPadY) * NumParams;
    CalParams->Set(array_size);
    CalParams = fCal_Par->GetPadCalParams(); // Array with the Cal parameters
    return;
}

// -----   Public method Init   --------------------------------------------
InitStatus R3BMwpc1Mapped2Cal::Init()
{
    R3BLOG(info, "");

    // INPUT DATA
    FairRootManager* rootManager = FairRootManager::Instance();
    if (!rootManager)
    {
        R3BLOG(fatal, "FairRootManager not found");
        return kFATAL;
    }

    fMwpcMappedDataCA = dynamic_cast<TClonesArray*>(rootManager->GetObject("Mwpc1MappedData"));
    if (!fMwpcMappedDataCA)
    {
        R3BLOG(fatal, "Mwpc1MappedData not found");
        return kFATAL;
    }

    // OUTPUT DATA
    // Calibrated data
    fMwpcCalDataCA = new TClonesArray("R3BMwpcCalData");
    rootManager->Register("Mwpc1CalData", "MWPC1 Cal", fMwpcCalDataCA, !fOnline);
    fMwpcCalDataCA->Clear();

    SetParameter();
    return kSUCCESS;
}

// -----   Public method ReInit   ----------------------------------------------
InitStatus R3BMwpc1Mapped2Cal::ReInit()
{
    SetParContainers();
    SetParameter();
    return kSUCCESS;
}

// -----   Public method Execution   --------------------------------------------
void R3BMwpc1Mapped2Cal::Exec(Option_t* option)
{
    // Reset entries in output arrays, local arrays
    Reset();

    // Reading the Input -- Mapped Data --
    Int_t nHits = fMwpcMappedDataCA->GetEntriesFast();
    if (nHits > (NumPadX + NumPadY) && nHits > 0)
    {
        R3BLOG(warn, "nHits>(NumPadX+NumPadY)");
    }
    if (nHits == 0)
        return;

    R3BMwpcMappedData** mappedData;
    mappedData = new R3BMwpcMappedData*[nHits];
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
            R3BLOG(error, "Plane " << planeId << " does not exist in MWPC1");

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

// -----   Public method Reset   ------------------------------------------------
void R3BMwpc1Mapped2Cal::Reset()
{
    R3BLOG(debug1, "Clearing Mwpc1CalData Structure");
    if (fMwpcCalDataCA)
        fMwpcCalDataCA->Clear();
}

// -----   Private method AddCalData  --------------------------------------------
R3BMwpcCalData* R3BMwpc1Mapped2Cal::AddCalData(Int_t plane, Int_t pad, Float_t charge)
{
    // It fills the R3BMwpcCalData
    TClonesArray& clref = *fMwpcCalDataCA;
    Int_t size = clref.GetEntriesFast();
    return new (clref[size]) R3BMwpcCalData(plane, pad, charge);
}

ClassImp(R3BMwpc1Mapped2Cal);
