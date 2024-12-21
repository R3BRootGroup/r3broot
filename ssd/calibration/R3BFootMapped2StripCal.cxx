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

// -----------------------------------------------------------------
// -----         R3BFootMapped2StripCal source file            -----
// -----      Created 08/09/21 by J.L. Rodriguez-Sanchez       -----
// -----------------------------------------------------------------

// ROOT headers
#include <TClonesArray.h>
#include <TMath.h>

// FAIR headers
#include <FairLogger.h>
#include <FairRootManager.h>
#include <FairRunAna.h>
#include <FairRuntimeDb.h>

// FOOT headers
#include "R3BFootCalData.h"
#include "R3BFootCalPar.h"
#include "R3BFootMapped2StripCal.h"
#include "R3BFootMappedData.h"
#include "R3BLogger.h"

// R3BFootMapped2StripCal: Default Constructor --------------------------
R3BFootMapped2StripCal::R3BFootMapped2StripCal()
    : R3BFootMapped2StripCal("R3BFootMapped2StripCal", 1)
{
}

// R3BFootMapped2StripCalPar: Standard Constructor --------------------------
R3BFootMapped2StripCal::R3BFootMapped2StripCal(const TString& name, Int_t iVerbose)
    : FairTask(name, iVerbose)
{
}
// Virtual R3BFootMapped2StripCal: Destructor
R3BFootMapped2StripCal::~R3BFootMapped2StripCal()
{
    R3BLOG(debug1, "");
    if (fFootCalData)
        delete fFootCalData;
}

void R3BFootMapped2StripCal::SetParContainers()
{
    // Parameter Container
    // Reading footCalPar from FairRuntimeDb
    FairRuntimeDb* rtdb = FairRuntimeDb::instance();
    R3BLOG_IF(error, !rtdb, "FairRuntimeDb not found");

    fCal_Par = dynamic_cast<R3BFootCalPar*>(rtdb->getContainer("footCalPar"));
    if (!fCal_Par)
    {
        R3BLOG(error, "footCalPar container not found");
    }
    else
    {
        R3BLOG(info, "footCalPar found");
    }
}

void R3BFootMapped2StripCal::SetParameter()
{
    //--- Parameter Container ---
    if (fCal_Par == nullptr)
    {
        R3BLOG(error, "footCalPar not found");
        return;
    }

    fNDets = fCal_Par->GetNumDets();        // Number of Detectors
    fNStrip = fCal_Par->GetNumStrips();     // Number of Strips
    fNumParams = fCal_Par->GetNumParsFit(); // Number of Parameters

    R3BLOG(info, "Nb detectors: " << fNDets);
    R3BLOG(info, "Nb strips: " << fNStrip);
    R3BLOG(info, "Nb parameters from pedestal fit: " << fNumParams);

    CalParams = new TArrayF();
    Int_t array_size = fNDets * fNStrip * fNumParams;
    CalParams->Set(array_size);
    CalParams = fCal_Par->GetStripCalParams(); // Array with the Cal parameters

    // Count the number of dead strips per AMS detector
    for (Int_t d = 0; d < fNDets; d++)
    {
        Int_t numdeadstrips = 0;
        for (Int_t i = 0; i < fNStrip; i++)
        {
            if (CalParams->GetAt(fNumParams * i + 1 + fNStrip * d * fNumParams) == -1)
                numdeadstrips++;
        }
        R3BLOG(info, "Nb of dead strips in FOOT detector " << d + 1 << ": " << numdeadstrips);
    }
}

// -----   Public method Init   --------------------------------------------
InitStatus R3BFootMapped2StripCal::Init()
{
    R3BLOG(info, "");

    FairRootManager* rmg = FairRootManager::Instance();
    R3BLOG_IF(fatal, rmg == nullptr, "FairRootManager not found");

    // INPUT DATA
    fFootMappedData = dynamic_cast<TClonesArray*>(rmg->GetObject("FootMappedData"));
    R3BLOG_IF(fatal, fFootMappedData == nullptr, "FootMappedData not found");

    // OUTPUT DATA
    fFootCalData = new TClonesArray("R3BFootCalData");
    rmg->Register("FootCalData", "FOOT strip Cal", fFootCalData, !fOnline);
    fFootCalData->Clear();

    SetParameter();
    return kSUCCESS;
}

// -----   Public method ReInit   ----------------------------------------------
InitStatus R3BFootMapped2StripCal::ReInit()
{
    SetParContainers();
    SetParameter();
    return kSUCCESS;
}

// -----   Public method Execution   --------------------------------------------
void R3BFootMapped2StripCal::Exec(Option_t* /*option*/)
{
    // Reset entries in output arrays, local arrays
    Reset();

    // Reading the Input -- Mapped Data --
    int nHits = fFootMappedData->GetEntriesFast();
    if (nHits == 0)
        return;

    R3BFootMappedData** mappedData = new R3BFootMappedData*[nHits];
    int8_t detId;
    int16_t stripId;
    double energy = 0.;
    double pedestal = 0.;
    double sigma = 0.;

    const auto detindex = fCal_Par->GetNumDets();
    const auto asics_index = 10;

    int StripNAve[detindex];                     // Number of Strips to compute average for each detector
    double Ave[detindex];                        // Average ADC value for each detector
    double StripNAveASIC[detindex][asics_index]; // Number of Strips to compute average for each ASICS
    // double StripNAveASIC_fine[detindex][asics_index];
    double AveASIC[detindex][asics_index]; // Average ADC value for each ASICS
    // double asic_offset_fine[detindex][asics_index];
    int StripCounter[detindex]; // Counter to disregard events with baseline jumps
    for (int i = 0; i < detindex; i++)
    {
        StripCounter[i] = 0;
        Ave[i] = 0.;
        StripNAve[i] = 0;
        for (int j = 0; j < asics_index; j++)
        {
            AveASIC[i][j] = 0.;
            StripNAveASIC[i][j] = 0;
        }
    }

    // Pedestal substraction
    for (Int_t i = 0; i < nHits; i++)
    {
        mappedData[i] = dynamic_cast<R3BFootMappedData*>(fFootMappedData->At(i));
        detId = mappedData[i]->GetDetId() - 1;
        stripId = mappedData[i]->GetStripId() - 1;

        if (CalParams)
        {
            pedestal = CalParams->GetAt(fNumParams * stripId + detId * fNumParams * fNStrip);
            sigma = CalParams->GetAt(fNumParams * stripId + 1 + detId * fNumParams * fNStrip);
        }

        energy = mappedData[i]->GetEnergy() - pedestal;

        // only use strips with signal below threshold for correction
        if (energy < fTimesSigma * sigma)
        {
            Ave[detId] += energy;
            StripNAve[detId]++;
        }
    }

    for (Int_t i = 0; i < fCal_Par->GetNumDets(); i++)
    {
        Ave[i] = Ave[i] / StripNAve[i];
    }

    // Average correction
    for (Int_t i = 0; i < nHits; i++)
    {
        detId = mappedData[i]->GetDetId() - 1;
        stripId = mappedData[i]->GetStripId() - 1;

        if (CalParams != nullptr)
        {
            pedestal = CalParams->GetAt(fNumParams * stripId + detId * fNumParams * fNStrip);
            sigma = CalParams->GetAt(fNumParams * stripId + 1 + detId * fNumParams * fNStrip);
        }

        energy = mappedData[i]->GetEnergy() - pedestal - Ave[detId];

        // only use strips with signal below threshold for correction
        if (energy < fTimesSigma * sigma)
        {
            Int_t ASIC = (Double_t)stripId / 64.;
            AveASIC[detId][ASIC] += energy;
            StripNAveASIC[detId][ASIC]++;
        }
    }

    for (Int_t i = 0; i < fCal_Par->GetNumDets(); i++)
    {
        for (Int_t j = 0; j < 10; j++)
        {
            AveASIC[i][j] = AveASIC[i][j] / StripNAveASIC[i][j];
        }
    }

    // ASIC Average correction
    for (Int_t i = 0; i < nHits; i++)
    {

        detId = mappedData[i]->GetDetId() - 1;
        stripId = mappedData[i]->GetStripId() - 1;

        if (CalParams)
        {
            pedestal = CalParams->GetAt(fNumParams * stripId + detId * fNumParams * fNStrip);
            sigma = CalParams->GetAt(fNumParams * stripId + 1 + detId * fNumParams * fNStrip);
        }

        auto ASIC1 = static_cast<int>(stripId / 64);

        // energy = mappedData[i]->GetEnergy() - pedestal - Ave[detId] - AveASIC[detId][ASIC1];

        energy = mappedData[i]->GetEnergy() - pedestal - fTimesSigma * sigma - Ave[detId] - AveASIC[detId][ASIC1];

        if (energy > 0. && pedestal != -1)
        {
            StripCounter[detId]++;
        }
    }

    for (Int_t i = 0; i < nHits; i++)
    {

        detId = mappedData[i]->GetDetId() - 1;
        stripId = mappedData[i]->GetStripId() - 1;

        if (CalParams)
        {
            pedestal = CalParams->GetAt(fNumParams * stripId + detId * fNumParams * fNStrip);
            sigma = CalParams->GetAt(fNumParams * stripId + 1 + detId * fNumParams * fNStrip);
        }

        auto ASIC3 = static_cast<int>(stripId / 64);

        // energy = mappedData[i]->GetEnergy() - pedestal - Ave[detId] - AveASIC[detId][ASIC3];
        energy = mappedData[i]->GetEnergy() - pedestal - fTimesSigma * sigma - Ave[detId] - AveASIC[detId][ASIC3];

        if (energy > 0 && pedestal != -1 && StripCounter[detId] < fNStrip)
        {
            AddCalData(detId + 1, stripId + 1, energy);
        }
    }
    if (mappedData)
        delete[] mappedData;
    return;
}

// -----   Public method Reset   ------------------------------------------------
void R3BFootMapped2StripCal::Reset()
{
    R3BLOG(debug1, "Clearing StripCalData Structure");
    if (fFootCalData)
        fFootCalData->Clear();
}

// -----   Private method AddCalData  --------------------------------------------
R3BFootCalData* R3BFootMapped2StripCal::AddCalData(int8_t detid, int16_t stripid, double energy)
{
    // It fills the R3BFootCalData
    TClonesArray& clref = *fFootCalData;
    Int_t size = clref.GetEntriesFast();
    return new (clref[size]) R3BFootCalData(detid, stripid, energy);
}

ClassImp(R3BFootMapped2StripCal)
