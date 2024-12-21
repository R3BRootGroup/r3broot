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
// -----             R3BTttxMapped2Cal source file                     -----
// -----               Created 18/01/24  by XXX                        -----
// -------------------------------------------------------------------------

// ROOT headers
#include <TClonesArray.h>
#include <TMath.h>
#include <TRandom.h>
#include <iomanip>

// Fair headers
#include <FairLogger.h>
#include <FairRootManager.h>
#include <FairRunAna.h>
#include <FairRuntimeDb.h>

// R3B headers
#include "R3BLogger.h"
#include "R3BTttxCalData.h"
#include "R3BTttxMapped2Cal.h"
#include "R3BTttxMappedData.h"
#include "R3BTttxStripCalPar.h"

constexpr int ov_mask = 0x80000;
constexpr int pu_mask = 0x40000;

// R3BTttxMapped2Cal: Default Constructor --------------------------
R3BTttxMapped2Cal::R3BTttxMapped2Cal()
    : R3BTttxMapped2Cal("R3BTttxMapped2Cal", 1)
{
}

// R3BTttxMapped2CalPar: Standard Constructor --------------------------
R3BTttxMapped2Cal::R3BTttxMapped2Cal(const TString& name, int iVerbose)
    : FairTask(name, iVerbose)
{
}

// Virtual R3BTttxMapped2Cal: Destructor
R3BTttxMapped2Cal::~R3BTttxMapped2Cal()
{
    if (fTttxCalData)
        delete fTttxCalData;
}

void R3BTttxMapped2Cal::SetParContainers()
{
    // Parameter Containers
    FairRuntimeDb* rtdb = FairRuntimeDb::instance();
    R3BLOG_IF(error, rtdb == nullptr, "FairRuntimeDb not found");

    fCal_Par = dynamic_cast<R3BTttxStripCalPar*>(rtdb->getContainer("tttxStripCalPar"));
    if (!fCal_Par)
    {
        R3BLOG(error, "Couldn't get handle on tttxStripCalPar container");
    }
    else
    {
        R3BLOG(info, "tttxStripCalPar container open");
    }
}

void R3BTttxMapped2Cal::SetParameter()
{
    //--- Parameter Container ---
    NumDets = fCal_Par->GetNumDets();      // Number of Detectors
    NumStrips = fCal_Par->GetNumStrips();  // Number of Strips
    NumParams = fCal_Par->GetNumParsFit(); // Number of Parameters

    R3BLOG(info, "Nb detectors: " << NumDets);
    R3BLOG(info, "Nb strips: " << NumStrips);
    R3BLOG(info, "Nb parameters from pedestal fit: " << NumParams);

    CalParams = new TArrayF();
    int array_size = NumDets * NumStrips * NumParams;
    CalParams->Set(array_size);
    CalParams = const_cast<TArrayF*>(fCal_Par->GetStripCalParams()); // Array with the Cal parameters

    for (int d = 0; d < NumDets; d++)
    {
        int numdeadstrips = 0;
        for (int i = 0; i < NumStrips; i++)
        {
            if (fCal_Par->GetInUse(d + 1, i + 1) != 1)
                numdeadstrips++;
        }
        R3BLOG(info, "Nb of dead strips in TTTX detector " << d << ": " << numdeadstrips);
    }
}

// -----   Public method Init   --------------------------------------------
InitStatus R3BTttxMapped2Cal::Init()
{
    R3BLOG(info, "");

    // INPUT DATA
    FairRootManager* rootManager = FairRootManager::Instance();
    R3BLOG_IF(fatal, rootManager == nullptr, "FairRootManager not found");

    fTttxMapData = dynamic_cast<TClonesArray*>(rootManager->GetObject("tttxMappedData"));
    R3BLOG_IF(fatal, fTttxMapData == nullptr, "tttxMappedData not found");

    // OUTPUT: Calibrated data
    fTttxCalData = new TClonesArray("R3BTttxCalData");
    rootManager->Register("tttxCalData", "TTTX strip Cal", fTttxCalData, fOnline ? false : true);

    SetParameter();
    //
    fTttx.resize(NumDets);
    for (auto& val : fTttx)
    {
        val.ReInit(NumStrips);
    }
    return kSUCCESS;
}

// -----   Public method ReInit   ----------------------------------------------
InitStatus R3BTttxMapped2Cal::ReInit()
{
    SetParContainers();
    SetParameter();
    return kSUCCESS;
}

// -----   Public method Execution   --------------------------------------------
void R3BTttxMapped2Cal::Exec(Option_t* /*option*/)
{
    // Reset entries in output arrays, local arrays
    Reset();

    R3BLOG_IF(error, !fCal_Par, "NO Container Parameter");

    // Reading the Input -- Mapped Data --
    auto nHits = fTttxMapData->GetEntriesFast();
    if (!nHits)
        return;
    for (int ihit = 0; ihit < nHits; ihit++)
    {
        auto* hit = dynamic_cast<R3BTttxMappedData*>(fTttxMapData->At(ihit));
        auto idet = static_cast<int>(hit->GetDetID()) - 1;
        auto istrip = static_cast<int>(hit->GetStripID()) - 1;
        if (istrip < NumStrips)
        {
            fTttx.at(idet).GetStrip(istrip).SetRawValues(hit);
        }
        else if (istrip == fch_tref - 1)
        {
            fTttx.at(idet).GetTref().SetRawValues(hit);
        }
        else if (istrip == fch_trig - 1)
        {
            fTttx.at(idet).GetTrig().SetRawValues(hit);
        }
    }
    // End of reading data

    // Apply calibrations and addcaldata
    for (int idet = 0; idet < NumDets; idet++)
    {
        auto fTrig = fTttx.at(idet).GetTrig();
        R3BLOG(debug, "Detector" << idet + 1 << "Multiplicity of Trigger: " << fTrig.GetMult());
        R3BLOG_IF(
            warning, fTrig.GetMult() > 1, "Multiple Trigger hits: " << fTrig.GetMult() << ". Using the first hit.");
        for (int istrip = 0; istrip < NumStrips; istrip++)
        {
            if (fCal_Par->GetInUse(idet + 1, istrip + 1) != 1)
                continue;
            CalculateStrip(idet, istrip, fTttx.at(idet).GetStrip(istrip), fTrig);
        }
        CalculateStrip(idet, fch_tref, fTttx.at(idet).GetTref(), fTrig); // Store Tref data also
    }
    return;
}

// -----   Public method Reset   ------------------------------------------------
void R3BTttxMapped2Cal::Reset()
{
    R3BLOG(debug, "Clearing tttxCalData structure");
    if (fTttxCalData)
        fTttxCalData->Clear();
    for (auto& val : fTttx)
    {
        val.Init();
    }
}

void R3BTttxMapped2Cal::CalculateStrip(int idet, int istrip, CalStrip& fStrip, CalStrip& fTrig)
{
    for (int imult = 0; imult < fStrip.GetMult(); imult++)
    {
        if (fTrig.GetMult() == 0)
        {
            R3BLOG(warning, "No Trigger hit found. Skip this event.");
            return;
        }

        auto itime = fStrip.GetT(imult) - fTrig.GetT(0);
        auto dtime = fTimeResolution * static_cast<double>(itime);
        if (dtime < fTimeMin || dtime > fTimeMax)
            return;
        auto rawE = static_cast<int>(fStrip.GetE(imult));
        double energy = 0.;
        if ((ov_mask & rawE) != 0x0)
        {
            energy = TMath::Infinity();
        }
        else if ((pu_mask & rawE) != 0x0)
        {
            energy = std::nan("");
        }
        else if (rawE != 0 && istrip < NumStrips) // for the case of tref, rawE is 0
        {
            for (int power = 0; power < NumParams; power++)
            {
                energy += CalParams->GetAt((idet * NumStrips + istrip) * NumParams + power) *
                          TMath::Power(static_cast<double>(rawE), power);
            }
        }
        R3BLOG(debug,
               "idet = " << idet << ", istrip = " << istrip << ", time = " << dtime << " ns, energy = " << energy
                         << " keV");
        AddCalData(idet + 1, istrip + 1, dtime, energy);
    }
}

// -----   Private method AddCalData  --------------------------------------------
R3BTttxCalData* R3BTttxMapped2Cal::AddCalData(uint8_t detid, uint8_t stripid, double time, double energy)
{
    // It fills the R3BTttxCalData
    TClonesArray& clref = *fTttxCalData;
    int size = clref.GetEntriesFast();
    return new (clref[size]) R3BTttxCalData(detid, stripid, time, energy);
}

ClassImp(R3BTttxMapped2Cal)
