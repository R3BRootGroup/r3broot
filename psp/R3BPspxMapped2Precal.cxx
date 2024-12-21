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
// -----                   R3BPspxMapped2Precal               -----
// -----            Created  13-03-2017 by I. Syndikus        -----
// -----              Modified  Dec 2019  by M. Holl          -----
// ----------------------------------------------------------------

// Fair headers
#include "FairLogger.h"
#include "FairRootManager.h"
#include "FairRunAna.h"
#include "FairRunOnline.h"
#include "FairRuntimeDb.h"

// PSP headers
#include "R3BEventHeader.h"
#include "R3BPspxMapped2Precal.h"
#include "R3BPspxMappedData.h"
#include "R3BPspxPrecalData.h"
#include "R3BPspxPrecalPar.h"

#include "TClonesArray.h"
#include <iostream>
#include <limits>

R3BPspxMapped2Precal::R3BPspxMapped2Precal()
    : FairTask("PspxMapped2Precal", 1)
    , fMappedItems()
    , fPrecalItems()
    , fOnline(kFALSE)
{
}

R3BPspxMapped2Precal::R3BPspxMapped2Precal(const char* name, Int_t iVerbose)
    : FairTask(name, iVerbose)
    , fMappedItems()
    , fPrecalItems()
    , fOnline(kFALSE)
{
}

R3BPspxMapped2Precal::~R3BPspxMapped2Precal()
{

    LOG(info) << "R3BPspxMapped2Precal: Delete instance";
    for (Int_t i = 0; i < fMappedItems.size(); i++)
    {
        delete fMappedItems[i];
    }
    for (Int_t i = 0; i < fPrecalItems.size(); i++)
    {
        delete fPrecalItems[i];
    }
}

void R3BPspxMapped2Precal::SetParameters()
{
    LOG(info) << "R3BPspxMapped2Precal::SetParameters()";
    //--- Parameter Container ---
    Int_t nDet = fPrecalPar->GetNumDetectors(); // Number of Detectors/Faces
    LOG(info) << nDet;
    gain.resize(nDet);
    threshold1.resize(nDet);
    threshold2.resize(nDet);
    for (Int_t d = 0; d < nDet; d++)
    {
        Int_t nStrips = fPrecalPar->GetNumStrips().At(d); // Number of Strips
        gain[d].resize(nStrips);
        threshold1[d].resize(nStrips);
        threshold2[d].resize(nStrips);
        Int_t parOffset =
            d * nStrips * 4 +
            (d + 1) * 3; // Position in parameter list. 4 parameters per strip + 3 "header" parameters per detector.
        for (Int_t s = 0; s < nStrips; s++)
        {
            TArrayF par = fPrecalPar->GetPrecalPar(); // Array with the parameters
            gain[d][s] = par.At(parOffset + 1);
            threshold1[d][s] = par.At(parOffset + 2);
            threshold2[d][s] = par.At(parOffset + 3);
            parOffset += 4; // move to next line in parameter file.
            LOG(info) << "Det: " << d << "\tstr: " << s << "\tgain: " << gain[d][s] << "\tthr1: " << threshold1[d][s]
                      << "\tthr2: " << threshold2[d][s];
        }
    }
}

InitStatus R3BPspxMapped2Precal::Init()
{
    /**
     * Initialize output data. Read input data and parameters.
     * The parameters get saved in dedicated arrays.
     * Print parameters, if verbosity is set to info.
     */

    LOG(info) << "R3BPspxMapped2Precal::Init()";
    FairRootManager* fMan = FairRootManager::Instance();
    if (!fMan)
    {
        LOG(error) << "R3BPspxMapped2Precal::Init() Root-manager not found.";
        return kFATAL;
    }

    // R3BEventHeader for trigger information, if needed!
    fHeader = dynamic_cast<R3BEventHeader*>(fMan->GetObject("R3BEventHeader"));

    const char xy[2] = { 'x', 'y' }; // orientation of detector face
    // Figure out how many detectors were registered by the reader
    for (Int_t d = 0; d < fPrecalPar->GetNumDetectors(); d++)
    {
        TClonesArray* tmp[2];

        for (Int_t f = 0; f < 2; f++)
        {
            tmp[f] = dynamic_cast<TClonesArray*>(
                fMan->GetObject(Form("Pspx%d_%cMapped", d + 1, xy[f]))); // = branch name in TTree
        }
        if (tmp[0] == NULL && tmp[1] == NULL)
        {
            if (d == 0)
            {
                LOG(error) << "R3BPspxMapped2Precal::Init() Couldn't get handle on PSPX-mapped items.";
                return kFATAL;
            }
            break;
        }
        for (Int_t f = 0; f < 2; f++)
        {
            fMappedItems.push_back(tmp[f]);
            fPrecalItems.push_back(new TClonesArray("R3BPspxPrecalData"));
            if (!fOnline)
            {
                fMan->Register(
                    Form("Pspx%d_%cPrecal", d + 1, xy[f]), Form("Pspx%d_%c", d + 1, xy[f]), fPrecalItems.back(), kTRUE);
            }
            else
            {
                fMan->Register(Form("Pspx%d_%cPrecal", d + 1, xy[f]),
                               Form("Pspx%d_%c", d + 1, xy[f]),
                               fPrecalItems.back(),
                               kFALSE);
            }
        }
    }

    SetParameters();
    return kSUCCESS;
}

void R3BPspxMapped2Precal::SetParContainers()
{
    /**
     * Initialize/Reads parameter file for conversion.
     */

    FairRuntimeDb* rtdb = FairRuntimeDb::instance();
    if (!rtdb)
    {
        LOG(error) << "R3BPspxMapped2Precal::FairRuntimeDb not opened!";
        return;
    }
    else
    {
        LOG(info) << "R3BPspxMapped2Precal::SetParContainers()";
    }

    fPrecalPar = dynamic_cast<R3BPspxPrecalPar*>(rtdb->getContainer("R3BPspxPrecalPar"));
    if (!fPrecalPar)
    {
        LOG(error) << "R3BPspxMapped2Precal::Could not get access to R3BPspxPrecalPar-Container.";
        return;
    }

    fPrecalPar->printParams();
}

InitStatus R3BPspxMapped2Precal::ReInit()
{
    /**
     * Initialize/Reads parameter file for conversion.
     */

    LOG(info) << "R3BPspxMapped2Precal::ReInit()";
    /*
        fPrecalPar = (R3BPspxPrecalPar*)FairRuntimeDb::instance()->getContainer("R3BPspxPrecalPar");

        if (!fPrecalPar)
        {
            LOG(error) << "Could not get access to R3BPspxPrecalPar-Container.";
            return kFATAL;
        }*/

    SetParContainers();
    SetParameters();
    return kSUCCESS;
}

void R3BPspxMapped2Precal::Exec(Option_t* option)
{
    /**
     * Does the conversion from Mapped to Precal level. It is called for every event.
     * Energies, which are below a channel specific threshold, will be ignored.
     * Applies (strip specific) gains to the energy entries of side 2 of the strip. This is
     * the first calibration step for the position reconstruction.
     */
    for (Int_t d = 0; d < fMappedItems.size(); d++)
    {
        if (!fMappedItems[d])
        {
            printf("Cannot access PSPX%d_%d mapped items\n", (d / 2) + 1, (d % 2) + 1);
            return;
        }
        if (fPrecalPar->GetNumStrips().At(d) == 0)
            continue;

        Int_t nMapped = fMappedItems[d]->GetEntries();

        for (Int_t i = 0; i < nMapped; i++)
        {

            R3BPspxMappedData* mappedData = dynamic_cast<R3BPspxMappedData*>(fMappedItems[d]->At(i));

            // get rid of error message from Febex (GSI firmware)
            if (mappedData->GetEnergy1() == 3075811 || mappedData->GetEnergy1() == 3075810)
                continue;
            if (mappedData->GetEnergy2() == 3075811 || mappedData->GetEnergy2() == 3075810)
                continue;

            if (mappedData->GetStrip1() == mappedData->GetStrip2())
            {
                Int_t strip = mappedData->GetStrip1();
                if (TMath::Abs(mappedData->GetEnergy1()) > threshold1[d][strip - 1] &&
                    TMath::Abs(mappedData->GetEnergy2()) > threshold2[d][strip - 1])
                { // strips
                    Float_t energy1 = mappedData->GetEnergy1();
                    Float_t energy2 = mappedData->GetEnergy2() * gain[d][strip - 1];
                    new ((*fPrecalItems[d])[fPrecalItems[d]->GetEntriesFast()])
                        R3BPspxPrecalData(strip, energy1, energy2);
                }
            }
        }
    }
}

void R3BPspxMapped2Precal::FinishEvent()
{
    for (Int_t i = 0; i < fMappedItems.size(); i++)
    {
        fMappedItems[i]->Clear();
    }
    for (Int_t i = 0; i < fPrecalItems.size(); i++)
    {
        fPrecalItems[i]->Clear();
    }
}

void R3BPspxMapped2Precal::FinishTask() {}

ClassImp(R3BPspxMapped2Precal)
