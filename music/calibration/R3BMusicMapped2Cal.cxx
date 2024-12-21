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

// ---------------------------------------------------------------------
// -----            R3BMusicMapped2Cal source file                 -----
// -----        Created 24/11/19  by J.L. Rodriguez-Sanchez        -----
// ---------------------------------------------------------------------

// ROOT headers
#include "TClonesArray.h"
#include "TMath.h"
#include "TRandom.h"

// Fair headers
#include "FairLogger.h"
#include "FairRootManager.h"
#include "FairRunAna.h"
#include "FairRuntimeDb.h"

// Music headers
#include "R3BMusicCalData.h"
#include "R3BMusicCalPar.h"
#include "R3BMusicMapped2Cal.h"
#include "R3BMusicMappedData.h"

#include <iomanip>

// R3BMusicMapped2Cal: Default Constructor --------------------------
R3BMusicMapped2Cal::R3BMusicMapped2Cal()
    : R3BMusicMapped2Cal("R3BMusic Calibrator", 1)
{
}

// R3BMusicMapped2CalPar: Standard Constructor --------------------------
R3BMusicMapped2Cal::R3BMusicMapped2Cal(const char* name, Int_t iVerbose)
    : FairTask(name, iVerbose)
    , fNumAnodes(MAX_NB_MUSICANODE)   // 8 anodes
    , fNumAnodesRef(MAX_NB_MUSICTREF) // 1 anode for TREF + 1 for trigger
    , fMaxMult(MAX_MULT_MUSIC_CAL)
    , fNumParams(3)
    , fNumPosParams(2)
    , fMaxSigma(200)
    , CalParams(NULL)
    , PosParams(NULL)
    , fCal_Par(NULL)
    , fMusicMappedDataCA(NULL)
    , fMusicCalDataCA(NULL)
    , fOnline(kFALSE)
{
}

// Virtual R3BMusicMapped2Cal: Destructor
R3BMusicMapped2Cal::~R3BMusicMapped2Cal()
{
    LOG(info) << "R3BMusicMapped2Cal: Delete instance";
    if (fMusicMappedDataCA)
        delete fMusicMappedDataCA;
    if (fMusicCalDataCA)
        delete fMusicCalDataCA;
}

void R3BMusicMapped2Cal::SetParContainers()
{
    // Parameter Container
    // Reading musicCalPar from FairRuntimeDb
    FairRuntimeDb* rtdb = FairRuntimeDb::instance();
    if (!rtdb)
    {
        LOG(error) << "FairRuntimeDb not opened!";
    }

    fCal_Par = dynamic_cast<R3BMusicCalPar*>(rtdb->getContainer("musicCalPar"));
    if (!fCal_Par)
    {
        LOG(error) << "R3BMusicMapped2Cal::SetParContainers() Couldn't get handle on musicCalPar container";
    }
    else
    {
        LOG(info) << "R3BMusicMapped2Cal::SetParContainers() musicCalPar container open";
    }
}

void R3BMusicMapped2Cal::SetParameters()
{
    //--- Parameter Container ---
    fNumAnodes = fCal_Par->GetNumAnodes();          // Number of anodes
    fNumParams = fCal_Par->GetNumParamsEFit();      // Number of Parameters
    fNumPosParams = fCal_Par->GetNumParamsPosFit(); // Number of Parameters

    LOG(info) << "R3BMusicMapped2Cal::SetParameters() Nb anodes: " << fNumAnodes;
    LOG(info) << "R3BMusicMapped2Cal::SetParameters() Nb parameters from pedestal fit: " << fNumParams;

    CalParams = new TArrayF();
    Int_t array_size = fNumAnodes * fNumParams;
    CalParams->Set(array_size);
    CalParams = fCal_Par->GetAnodeCalParams(); // Array with the Cal parameters

    LOG(info) << "R3BMusicMapped2Cal::SetParameters() Nb parameters for position fit: " << fNumPosParams;
    PosParams = new TArrayF();
    Int_t array_pos = fNumAnodes * fNumPosParams;
    PosParams->Set(array_pos);
    PosParams = fCal_Par->GetPosParams(); // Array with the Cal parameters

    // Count the number of dead anodes
    Int_t numdeadanodes = 0;
    for (Int_t i = 0; i < fNumAnodes; i++)
        if (CalParams->GetAt(fNumParams * i + 1) == -1)
            numdeadanodes++;
    LOG(info) << "R3BMusicMapped2Cal::SetParameters() Nb of dead anodes in MUSIC : " << numdeadanodes;
}

// -----   Public method Init   --------------------------------------------
InitStatus R3BMusicMapped2Cal::Init()
{
    LOG(info) << "R3BMusicMapped2Cal::Init()";

    // INPUT DATA
    FairRootManager* rootManager = FairRootManager::Instance();
    if (!rootManager)
    {
        LOG(error) << "R3BMusicMapped2Cal::Init() Root-manager not found.";
        return kFATAL;
    }

    fMusicMappedDataCA = dynamic_cast<TClonesArray*>(rootManager->GetObject("MusicMappedData"));
    if (!fMusicMappedDataCA)
    {
        LOG(error) << "R3BMusicMapped2Cal::Init() MusicMappedData not found.";
        return kFATAL;
    }

    // OUTPUT DATA
    // Calibrated data
    fMusicCalDataCA = new TClonesArray("R3BMusicCalData", MAX_MULT_MUSIC_CAL * (fNumAnodes + fNumAnodesRef));
    rootManager->Register("MusicCalData", "MUSIC Cal", fMusicCalDataCA, !fOnline);

    SetParameters();
    return kSUCCESS;
}

// -----   Public method ReInit   ----------------------------------------------
InitStatus R3BMusicMapped2Cal::ReInit()
{
    SetParContainers();
    SetParameters();
    return kSUCCESS;
}

// -----   Public method Execution   --------------------------------------------
void R3BMusicMapped2Cal::Exec(Option_t*)
{
    // Reset entries in output arrays, local arrays
    Reset();

    // Reading the Input -- Mapped Data --
    Int_t nHits = fMusicMappedDataCA->GetEntriesFast();
    // if (nHits != fNumAnodes && nHits > 0)
    //  LOG(warn) << "R3BMusicMapped2Cal: nHits!=" << nHits << " NumAnodes:NumDets" << fNumAnodes << ":" << fNumDets;
    if (nHits == 0)
        return;

    R3BMusicMappedData** mappedData = new R3BMusicMappedData*[nHits];
    UShort_t anodeId = 0;
    Double_t pedestal = 0.;
    Double_t slope = 1.;
    Double_t sigma = 0.;

    for (Int_t i = 0; i < (fNumAnodes + fNumAnodesRef); i++)
    {
        mulanode[i] = 0;
        for (Int_t j = 0; j < fMaxMult; j++)
        {
            energy[j][i] = 0.;
            dtime[j][i] = 0.;
        }
    }

    for (Int_t i = 0; i < nHits; i++)
    {
        mappedData[i] = dynamic_cast<R3BMusicMappedData*>(fMusicMappedDataCA->At(i));
        anodeId = mappedData[i]->GetAnodeID();

        if (anodeId < fNumAnodes && fCal_Par->GetInUse(anodeId + 1) == 1)
        {
            pedestal = CalParams->GetAt(fNumParams * anodeId);
            slope = CalParams->GetAt(fNumParams * anodeId + 1);
            // sigma=CalParams->GetAt(fNumParams*anodeId+2);
            // LOG(info) << detId << " " << anodeId<<" "<< mappedData[i]->GetEnergy()<< " " << pedestal;
            energy[mulanode[anodeId]][anodeId] = pedestal + slope * mappedData[i]->GetEnergy();
            dtime[mulanode[anodeId]][anodeId] = mappedData[i]->GetTime();
            mulanode[anodeId]++;
        }
        else if (anodeId >= fNumAnodes)
        {
            // LOG(info) <<"a="<< anodeId<<" e="<< mappedData[i]->GetEnergy()<< "  t=" << mappedData[i]->GetTime();
            dtime[mulanode[anodeId]][anodeId] = mappedData[i]->GetTime(); // Ref. Time
            mulanode[anodeId]++;
        }
    }

    // Fill data only if there are trigger and TREF signals
    if (mulanode[fNumAnodes] == 1 && mulanode[fNumAnodes + 1] == 1)
    {
        for (Int_t i = 0; i < fNumAnodes; i++)
        {
            Float_t a0 = PosParams->GetAt(fNumPosParams * i);
            Float_t a1 = PosParams->GetAt(fNumPosParams * i + 1);
            for (Int_t j = 0; j < mulanode[fNumAnodes]; j++)
                for (Int_t k = 0; k < mulanode[i]; k++)
                {
                    if (energy[k][i] > 0.)
                        AddCalData(i, a0 + a1 * (dtime[k][i] - dtime[j][fNumAnodes]), energy[k][i]);
                }
        }
    }
    if (mappedData)
        delete[] mappedData;
    return;
}

// -----   Public method Reset   ------------------------------------------------
void R3BMusicMapped2Cal::Reset()
{
    LOG(debug) << "Clearing MusicCalData Structure";
    if (fMusicCalDataCA)
        fMusicCalDataCA->Clear();
}

// -----   Private method AddCalData  --------------------------------------------
R3BMusicCalData* R3BMusicMapped2Cal::AddCalData(UShort_t aid, Double_t dt, Double_t e)
{
    // It fills the R3BMusicCalData
    TClonesArray& clref = *fMusicCalDataCA;
    Int_t size = clref.GetEntriesFast();
    return new (clref[size]) R3BMusicCalData(aid, dt, e);
}

ClassImp(R3BMusicMapped2Cal);
