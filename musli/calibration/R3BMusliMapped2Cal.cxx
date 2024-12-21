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

// ROOT headers
#include "TClonesArray.h"
#include "TMath.h"
#include "TRandom.h"

// Fair headers
#include "FairLogger.h"
#include "FairRootManager.h"
#include "FairRunAna.h"
#include "FairRuntimeDb.h"

#include "R3BLogger.h"

// Musli headers
#include "R3BEventHeader.h"
#include "R3BMusliCalData.h"
#include "R3BMusliCalPar.h"
#include "R3BMusliMapped2Cal.h"
#include "R3BMusliMappedData.h"

#include <iomanip>

// R3BMusliMapped2Cal: Default Constructor --------------------------
R3BMusliMapped2Cal::R3BMusliMapped2Cal()
    : R3BMusliMapped2Cal("R3BMusli Calibrator", 1)
{
}

// R3BMusliMapped2CalPar: Standard Constructor --------------------------
R3BMusliMapped2Cal::R3BMusliMapped2Cal(const char* name, Int_t iVerbose)
    : FairTask(name, iVerbose)
    , fNumSignals(18)
    , fNumGroupsAnodes(15)
    , fNumParamsEneFit(2)
    , fNumParamsPosFit(2)
    , fNumParamsMultHit(2)
    , fMaxMult(20)
    , fEneCalParams(NULL)
    , fPosCalParams(NULL)
    , fMultHitCalParams(NULL)
    , fCal_Par(NULL)
    , fMusliMappedDataCA(NULL)
    , fMusliCalDataCA(NULL)
    , fOnline(kFALSE)
    , fHeader(NULL)
    , winL(0.)
    , winR(0.)
    , fUseMultHit(kFALSE)
{
}

// Virtual R3BMusliMapped2Cal: Destructor
R3BMusliMapped2Cal::~R3BMusliMapped2Cal()
{
    R3BLOG(debug1, "R3BMusliMapped2Cal: Delete instance");
    if (fMusliMappedDataCA)
        delete fMusliMappedDataCA;
    if (fMusliCalDataCA)
        delete fMusliCalDataCA;
}

void R3BMusliMapped2Cal::SetParContainers()
{
    // Parameter Container
    // Reading musicCalPar from FairRuntimeDb
    FairRuntimeDb* rtdb = FairRuntimeDb::instance();
    R3BLOG_IF(fatal, !rtdb, "FairRuntimeDb not found");

    fCal_Par = dynamic_cast<R3BMusliCalPar*>(rtdb->getContainer("musliCalPar"));
    if (!fCal_Par)
    {
        R3BLOG(error, "R3BMusliMapped2Cal::SetParContainers() Couldn't get handle on musliCalPar container");
    }
    else
    {
        R3BLOG(info, "R3BMusliMapped2Cal::SetParContainers() musliCalPar container open");
    }
}

void R3BMusliMapped2Cal::SetParameters()
{
    //--- Parameter Container ---
    fNumSignals = fCal_Par->GetNumSignals();           // Number of signals at Mapped
    fNumGroupsAnodes = fCal_Par->GetNumGroupsAnodes(); // Number of groups of anodes
    fNumParamsEneFit = fCal_Par->GetNumParamsEneFit(); // Number of ene calib params per group of anodes
    fNumParamsPosFit = fCal_Par->GetNumParamsPosFit(); // Number of pos calib params per group of anodes
    fMaxMult = fCal_Par->GetMaxMult();                 // maximum number of hit per output channel

    Int_t array_ene = fNumGroupsAnodes * fNumParamsEneFit; // array of ene calibration parameters
    fEneCalParams = new TArrayF();
    fEneCalParams->Set(array_ene);
    fEneCalParams = fCal_Par->GetEneCalParams();

    Int_t array_pos = fNumGroupsAnodes * fNumParamsPosFit; // arrays of pos calibration parameters
    fPosCalParams = new TArrayF();
    fPosCalParams->Set(array_pos);
    fPosCalParams = fCal_Par->GetPosCalParams();

    Int_t array_mult = fNumGroupsAnodes * fNumParamsMultHit; // array of dtraw limits for map_multi-hit analysis
    fMultHitCalParams = new TArrayF();
    fMultHitCalParams->Set(array_mult);
    fMultHitCalParams = fCal_Par->GetMultHitCalParams();
}

// -----   Public method Init   --------------------------------------------
InitStatus R3BMusliMapped2Cal::Init()
{
    R3BLOG(info, "R3BMusliMapped2Cal::Init()");

    // INPUT DATA
    FairRootManager* rootManager = FairRootManager::Instance();
    if (!rootManager)
    {
        R3BLOG(fatal, "R3BMusliMapped2Cal::Init() FairRootManager not found.");
        return kFATAL;
    }

    fHeader = dynamic_cast<R3BEventHeader*>(rootManager->GetObject("EventHeader."));

    fMusliMappedDataCA = dynamic_cast<TClonesArray*>(rootManager->GetObject("MusliMappedData"));
    if (!fMusliMappedDataCA)
    {
        R3BLOG(fatal, "R3BMusliMapped2Cal::Init() MusliMappedData not found.");
        return kFATAL;
    }

    // OUTPUT DATA
    fMusliCalDataCA = new TClonesArray("R3BMusliCalData");
    rootManager->Register("MusliCalData", "MUSLI Cal", fMusliCalDataCA, !fOnline);
    fMusliCalDataCA->Clear();

    SetParameters();
    return kSUCCESS;
}

// -----   Public method ReInit   ----------------------------------------------
InitStatus R3BMusliMapped2Cal::ReInit()
{
    SetParContainers();
    SetParameters();
    return kSUCCESS;
}

// -----   Public method Execution   --------------------------------------------
void R3BMusliMapped2Cal::Exec(Option_t* option)
{
    // Reset entries in output arrays, local arrays
    Reset();

    // variables to read map data
    UInt_t sig = 0;
    UInt_t map_mult[fNumSignals];
    UInt_t map_sig[fNumSignals * fMaxMult];
    Double_t map_eraw[fNumSignals * fMaxMult];
    Double_t map_traw[fNumSignals * fMaxMult];
    Double_t pospar[fNumParamsPosFit];
    Double_t enepar[fNumParamsEneFit];
    Double_t dtraw, dtcal, ecal;

    // Reading the Input -- Mapped Data --
    UInt_t nHits = fMusliMappedDataCA->GetEntries();
    if (nHits == 0)
        return;

    // init tmp arrays
    for (Int_t i = 0; i < fNumSignals; i++)
    {
        map_mult[i] = 0;
        for (Int_t j = 0; j < fMaxMult; j++)
        {
            map_sig[i * fMaxMult + j] = 0;
            map_eraw[i * fMaxMult + j] = 0.;
            map_traw[i * fMaxMult + j] = 0.;
        }
    }

    // loop over the entries of the fMapped TClonesArray
    for (UInt_t i = 0; i < nHits; i++)
    {
        R3BMusliMappedData* map = (R3BMusliMappedData*)fMusliMappedDataCA->At(i);

        sig = map->GetSignal() - 1; // 0-based signal Id

        if (0 <= sig && sig < fNumSignals)
        {
            if (map_mult[sig] < fMaxMult)
            {
                map_sig[sig * fMaxMult + map_mult[sig]] = sig + 1;
                map_eraw[sig * fMaxMult + map_mult[sig]] = (Double_t)map->GetEnergy();
                map_traw[sig * fMaxMult + map_mult[sig]] = (Double_t)map->GetTime();
            }
            map_mult[sig]++;
        }
        else
            LOG(info) << "R3BMusliMapped2Cal::Exec WARNING, signal id is " << map->GetSignal();
    }

    if (fNumGroupsAnodes <= fNumSignals)
    {
        for (Int_t i = 0; i < fNumGroupsAnodes; i++)
        {
            if (map_mult[i] == 0) // no data
                continue;
            // else if (map_mult[16] == 1 && map_mult[i] == 1 && map_traw[i * fMaxMult] > 0 &&
            // hot fix for s091 and s118 due to a bad cable first Tref hit is the good one
            else if (0 < map_mult[16] && map_mult[16] <= 2 && map_mult[i] == 1 && map_traw[i * fMaxMult] > 0 &&
                     map_eraw[i * fMaxMult] > 0) // one Tref and one [i] signal with real values
            {
                dtcal = 0.;
                for (Int_t k = 0; k < fNumParamsPosFit; k++)
                {
                    dtraw = map_traw[i * fMaxMult] - map_traw[16 * fMaxMult];
                    pospar[k] = fPosCalParams->GetAt(fNumParamsPosFit * i + k);
                    dtcal += pospar[k] * pow(dtraw, k);
                }
                ecal = 0;
                for (Int_t k = 0; k < fNumParamsEneFit; k++)
                {
                    enepar[k] = fEneCalParams->GetAt(fNumParamsEneFit * i + k);
                    ecal += enepar[k] * pow(map_eraw[i * fMaxMult], k);
                }
                if (ecal > 0.)
                {
                    AddCalData(map_sig[i * fMaxMult], dtcal, ecal);
                }
            }
            else if (map_mult[16] == 2 && map_mult[i] == 1 && map_traw[i * fMaxMult] > 0 &&
                     map_eraw[i * fMaxMult] > 0) // first Tref and one [i] signal with real values
            {
                dtcal = 0.;
                for (Int_t k = 0; k < fNumParamsPosFit; k++)
                {
                    dtraw = map_traw[i * fMaxMult] - map_traw[16 * fMaxMult];
                    pospar[k] = fPosCalParams->GetAt(fNumParamsPosFit * i + k);
                    dtcal += pospar[k] * pow(dtraw, k);
                }
                ecal = 0;
                for (Int_t k = 0; k < fNumParamsEneFit; k++)
                {
                    enepar[k] = fEneCalParams->GetAt(fNumParamsEneFit * i + k);
                    ecal += enepar[k] * pow(map_eraw[i * fMaxMult], k);
                }
                if (ecal > 0.)
                    AddCalData(map_sig[i * fMaxMult], dtcal, ecal);
            }
            else if ((map_mult[16] > 1 || map_mult[i] > 1) && fUseMultHit)
            {
                LOG(info) << "enter in multi hit case";
                Double_t good_tref = -1.;
                Int_t no_of_tref = 0;
                if (map_mult[16] == 1)
                {
                    no_of_tref = 1;
                    good_tref = map_traw[16 * fMaxMult];
                }
                else if (map_mult[16] > 1)
                {
                    good_tref = map_traw[16 * fMaxMult];
                    for (Int_t j = 0; j < map_mult[16]; j++)
                    {
                        // TO DO : need to add some comments to explain this condition
                        Double_t cfd_t = ((map_traw[16 * fMaxMult + j] * 25. / 256.) - fHeader->GetTStart());
                        if (cfd_t > winL && cfd_t < winR)
                        {
                            no_of_tref++;
                            good_tref = map_traw[16 * fMaxMult + j];
                        }
                    }
                }
                if (no_of_tref == 1 && good_tref > 0 && map_mult[i] < fMaxMult)
                {
                    for (Int_t j = 0; j < map_mult[i]; j++)
                    {
                        if (map_traw[i * fMaxMult + j] > 0 && map_eraw[i * fMaxMult + j] > 0)
                        {
                            dtraw = map_traw[i * fMaxMult + j] - good_tref;
                            if (dtraw > fMultHitCalParams->GetAt(fNumParamsMultHit * i + 0) &&
                                dtraw < fMultHitCalParams->GetAt(fNumParamsMultHit * i + 1))
                            {
                                dtcal = 0.;
                                for (Int_t k = 0; k < fNumParamsPosFit; k++)
                                {
                                    pospar[k] = fPosCalParams->GetAt(fNumParamsPosFit * i + k);
                                    dtcal += pospar[k] * pow(dtraw, k);
                                }
                                ecal = 0;
                                for (Int_t k = 0; k < fNumParamsEneFit; k++)
                                {
                                    enepar[k] = fEneCalParams->GetAt(fNumParamsEneFit * i + k);
                                    ecal += enepar[k] * pow(map_eraw[i * fMaxMult + j], k);
                                }
                                if (ecal > 0.)
                                    AddCalData(map_sig[i * fMaxMult + j], dtcal, ecal);
                            }
                        }
                    }
                }
            }
        }
    }
    return;
}

// -----   Public method Reset   ------------------------------------------------
void R3BMusliMapped2Cal::Reset()
{
    R3BLOG(debug1, "Clearing MusliCalData Structure");
    if (fMusliCalDataCA)
        fMusliCalDataCA->Clear();
}

// -----   Private method AddCalData  --------------------------------------------
R3BMusliCalData* R3BMusliMapped2Cal::AddCalData(UInt_t sig, Double_t dt, Double_t e)
{
    // It fills the R3BMusliCalData
    TClonesArray& clref = *fMusliCalDataCA;
    Int_t size = clref.GetEntriesFast();
    return new (clref[size]) R3BMusliCalData(sig, dt, e);
}

ClassImp(R3BMusliMapped2Cal);
