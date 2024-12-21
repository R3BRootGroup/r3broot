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
#include "TDecompSVD.h"
#include "TMath.h"
#include "TRandom.h"
#include "TVectorD.h"

// Fair headers
#include "FairLogger.h"
#include "FairRootManager.h"
#include "FairRunAna.h"
#include "FairRuntimeDb.h"

// Musli headers
#include "R3BLogger.h"
#include "R3BMusliCal2Hit.h"
#include "R3BMusliHitPar.h"
#include "R3BTGeoPar.h"

#include <iomanip>

// R3BMusliCal2Hit: Default Constructor --------------------------
R3BMusliCal2Hit::R3BMusliCal2Hit()
    : R3BMusliCal2Hit("R3BMusli Calibrator", 1)
{
}

// R3BMusliCal2HitPar: Standard Constructor --------------------------
R3BMusliCal2Hit::R3BMusliCal2Hit(const char* name, Int_t iVerbose)
    : FairTask(name, iVerbose)
    , fNumGroups(MAX_NUM_GROUPS_CAL)
    , fNumTypes(MAX_NUM_TYPES_HIT)
    , fMaxMult(MAX_MULT_MUSLI)
    , fNumEcorrBeta(2)
    , fNumZ(3)
    , fEaveVsBetaHitParams(NULL)
    , fEcorrBetaHitParams(NULL)
    , fZHitParams(NULL)
    , fHit_Par(NULL)
    , fMusliCalDataCA(NULL)
    , fMusliHitDataCA(NULL)
    , fFrsDataCA(NULL)
    , fFrsSciTofCal(kFALSE)
    , fOnline(kFALSE)
    , fExpId(522) // s522 at 1.2 A.GeV, s509 at 400 A.MeV
    , fDirectBeta(0.876)
    , fIdS2(1)
    , fIdCaveC(2)
{
}

// Virtual R3BMusliCal2Hit: Destructor
R3BMusliCal2Hit::~R3BMusliCal2Hit()
{
    R3BLOG(debug1, "R3BMusliCal2Hit: Delete instance");
    if (fMusliCalDataCA)
        delete fMusliCalDataCA;
    if (fMusliHitDataCA)
        delete fMusliHitDataCA;
    if (fFrsDataCA)
        delete fFrsDataCA;
}

void R3BMusliCal2Hit::SetParContainers()
{
    // Parameter Container
    // Reading musicHitPar from FairRuntimeDb
    FairRuntimeDb* rtdb = FairRuntimeDb::instance();
    R3BLOG_IF(fatal, !rtdb, "R3BMusliCal2Hit::SetParContainers(), FairRuntimeDb not found");

    fHit_Par = dynamic_cast<R3BMusliHitPar*>(rtdb->getContainer("musliHitPar"));
    if (!fHit_Par)
    {
        R3BLOG(error, "R3BMusliCal2Hit::SetParContainers() Couldn't get handle on musliHitPar container");
    }
    else
    {
        R3BLOG(info, "R3BMusliCal2Hit::SetParContainers() musliHitPar container open");
    }

    // Reading the TGeoPar from the FairRun
    fMusliGeo_Par = dynamic_cast<R3BTGeoPar*>(rtdb->getContainer("MusliGeoPar"));
    R3BLOG_IF(error, !fMusliGeo_Par, "R3BMusliCal2Hit::SetParContainers() Couldn´t access to MusliGeoPar container.");
}

void R3BMusliCal2Hit::SetParameters()
{
    //--- Parameter Container ---
    fNumGroups = fHit_Par->GetNumGroupsAnodes();
    fNumTypes = fHit_Par->GetNumTypes();
    fMaxMult = fHit_Par->GetMaxMult();
    fNumEcorrBeta = fHit_Par->GetNumParamsEcorrBetaFit();
    fNumZ = fHit_Par->GetNumParamsZFit();

    fEaveVsBetaHitParams = new TArrayD();
    fEaveVsBetaHitParams->Set(fNumTypes);
    fEaveVsBetaHitParams = fHit_Par->GetEaveVsBetaHitParams();

    fEcorrBetaHitParams = new TArrayD();
    fEcorrBetaHitParams->Set(fNumTypes);
    fEcorrBetaHitParams = fHit_Par->GetEcorrBetaHitParams();

    fZHitParams = new TArrayD();
    fZHitParams->Set(fNumTypes);
    fZHitParams = fHit_Par->GetZHitParams();
}

// -----   Public method Init   --------------------------------------------
InitStatus R3BMusliCal2Hit::Init()
{
    R3BLOG(info, "R3BMusliCal2Hit::Init()");

    // INPUT DATA
    FairRootManager* rootManager = FairRootManager::Instance();
    if (!rootManager)
    {
        R3BLOG(fatal, "R3BMusliCal2Hit::Init() FairRootManager not found.");
        return kFATAL;
    }

    fMusliCalDataCA = dynamic_cast<TClonesArray*>(rootManager->GetObject("MusliCalData"));
    if (!fMusliCalDataCA)
    {
        R3BLOG(info, "R3BMusliCal2Hit::Init() MusliCalData not found.");
    }

    fFrsDataCA = dynamic_cast<TClonesArray*>(rootManager->GetObject("FrsSciTofCalData"));
    if (fFrsDataCA)
    {
        fFrsSciTofCal = kTRUE;
    }
    else
    {
        fFrsSciTofCal = kFALSE;
        fFrsDataCA = dynamic_cast<TClonesArray*>(rootManager->GetObject("FrsData"));
        if (!fFrsDataCA)
        {
            R3BLOG(info, "R3BMusliCal2Hit::Init() FrsSciTofCalData not found, try FrsData");
            R3BLOG(info, "R3BMusliCal2Hit::Init() FrsData not found neither");
        }
    }

    // OUTPUT DATA
    fMusliHitDataCA = new TClonesArray("R3BMusliHitData");
    rootManager->Register("MusliHitData", "MUSLI Hit", fMusliHitDataCA, !fOnline);
    fMusliHitDataCA->Clear();

    SetParameters();
    return kSUCCESS;
}

// -----   Public method ReInit   ----------------------------------------------
InitStatus R3BMusliCal2Hit::ReInit()
{
    SetParContainers();
    SetParameters();
    return kSUCCESS;
}

// -----   Public method Execution   --------------------------------------------
void R3BMusliCal2Hit::Exec(Option_t* option)
{
    // Reset entries in output arrays, local arrays
    Reset();

    // Get the Beta of the incoming beam
    Double_t beta = fDirectBeta;
    if (fFrsSciTofCal == kTRUE && fFrsDataCA && fFrsDataCA->GetEntriesFast() > 0)
    {
        UInt_t nHits = fFrsDataCA->GetEntriesFast();
        for (UInt_t ihit = 0; ihit < nHits; ihit++)
        {
            R3BFrsSciTofCalData* hittofcal = dynamic_cast<R3BFrsSciTofCalData*>(fFrsDataCA->At(ihit));
            if (!hittofcal)
                continue;
            if (hittofcal->GetDetIdSta() == fIdS2 && hittofcal->GetDetIdSto() == fIdCaveC)
            {
                beta = hittofcal->GetBeta();
            }
        }
    }
    else if (fFrsDataCA && fFrsDataCA->GetEntriesFast() == 1)
    {
        R3BFrsData* hit = dynamic_cast<R3BFrsData*>(fFrsDataCA->At(0));
        if (hit)
        {
            beta = hit->GetBeta();
        }
        else
        {
            beta = 0.;
        }
    }

    // Reading the Input -- Cal Data --
    for (Int_t i = 0; i < fNumGroups; i++)
    {
        mult_cal[i] = 0;
        for (Int_t j = 0; j < fMaxMult; j++)
        {
            dt_cal[j][i] = -1000;
            e_cal[j][i] = 0;
        }
    }

    Int_t signal_cal; // 0-based [0..14]
    Int_t nHits = fMusliCalDataCA->GetEntriesFast();
    if (nHits == 0)
        return;
    R3BMusliCalData** calData = new R3BMusliCalData*[nHits];
    for (Int_t i = 0; i < nHits; i++)
    {
        calData[i] = dynamic_cast<R3BMusliCalData*>(fMusliCalDataCA->At(i));
        signal_cal = calData[i]->GetSignal() - 1;
        dt_cal[mult_cal[signal_cal]][signal_cal] = calData[i]->GetDT();
        e_cal[mult_cal[signal_cal]][signal_cal] = calData[i]->GetEnergy();
        mult_cal[signal_cal]++;
    }

    // Calculate and fill the Output -- Hit Data --
    for (Int_t i = 0; i < fNumTypes; i++)
    {
        nba[i] = 0;
        Esum[i] = 0;
        eave[i] = 0;
        e_hit[i] = 0;
        z_hit[i] = 0;
        x_hit[i] = -1000;
        theta_hit[i] = -1000;
    }

    // Calculate the energy sum for fNbAnodes = 2: fSignal = [1..8]
    for (Int_t i = 0; i < 8; i++)
    {
        if (mult_cal[i] == 1 && e_cal[0][i] > 0)
        {
            Esum[0] += e_cal[0][i];
            nba[0]++;
        }
        if (mult_cal[3] == 1 && mult_cal[4] == 1)
            x_hit[0] = 0.5 * (dt_cal[0][3] + dt_cal[0][4]);
    }

    // Calculate the energy sum for fNbAnodes = 4: fSignal = [9..12]
    for (Int_t i = 8; i < 12; i++)
    {
        if (mult_cal[i] == 1 && e_cal[0][i] > 0)
        {
            Esum[1] += e_cal[0][i];
            nba[1]++;
        }
        if (mult_cal[9] == 1 && mult_cal[10] == 1)
            x_hit[1] = 0.5 * (dt_cal[0][9] + dt_cal[0][10]);
    }

    // Calculate the energy sum for fNbAnodes = 8: fSignal = [13,14]
    for (Int_t i = 12; i < 14; i++)
    {
        if (mult_cal[i] == 1 && e_cal[0][i] > 0)
        {
            Esum[2] += e_cal[0][i];
            nba[2]++;
        }
        if (mult_cal[12] == 1 && mult_cal[13] == 1)
            x_hit[2] = 0.5 * (dt_cal[0][12] + dt_cal[0][13]);
    }

    // Calculate the energy sum for fNbAnodes = 16: fSignal = [15]
    if (mult_cal[14] == 1 && e_cal[0][14] > 0)
    {
        Esum[3] = e_cal[0][14];
        nba[3] = 1;
        x_hit[3] = dt_cal[0][14];
    }

    for (Int_t type = 1; type <= fNumTypes; type++)
    {
        // Calculate theta value only if type = 1, exclude first and last pairs of anodes
        if (type == 1 && nba[type - 1] == 8)
        { // one ion only along the full Z path in the music

            Double_t musliZpos[6];
            Double_t musliXpos[6];
            for (int i = 0; i < 6; i++)
            {
                // GetPosZ() gives the center of the chamber in cm
                // width of one anode is 25. mm -> one pair is 50.
                // starts the theta calculation in the middle of the second pair
                //                              i.e.  125 mm. upstream GetPosZ()
                musliZpos[i] = fMusliGeo_Par->GetPosZ() * 10. - 125. + i * 50.; // [mm]
                musliXpos[i] = dt_cal[0][i + 1];                                // only mult_cal[i] == 1 is selected
            }
            TVectorD VectX;
            VectX.Use(6, musliXpos);
            TVectorD VectZ;
            VectZ.Use(6, musliZpos);
            TMatrixD MatXZ(6, 2);
            TMatrixDColumn(MatXZ, 0) = 1.0;
            TMatrixDColumn(MatXZ, 1) = VectZ;
            TDecompSVD svdFit(MatXZ);
            Bool_t bFit;
            const TVectorD VectSvdFit = svdFit.Solve(VectX, bFit);
            theta_hit[type - 1] = TMath::ATan(VectSvdFit[1]);
        }

        if (nba[type - 1] > 0)
        {
            eave[type - 1] = Esum[type - 1] / nba[type - 1];
            e_hit[type - 1] = eave[type - 1] * fEaveVsBetaHitParams->GetAt(type - 1) /
                              BetaCorr_std(beta,
                                           fEcorrBetaHitParams->GetAt((type - 1) * fNumEcorrBeta),
                                           fEcorrBetaHitParams->GetAt((type - 1) * fNumEcorrBeta + 1));
        }
        if (e_hit[type - 1] > 0)
        {
            z_hit[type - 1] = SqrtE2Z(sqrt(e_hit[type - 1]),
                                      fZHitParams->GetAt((type - 1) * fNumZ),
                                      fZHitParams->GetAt((type - 1) * fNumZ + 1),
                                      fZHitParams->GetAt((type - 1) * fNumZ + 2));
            AddHitData(type, e_hit[type - 1], z_hit[type - 1], x_hit[type - 1], theta_hit[type - 1]);
        }
    }

    if (calData)
        delete[] calData;

    return;
}

// -----   Public method Reset   ------------------------------------------------
void R3BMusliCal2Hit::Reset()
{
    R3BLOG(debug1, "Clearing MusliHitData Structure");
    if (fMusliHitDataCA)
        fMusliHitDataCA->Clear();
}

// -----   Private method AddCalData  --------------------------------------------
R3BMusliHitData* R3BMusliCal2Hit::AddHitData(UInt_t type, Double_t e, Double_t z, Double_t x, Double_t theta)
{
    // It fills the R3BMusliCalData
    TClonesArray& clref = *fMusliHitDataCA;
    Int_t size = clref.GetEntriesFast();
    return new (clref[size]) R3BMusliHitData(type, e, z, x, theta);
}

ClassImp(R3BMusliCal2Hit);
