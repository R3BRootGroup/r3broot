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

// -------------------------------------------------------------
// -----         R3BMusicCal2Hit source file               -----
// -----    Created 30/11/19  by J.L. Rodriguez-Sanchez    -----
// -------------------------------------------------------------

// ROOT headers
#include "TClonesArray.h"
#include "TDecompSVD.h"
#include "TMath.h"
#include "TMatrixD.h"
#include "TRandom.h"
#include "TVectorD.h"

// Fair headers
#include "FairLogger.h"
#include "FairRootManager.h"
#include "FairRunAna.h"
#include "FairRuntimeDb.h"

#include <iomanip>

// Music headers
#include "R3BMusicCal2Hit.h"
#include "R3BMusicCalData.h"
#include "R3BMusicHitData.h"
#include "R3BMusicHitPar.h"

// R3BMusicCal2Hit: Default Constructor --------------------------
R3BMusicCal2Hit::R3BMusicCal2Hit()
    : R3BMusicCal2Hit("R3BMusicCal2Hit", 1)
{
}

// R3BMusicCal2HitPar: Standard Constructor --------------------------
R3BMusicCal2Hit::R3BMusicCal2Hit(const char* name, Int_t iVerbose)
    : FairTask(name, iVerbose)
    , fNumAnodes(8)
    , fNumAnodesAngleFit(0)
    , fNumParams(2)
    , CalZParams(NULL)
    , fCal_Par(NULL)
    , fMusicHitDataCA(NULL)
    , fMusicCalDataCA(NULL)
    , fZ0(0.)
    , fZ1(0.)
    , fZ2(0.)
    , fOnline(kFALSE)
    , fx0_point(0.0)
    , fy0_point(0.0)
    , frot_ang(0.0)
    , fSim(false)
{
}

// Virtual R3BMusicCal2Hit: Destructor
R3BMusicCal2Hit::~R3BMusicCal2Hit()
{
    LOG(debug) << "R3BMusicCal2Hit::Destructor";
    if (fMusicHitDataCA)
        delete fMusicHitDataCA;
}

void R3BMusicCal2Hit::SetParContainers()
{
    // Parameter Container
    // Reading MusicCalPar from FairRuntimeDb
    FairRuntimeDb* rtdb = FairRuntimeDb::instance();
    if (!rtdb)
    {
        LOG(error) << "FairRuntimeDb not found";
    }

    fCal_Par = dynamic_cast<R3BMusicHitPar*>(rtdb->getContainer("musicHitPar"));
    if (!fCal_Par)
    {
        LOG(error) << "R3BMusicCal2HitPar::Init() Couldn't get handle on musicHitPar container";
    }
    else
    {
        LOG(info) << "R3BMusicCal2HitPar:: musicHitPar container open";
    }
}

void R3BMusicCal2Hit::SetParameter()
{
    //--- Parameter Container ---
    fNumAnodes = fCal_Par->GetNumAnodes(); // Number of anodes

    // Anodes that don't work set to zero
    for (Int_t i = 0; i < fNumAnodes; i++)
    {
        fStatusAnodes[i] = fCal_Par->GetInUse(i + 1);
    }

    fNumParams = fCal_Par->GetNumParZFit(); // Number of Parameters
    LOG(info) << "R3BMusicCal2Hit::SetParameter() Nb parameters for charge-Z: " << fNumParams;
    CalZParams = new TArrayF();
    CalZParams->Set(fNumParams);
    CalZParams = fCal_Par->GetZHitPar(); // Array with the Cal parameters

    CalAngCorParams = new TArrayF();
    CalAngCorParams->Set(3);
    CalAngCorParams = fCal_Par->GetAngCorPar(); // Array with the Cal parameters
    fx0_point = CalAngCorParams->GetAt(0);
    fy0_point = CalAngCorParams->GetAt(1);
    frot_ang = CalAngCorParams->GetAt(2);
    // Parameters detector
    if (fNumParams == 2)
    {
        LOG(info) << "R3BMusicCal2Hit parameters for charge-Z: " << CalZParams->GetAt(0) << " : "
                  << CalZParams->GetAt(1);
        fZ0 = CalZParams->GetAt(0);
        fZ1 = CalZParams->GetAt(1);
    }
    else if (fNumParams == 3)
    {
        LOG(info) << "R3BMusicCal2Hit parameters for charge-Z:" << CalZParams->GetAt(0) << " : " << CalZParams->GetAt(1)
                  << " : " << CalZParams->GetAt(2);
        fZ0 = CalZParams->GetAt(0);
        fZ1 = CalZParams->GetAt(1);
        fZ2 = CalZParams->GetAt(2);
    }
    else
        LOG(error) << "R3BMusicCal2Hit parameters for charge-Z cannot be used here, number of parameters must be < 4, "
                      "currently it is "
                   << fNumParams;
}

// -----   Public method Init   --------------------------------------------
InitStatus R3BMusicCal2Hit::Init()
{
    LOG(info) << "R3BMusicCal2Hit::Init()";

    // INPUT DATA
    FairRootManager* rootManager = FairRootManager::Instance();
    if (!rootManager)
    {
        LOG(error) << "R3BMusicCal2Hit::Init() Root-manager not found.";
        return kFATAL;
    }

    fMusicCalDataCA = dynamic_cast<TClonesArray*>(rootManager->GetObject("MusicCalData"));
    if (!fMusicCalDataCA)
    {
        LOG(error) << "R3BMusicCal2Hit::Init() MusicCalData not found.";
        return kFATAL;
    }

    // OUTPUT DATA
    fMusicHitDataCA = new TClonesArray("R3BMusicHitData");
    rootManager->Register("MusicHitData", "Music Hit", fMusicHitDataCA, !fOnline);

    SetParameter();
    return kSUCCESS;
}

// -----   Public method ReInit   ----------------------------------------------
InitStatus R3BMusicCal2Hit::ReInit()
{
    SetParContainers();
    SetParameter();
    return kSUCCESS;
}

// -----   Public method Execution   --------------------------------------------
void R3BMusicCal2Hit::Exec(Option_t*)
{
    // Reset entries in output arrays, local arrays
    Reset();

    Int_t nHits = fMusicCalDataCA->GetEntriesFast();
    if (nHits == 0)
    {
        return;
    }

    // This option is only for simulations
    if (fSim)
    {
        ExecSim(nHits);
        return;
    }

    // R3BMusicCalData* CalDat;
    R3BMusicCalData** CalDat = new R3BMusicCalData*[nHits];

    Int_t secId, anodeId;
    Double_t energyperanode[fNumAnodes], dt[fNumAnodes], good_dt[fNumAnodes];
    Int_t nbdet = 0;

    for (Int_t j = 0; j < fNumAnodes; j++)
    {
        energyperanode[j] = 0.;
        dt[j] = 0.;
    }

    for (Int_t j = 0; j < fNumAnodes; j++)
    {
        good_dt[j] = 0.;
    }

    for (Int_t i = 0; i < nHits; i++)
    {
        CalDat[i] = dynamic_cast<R3BMusicCalData*>(fMusicCalDataCA->At(i));
        anodeId = CalDat[i]->GetAnodeID();
        energyperanode[anodeId] = CalDat[i]->GetEnergy();
        dt[anodeId] = CalDat[i]->GetDTime();
    }

    Double_t nba = 0., theta = 0., Esum = 0.;
    // calculate truncated dE from 8 anodes, MUSIC
    fNumAnodesAngleFit = 0;
    for (Int_t j = 0; j < fNumAnodes; j++)
    {
        if (energyperanode[j] > 0. && fStatusAnodes[j] == 1)
        {
            Esum = Esum + energyperanode[j];
            good_dt[fNumAnodesAngleFit] = dt[j];
            fPosAnodes[fNumAnodesAngleFit] = fCal_Par->GetAnodePos(j + 1);
            fNumAnodesAngleFit++;
            nba++;
        }
    }

    // if(nba!=fNumAnodesAngleFit)LOG(error) << "R3BMusicCal2Hit::nba("<< nba<<") and
    // fNumAnodesAngleFit("<<fNumAnodesAngleFit <<") are different";

    if (fNumAnodesAngleFit > 2 && Esum / nba > 0.)
    {
        fPosZ.Use(fNumAnodesAngleFit, fPosAnodes);
        TMatrixD A(fNumAnodesAngleFit, 2);
        TMatrixDColumn(A, 0) = 1.0;
        TMatrixDColumn(A, 1) = fPosZ;

        // singular value decomposition of matrix A
        //(generalizes the eigendecomposition of a square normal matrix with an orthonormal eigenbasis
        // to any m × n matrix)
        TDecompSVD svd(A);
        Bool_t ok;
        TVectorD dt_r;
        dt_r.Use(fNumAnodesAngleFit, good_dt);
        TVectorD c_svd_r = svd.Solve(dt_r, ok);
        theta = c_svd_r[1];

        Double_t sqrtEsum_nba_rot = fy0_point + (theta * 1000 - fx0_point) * sin(frot_ang) +
                                    (TMath::Sqrt(Esum / nba) - fy0_point) * cos(frot_ang);
        Double_t Esum_nba_rot = sqrtEsum_nba_rot * sqrtEsum_nba_rot;
        Double_t zhit = fZ0 + fZ1 * sqrtEsum_nba_rot + fZ2 * Esum_nba_rot;
        if (zhit > 0)
            AddHitData(theta, zhit, Esum_nba_rot);
    }
    if (CalDat)
        delete[] CalDat;
    return;
}

void R3BMusicCal2Hit::ExecSim(int nHits)
{

    R3BMusicCalData** CalDat = new R3BMusicCalData*[nHits];

    Int_t secId, anodeId;
    Double_t energyperanode[fNumAnodes], dt[fNumAnodes], good_dt[fNumAnodes];
    Int_t nbdet = 0;

    for (Int_t j = 0; j < fNumAnodes; j++)
    {
        energyperanode[j] = 0.;
        dt[j] = 0.;
    }

    for (Int_t j = 0; j < fNumAnodes; j++)
    {
        good_dt[j] = 0.;
    }

    for (Int_t i = 0; i < nHits; i++)
    {
        CalDat[i] = dynamic_cast<R3BMusicCalData*>(fMusicCalDataCA->At(i));
        anodeId = CalDat[i]->GetAnodeID();
        energyperanode[anodeId] = CalDat[i]->GetEnergy();
        dt[anodeId] = CalDat[i]->GetDTime();
    }

    Double_t nba = 0., theta = 0., Esum = 0.;
    // calculate truncated dE from 8 anodes, MUSIC
    fNumAnodesAngleFit = 0;
    for (Int_t j = 0; j < fNumAnodes; j++)
    {
        if (energyperanode[j] > 0. && fStatusAnodes[j] == 1)
        {
            Esum = Esum + energyperanode[j];
            good_dt[fNumAnodesAngleFit] = dt[j];
            fPosAnodes[fNumAnodesAngleFit] = fCal_Par->GetAnodePos(j + 1);
            fNumAnodesAngleFit++;
            nba++;
        }
    }

    if (Esum / nba > 0.)
    {
        fPosZ.Use(fNumAnodesAngleFit, fPosAnodes);
        TMatrixD A(fNumAnodesAngleFit, 2);
        TMatrixDColumn(A, 0) = 1.0;
        TMatrixDColumn(A, 1) = fPosZ;
        TDecompSVD svd(A);
        Bool_t ok;
        TVectorD dt_r;
        dt_r.Use(fNumAnodesAngleFit, good_dt);
        TVectorD c_svd_r = svd.Solve(dt_r, ok);
        theta = c_svd_r[1];

        Double_t zhit = fZ0 + fZ1 * TMath::Sqrt(Esum / nba);

        if (zhit > 0)
            // AddHitData(theta, zhit, Esum / nba);
            AddHitData(theta, zhit, Esum / nba, good_dt[3]); // for tracking
    }
    if (CalDat)
        delete[] CalDat;

    return;
}

// -----   Public method Reset   ------------------------------------------------
void R3BMusicCal2Hit::Reset()
{
    LOG(debug) << "Clearing R3BMusicHitData Structure";
    if (fMusicHitDataCA)
        fMusicHitDataCA->Clear();
}
// -----   Private method AddHitData  --------------------------------------------
R3BMusicHitData* R3BMusicCal2Hit::AddHitData(Double_t theta, Double_t charge_z, Double_t ene_ave, Double_t good_dt)
{
    // It fills the R3BMusicHitData
    TClonesArray& clref = *fMusicHitDataCA;
    Int_t size = clref.GetEntriesFast();
    return new (clref[size]) R3BMusicHitData(theta, charge_z, ene_ave, good_dt);
}

ClassImp(R3BMusicCal2Hit);
