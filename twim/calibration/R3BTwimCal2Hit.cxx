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
// -----           R3BTwimCal2Hit source file              -----
// -----    Created 30/11/19 by J.L. Rodriguez-Sanchez     -----
// -------------------------------------------------------------

// ROOT headers
#include <TCanvas.h>
#include <TClonesArray.h>
#include <TDecompSVD.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TMath.h>
#include <TMatrixD.h>
#include <TRandom.h>
#include <TVector3.h>
#include <TVectorD.h>
#include <cmath>
#include <cstdint>

// FAIR headers
#include <FairLogger.h>
#include <FairRootManager.h>
#include <FairRunAna.h>
#include <FairRuntimeDb.h>

// Twim headers
#include "R3BEventHeader.h"
#include "R3BLogger.h"
#include "R3BTwimCal2Hit.h"
#include "R3BTwimCalData.h"
#include "R3BTwimHitData.h"
#include "R3BTwimHitPar.h"
#ifdef SOFIA
#include "R3BSofTofWHitData.h"
#endif

// R3BTwimCal2Hit: Default Constructor --------------------------
R3BTwimCal2Hit::R3BTwimCal2Hit()
    : R3BTwimCal2Hit("R3BTwimCal2Hit", 1)
{
}

// R3BTwimCal2HitPar: Standard Constructor --------------------------
R3BTwimCal2Hit::R3BTwimCal2Hit(const char* name, Int_t iVerbose)
    : FairTask(name, iVerbose)
    , fNumSec(4)
    , fNumAnodes(16)
    , fNumAnodesAngleFit(0)
    , fNumParamsTof(3)
    , fNumParams(5)
    , fMaxEnergyperanode(65535)
    , CalZParams(NULL)
    , CalZTofParams(NULL)
    , fCal_Par(NULL)
    , fTwimHitDataCA(NULL)
    , fTwimCalDataCA(NULL)
    , fHitItemsTofW(NULL)
    , fOnline(kFALSE)
    , fExpId(0)
    , fTpat(-1)
    , fDebug(false)
{
}

// Virtual R3BTwimCal2Hit: Destructor
R3BTwimCal2Hit::~R3BTwimCal2Hit()
{
    R3BLOG(debug1, "");
    if (fTwimHitDataCA)
        delete fTwimHitDataCA;
    if (fHitItemsTofW)
        delete fHitItemsTofW;
}

void R3BTwimCal2Hit::SetParContainers()
{
    // Parameter Container
    // Reading TwimCalPar from FairRuntimeDb
    auto* rtdb = FairRuntimeDb::instance();
    R3BLOG_IF(fatal, !rtdb, "FairRuntimeDb not found");

    fCal_Par = dynamic_cast<R3BTwimHitPar*>(rtdb->getContainer("twimHitPar"));
    if (!fCal_Par)
    {
        R3BLOG(error, "Couldn't get handle on twimHitPar container");
    }
    else
    {
        R3BLOG(info, "twimHitPar container open");
    }
    return;
}

void R3BTwimCal2Hit::SetParameter()
{
    //--- Parameter Container ---
    fNumSec = fCal_Par->GetNumSec();        // Number of Sections
    fNumAnodes = fCal_Par->GetNumAnodes();  // Number of anodes
    fNumParams = fCal_Par->GetNumParZFit(); // Number of Parameters

    R3BLOG(info, "Nb sections: " << fNumSec);
    R3BLOG(info, "Nb anodes: " << fNumAnodes);

    // Anodes that don't work set to zero
    for (Int_t s = 0; s < fNumSec; s++)
        for (Int_t i = 0; i < fNumAnodes; i++)
        {
            StatusAnodes[s][i] = fCal_Par->GetInUse(s + 1, i + 1);
        }

    if (fHitItemsTofW)
    {
        R3BLOG(info, "Nb parameters for charge-Z vs Tof correction: " << fNumParamsTof);
        CalZTofParams = new TArrayF();
        Int_t size_tof = fNumSec * fNumParamsTof;
        CalZTofParams->Set(size_tof);
        CalZTofParams = fCal_Par->GetZTofHitPar(); // Array with the Cal parameters
    }

    R3BLOG(info, "Nb parameters for charge-Z: " << fNumParams);
    CalZParams = new TArrayF();
    Int_t array_size = fNumSec * fNumParams;
    CalZParams->Set(array_size);
    CalZParams = fCal_Par->GetZHitPar(); // Array with the Cal parameters

    // Parameters detector
    for (Int_t s = 0; s < fNumSec; s++)
    {
        // Parameters detector
        fEmean_tof[s] = fCal_Par->GetEmean_tof(s + 1);
        fEmean_dt[s] = fCal_Par->GetEmean_dt(s + 1);
        if (fNumParams == 2)
        {
            R3BLOG(info,
                   "Parameters for charge-Z:" << CalZParams->GetAt(s * fNumParams) << " : "
                                              << CalZParams->GetAt(s * fNumParams + 1));
            fZ0[s] = CalZParams->GetAt(s * fNumParams);
            fZ1[s] = CalZParams->GetAt(s * fNumParams + 1);
        }
        else if (fNumParams == 3)
        {
            if (fHitItemsTofW)
            {
                R3BLOG(info,
                       "Parameters for charge-Z vs tof:" << CalZParams->GetAt(s * fNumParams) << " : "
                                                         << CalZParams->GetAt(s * fNumParams + 1) << " : "
                                                         << CalZParams->GetAt(s * fNumParams + 2));
            }
            fZ0[s] = CalZParams->GetAt(s * fNumParams);
            fZ1[s] = CalZParams->GetAt(s * fNumParams + 1);
            fZ2[s] = CalZParams->GetAt(s * fNumParams + 2);
        }
        else if (fNumParams == 4)
        {
            if (fHitItemsTofW)
            {
                R3BLOG(info,
                       "Parameters for charge-Z vs tof:" << CalZParams->GetAt(s * fNumParams) << " : "
                                                         << CalZParams->GetAt(s * fNumParams + 1) << " : "
                                                         << CalZParams->GetAt(s * fNumParams + 2) << " : "
                                                         << CalZParams->GetAt(s * fNumParams + 3));
            }
            fZ0[s] = CalZParams->GetAt(s * fNumParams);
            fZ1[s] = CalZParams->GetAt(s * fNumParams + 1);
            fZ2[s] = CalZParams->GetAt(s * fNumParams + 2);
            fZ3[s] = CalZParams->GetAt(s * fNumParams + 3);
        }
        else if (fNumParams == 5)
        {
            if (fHitItemsTofW)
            {
                R3BLOG(info,
                       "Parameters for charge-Z vs tof:"
                           << CalZParams->GetAt(s * fNumParams) << " : " << CalZParams->GetAt(s * fNumParams + 1)
                           << " : " << CalZParams->GetAt(s * fNumParams + 2) << " : "
                           << CalZParams->GetAt(s * fNumParams + 3) << " : " << CalZParams->GetAt(s * fNumParams + 4));
            }
            fZ0[s] = CalZParams->GetAt(s * fNumParams);
            fZ1[s] = CalZParams->GetAt(s * fNumParams + 1);
            fZ2[s] = CalZParams->GetAt(s * fNumParams + 2);
            fZ3[s] = CalZParams->GetAt(s * fNumParams + 3);
            fZ4[s] = CalZParams->GetAt(s * fNumParams + 4);
        }
        else
            R3BLOG(warn, "Parameters for charge-Z cannot be used here, number of parameters: " << fNumParams);
    }
    return;
}

// -----   Public method Init   --------------------------------------------
InitStatus R3BTwimCal2Hit::Init()
{
    R3BLOG(info, "");
    auto* frm = FairRootManager::Instance();
    R3BLOG_IF(fatal, !frm, "FairRootManager not found.");

    header = dynamic_cast<R3BEventHeader*>(frm->GetObject("EventHeader."));
    R3BLOG_IF(error, !header, "EventHeader. not found.");

    if (fExpId == 0) // Obtain global ExpId if it's not set locally.
    {
        fExpId = header->GetExpId();
        R3BLOG(info, "fExpId: " << fExpId);
    }

    // INPUT DATA
    // get access to cal data of the Twim
    fTwimCalDataCA = dynamic_cast<TClonesArray*>(frm->GetObject("TwimCalData"));
    R3BLOG_IF(fatal, !fTwimCalDataCA, "TwimCalData not found");

    // get access to hit data of the Tof-Wall
    fHitItemsTofW = dynamic_cast<TClonesArray*>(frm->GetObject("TofWHitData"));
    R3BLOG_IF(warn, !fHitItemsTofW, "TofWHitData not found");

    // OUTPUT DATA
    fTwimHitDataCA = new TClonesArray("R3BTwimHitData");
    frm->Register("TwimHitData", "TWIM_Hit", fTwimHitDataCA, !fOnline);
    fTwimHitDataCA->Clear();

    if (fDebug)
    {
        dx_canvas = new TCanvas*[fNumSec];
        dx_all_canvas = new TCanvas*[fNumSec];
        dx_vs_pos_canvas = new TCanvas*[fNumSec];
        char Name1c[255];
        char Name2c[255];
        char Name3c[255];

        for (Int_t i = 0; i < fNumSec; i++)
        {
            sprintf(Name1c, "dx_canvas_sec_%d", i + 1);
            dx_canvas[i] = new TCanvas(Name1c, Name1c, 10, 10, 800, 700);
            dx_canvas[i]->Divide(4, 4);

            sprintf(Name2c, "dx_vs_pos_canvas_sec_%d", i + 1);
            dx_vs_pos_canvas[i] = new TCanvas(Name2c, Name2c, 10, 10, 800, 700);
            dx_vs_pos_canvas[i]->Divide(4, 4);

            sprintf(Name3c, "dx_all_canvas_sec_%d", i + 1);
            dx_all_canvas[i] = new TCanvas(Name3c, Name3c, 10, 10, 800, 700);
        }

        char Name1[255];
        char Name2[255];
        char Name1_legend[255];
        char Name2_title[255];
        dx = new TH1F*[fNumSec * fNumAnodes];
        dx_vs_pos = new TH2F*[fNumSec * fNumAnodes];
        for (Int_t s = 0; s < fNumSec; s++)
        {
            for (Int_t i = 0; i < fNumAnodes; i++)
            {
                sprintf(Name1, "dx_sec_%d_anode_%d", s + 1, i + 1);
                sprintf(Name1_legend, "Anode %d", i + 1);
                sprintf(Name2, "dx_vs_pos_sec_%d_anode_%d", s + 1, i + 1);
                sprintf(Name2_title, "Section %d, Anode %d", s + 1, i + 1);
                dx[s * fNumAnodes + i] = new TH1F(Name1, Name1_legend, 1000, -1, 1);
                dx[s * fNumAnodes + i]->GetXaxis()->SetTitle("#deltax [mm]");
                dx[s * fNumAnodes + i]->GetYaxis()->SetTitle("Counts");
                if (s < 2)
                    dx_vs_pos[s * fNumAnodes + i] = new TH2F(Name2, Name2_title, 1000, 0, 100, 1000, -1, 1);
                else
                    dx_vs_pos[s * fNumAnodes + i] = new TH2F(Name2, Name2_title, 1000, -100, 0, 1000, -1, 1);
                dx_vs_pos[s * fNumAnodes + i]->GetXaxis()->SetTitle("X [mm]");
                dx_vs_pos[s * fNumAnodes + i]->GetYaxis()->SetTitle("#deltax [mm]");
                dx_canvas[s]->cd(i + 1);
                dx[s * fNumAnodes + i]->Draw("");
                dx_vs_pos_canvas[s]->cd(i + 1);
                dx_vs_pos[s * fNumAnodes + i]->Draw("ZCOL");
            }
        }

        for (Int_t s = 0; s < fNumSec; s++)
        {
            dx_all_canvas[s]->cd();
            dx[s * fNumAnodes + 3]->SetLineColor(46);
            dx[s * fNumAnodes + 3]->Draw("");
            for (Int_t i = 4; i < 13; i++)
            {
                dx[s * fNumAnodes + i]->SetLineColor(i - 3);
                dx[s * fNumAnodes + i]->Draw("SAME");
            }
        }
    }

    SetParameter();
    return kSUCCESS;
}

// -----   Public method ReInit   ----------------------------------------------
InitStatus R3BTwimCal2Hit::ReInit()
{
    SetParContainers();
    SetParameter();
    return kSUCCESS;
}

// -----   Public method Execution   --------------------------------------------
void R3BTwimCal2Hit::Exec(Option_t* option)
{
    if ((fTpat >= 0) && (header) && ((header->GetTpat() & fTpat) != fTpat))
        return;
    // At the moment we will use the expid to select the reconstruction
    // this should be changed in the future because expid is not necessary
    if (fExpId == 91 || fExpId == 118 || fExpId == 444 || fExpId == 467 || fExpId == 509 || fExpId == 522)
    {
        S467();
    }
    else if (fExpId == 4551)
    {
        S4551();
    }
    else if (fExpId == 455)
    {
        S455();
    }

    return;
}
void R3BTwimCal2Hit::S4551()
{
    // Reset entries in output arrays, local arrays
    Reset();

#ifdef SOFIA

    Int_t nHitTwim = fTwimCalDataCA->GetEntriesFast();
    Int_t nHitTofW = fHitItemsTofW->GetEntriesFast();

    if (nHitTofW == 2 && nHitTwim > 16 && nHitTwim < 33)
    {

        R3BTwimCalData** CalDat = new R3BTwimCalData*[nHitTwim];
        R3BSofTofWHitData** HitTofW = new R3BSofTofWHitData*[nHitTofW];

        // #ifdef SOFIA
        //  Fill TofW Hit data
        TVector3 postof[2];
        postof[0].SetXYZ(-1000., -1000., 0.);
        postof[1].SetXYZ(-1000., -1000., 0.);
        Double_t tof[4];

        for (Int_t i = 0; i < nHitTofW; i++)
        {
            HitTofW[i] = (R3BSofTofWHitData*)(fHitItemsTofW->At(i));
            if (i == 0)
            {
                postof[0].SetXYZ(HitTofW[i]->GetPaddle(), HitTofW[i]->GetY(), HitTofW[i]->GetTof());
            }
            else
            {
                postof[1].SetXYZ(HitTofW[i]->GetPaddle(), HitTofW[i]->GetY(), HitTofW[i]->GetTof());
            }
        }

        for (Int_t i = 0; i < fNumSec; i++)
        {
            tof[i] = 0;
        }

        if (postof[0].Y() < 0 && postof[1].Y() > 0)
        {
            tof[0] = postof[0].Z();
            tof[1] = postof[1].Z();
            tof[2] = postof[1].Z();
            tof[3] = postof[0].Z();
        }

        else if (postof[0].Y() > 0 && postof[1].Y() < 0)
        {
            tof[0] = postof[1].Z();
            tof[1] = postof[0].Z();
            tof[2] = postof[0].Z();
            tof[3] = postof[1].Z();
        }

        else if (postof[0].Y() > 0 && postof[1].Y() > 0)
        {
            if (postof[0].X() > postof[1].X())
            {
                tof[0] = postof[1].Z();
                tof[1] = postof[1].Z();
                tof[2] = postof[0].Z();
                tof[3] = postof[0].Z();
            }
            else if (postof[0].X() < postof[1].X())
            {
                tof[0] = postof[0].Z();
                tof[2] = postof[1].Z();
                tof[1] = postof[0].Z();
                tof[3] = postof[1].Z();
            }
        }

        else if (postof[0].Y() < 0 && postof[1].Y() < 0)
        {
            if (postof[0].X() > postof[1].X())
            {
                tof[0] = postof[1].Z();
                tof[1] = postof[1].Z();
                tof[2] = postof[0].Z();
                tof[3] = postof[0].Z();
            }
            else if (postof[0].X() < postof[1].X())
            {
                tof[0] = postof[0].Z();
                tof[1] = postof[0].Z();
                tof[2] = postof[1].Z();
                tof[3] = postof[1].Z();
            }
        }

        Int_t secId = 0, anodeId = 0;
        Double_t energyperanode[fNumSec][fNumAnodes];
        Double_t dt[fNumSec][fNumAnodes];
        double good_dt[fNumAnodes];
        Int_t nbdet = 0;

        for (Int_t j = 0; j < fNumAnodes; j++)
        {
            good_dt[j] = std::nan("");
            for (Int_t i = 0; i < fNumSec; i++)
            {
                energyperanode[i][j] = 0.;
                dt[i][j] = -1000.;
            }
        }

        for (Int_t i = 0; i < nHitTwim; i++)
        {
            CalDat[i] = (R3BTwimCalData*)(fTwimCalDataCA->At(i));
            secId = CalDat[i]->GetSecID() - 1;
            anodeId = CalDat[i]->GetAnodeID() - 1;
            if (energyperanode[secId][anodeId] == 0)
            {
                energyperanode[secId][anodeId] = CalDat[i]->GetEnergy();
                dt[secId][anodeId] = CalDat[i]->GetDTime();
            }
        }

        // calculate truncated dE from 16 anodes, Twim-MUSIC
        for (Int_t i = 0; i < fNumSec; i++)
        {
            Double_t nba = 0, offset = 0., theta = -5000., Esum = 0.;
            fNumAnodesAngleFit = 0;
            Double_t dt_ref = 0.;
            Double_t Xfit[fNumAnodes];

            for (Int_t j = 0; j < fNumAnodes; j++)
            {
                Xfit[j] = 0;
                fIndex[fNumAnodesAngleFit] = j;
                if (j > 2 && j < 13 && energyperanode[i][j] > 0 && energyperanode[i][j] < fMaxEnergyperanode &&
                    StatusAnodes[i][j] == 1 && dt[i][j] > -100. && dt[i][j] < 100.)
                {
                    good_dt[fNumAnodesAngleFit] = dt[i][j];
                    Esum += energyperanode[i][j];
                    if (dt[i][j] > -1000.)
                    {
                        if (j == 7)
                            dt_ref = dt[i][j];
                    }
                    fPosAnodes[fNumAnodesAngleFit] = fCal_Par->GetAnodePos(j + 1);
                    fNumAnodesAngleFit++;
                    nba++;
                }
            } // loop fNumAnodes
            if (nba > 3 && (Esum / nba) > 0.)
            {
                if (fNumAnodesAngleFit > 8)
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
                    offset = c_svd_r[0];
                    theta = c_svd_r[1];
                    for (Int_t y = 0; y < fNumAnodesAngleFit; y++)
                    {
                        Xfit[y] = c_svd_r[0] + c_svd_r[1] * fPosAnodes[y];
                        if (fDebug)
                        {
                            dx[i * fNumAnodes + fIndex[y]]->Fill(Xfit[y] - good_dt[y]);
                            dx_vs_pos[i * fNumAnodes + fIndex[y]]->Fill(Xfit[y], Xfit[y] - good_dt[y]);
                        }
                    }
                }
                if (CalZTofParams)
                {
                    Double_t Esum_mean = Esum / nba;
                    Esum_mean = fEmean_tof[i] * Esum_mean /
                                (CalZTofParams->GetAt(i * fNumParamsTof) +
                                 CalZTofParams->GetAt(i * fNumParamsTof + 1) * tof[i] +
                                 CalZTofParams->GetAt(i * fNumParamsTof + 2) * tof[i] * tof[i]);
                    TSpline3* spl = fCal_Par->GetSpline(i + 1);
                    Esum_mean = fEmean_dt[i] * Esum_mean / spl->Eval(dt_ref);
                    Double_t zhit = fZ0[i] + fZ1[i] * TMath::Sqrt(Esum_mean) + fZ2[i] * Esum_mean +
                                    fZ3[i] * TMath::Power(Esum_mean, 3. / 2.) + fZ4[i] * TMath::Power(Esum_mean, 2.);
                    if (zhit > 0)
                        AddHitData(i + 1, theta, zhit, dt_ref, offset, Esum_mean);
                }
                else
                {
                    Double_t Esum_mean = Esum / nba;
                    Double_t zhit = fZ0[i] + fZ1[i] * TMath::Sqrt(Esum_mean) + fZ2[i] * Esum_mean +
                                    fZ3[i] * TMath::Power(Esum_mean, 3. / 2.) + fZ4[i] * TMath::Power(Esum_mean, 2.);
                    if (zhit > 0)
                        AddHitData(i + 1, theta, zhit, dt_ref, offset, Esum_mean);
                }
            } // loop nba>8
        }     // loop NumSec
        if (CalDat)
            delete[] CalDat;
        if (HitTofW)
            delete[] HitTofW;
    } // if nHitTofW == 2

#endif
    return;
}

// -----   Private method s455       --------------------------------------------
void R3BTwimCal2Hit::S455()
{
    // Reset entries in output arrays, local arrays
    Reset();

    Int_t nHits = fTwimCalDataCA->GetEntriesFast();
    if (!nHits)
        return;

    R3BTwimCalData** CalDat = new R3BTwimCalData*[nHits];

    Int_t secId = 0, anodeId = 0;
    Double_t energyperanode[fNumSec][fNumAnodes];
    Double_t dt[fNumSec][fNumAnodes];
    double good_dt[fNumAnodes];
    Int_t nbdet = 0;

    for (Int_t j = 0; j < fNumAnodes; j++)
    {
        good_dt[j] = std::nan("");
        for (Int_t i = 0; i < fNumSec; i++)
        {
            energyperanode[i][j] = 0.;
            dt[i][j] = 0.;
        }
    }

    // std::cout<<"Event " <<std::endl;
    for (Int_t i = 0; i < nHits; i++)
    {
        CalDat[i] = dynamic_cast<R3BTwimCalData*>(fTwimCalDataCA->At(i));
        secId = CalDat[i]->GetSecID() - 1;
        anodeId = CalDat[i]->GetAnodeID() - 1;
        if (energyperanode[secId][anodeId] == 0)
        {
            energyperanode[secId][anodeId] = CalDat[i]->GetEnergy();
            dt[secId][anodeId] = CalDat[i]->GetDTime();
            // std::cout<< secId <<" "<<anodeId <<" "<< CalDat[i]->GetEnergy() <<std::endl;
        }
    }

#ifdef SOFIA
    // Fill TofW Hit data
    Double_t tof[2] = { 0., 0. };
    Int_t padid[2] = { 0, 0 };
    if (fHitItemsTofW && fHitItemsTofW->GetEntriesFast() > 0)
    {
        nHits = fHitItemsTofW->GetEntriesFast();
        for (Int_t ihit = 0; ihit < nHits; ihit++)
        {
            R3BSofTofWHitData* hit = (R3BSofTofWHitData*)fHitItemsTofW->At(ihit);
            if (!hit)
                continue;
            if (padid[0] == 0)
            {
                padid[0] = hit->GetPaddle();
                tof[0] = hit->GetTof();
            }
            else
            {
                if (hit->GetPaddle() > padid[0] && hit->GetPaddle() - padid[0] > 1)
                {
                    padid[1] = hit->GetPaddle();
                    tof[1] = hit->GetTof(); // right
                }
                else if (hit->GetPaddle() - padid[0] < -1)
                {
                    tof[1] = tof[0]; // right
                    padid[1] = padid[0];
                    tof[0] = hit->GetTof(); // new left
                    padid[0] = hit->GetPaddle();
                }
            }
        }
    }

    // calculate truncated dE from 16 anodes, Twim-MUSIC
    for (Int_t i = 0; i < fNumSec; i++)
    {
        Double_t nba = 0, offset = 0., theta = -5000., Esum = 0.;
        fNumAnodesAngleFit = 0;
        for (Int_t j = 0; j < fNumAnodes; j++)
        {
            // if(i==0){
            // if (((energyperanode[i][j] > 0 && energyperanode[i+1][j] == 0) || (energyperanode[i][j] == 0 &&
            // energyperanode[i+1][j] > 0)) && energyperanode[i][j] < 65535 && StatusAnodes[i][j] == 1 &&
            // energyperanode[i+1][j] < 65535 && StatusAnodes[i+1][j] == 1)

            if (energyperanode[i][j] > 0 && energyperanode[i][j] < fMaxEnergyperanode && StatusAnodes[i][j] == 1)
            {
                Esum += energyperanode[i][j];
                if (dt[i][j] > 0.)
                    good_dt[fNumAnodesAngleFit] = dt[i][j];
                else if (dt[i + 1][j] > 0.)
                    good_dt[fNumAnodesAngleFit] = dt[i + 1][j];
                fPosAnodes[fNumAnodesAngleFit] = fCal_Par->GetAnodePos(j + 1);
                fNumAnodesAngleFit++;
                nba++;
                // std::cout<< i <<" "<< j <<" "<< energyperanode[i][j] <<std::endl;
            }
            /*}else
            if ( energyperanode[i][j] > 0 && energyperanode[i][j] < 65535 && StatusAnodes[i][j] == 1 &&
            energyperanode[i+1][j] < 65535 && StatusAnodes[i+1][j] == 1)
            {
                Esum = Esum + energyperanode[i][j];// + energyperanode[i + 1][j];
                if (dt[i][j] > 0.)
                    good_dt[fNumAnodesAngleFit] = dt[i][j];
                else if (dt[i + 1][j] > 0.)
                    good_dt[fNumAnodesAngleFit] = dt[i + 1][j];
                fPosAnodes[fNumAnodesAngleFit] = fCal_Par->GetAnodePos(j + 1);
                fNumAnodesAngleFit++;
                nba++;
            }
            */
        }

        if (nba > 8 && (Esum / nba) > 0.)
        {
            if (fNumAnodesAngleFit > 8)
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
                offset = c_svd_r[0];
                theta = c_svd_r[1];
            }
            if (CalZTofParams)
            {
                Double_t Esum_mean = Esum / nba;
                if (i < 2)
                    Esum_mean = Esum_mean - (CalZTofParams->GetAt(i * fNumParams) +
                                             CalZTofParams->GetAt(i * fNumParams + 1) * tof[0] +
                                             CalZTofParams->GetAt(i * fNumParams + 2) * tof[0] * tof[0]);
                else
                    Esum_mean = Esum_mean - (CalZTofParams->GetAt(i * fNumParams) +
                                             CalZTofParams->GetAt(i * fNumParams + 1) * tof[1] +
                                             CalZTofParams->GetAt(i * fNumParams + 2) * tof[1] * tof[1]);

                Double_t zhit = fZ0[i] + fZ1[i] * TMath::Sqrt(Esum_mean) + fZ2[i] * Esum_mean;
                // if (zhit > 0 && theta > -5000.)
                if (zhit > 0 && !TMath::IsNaN(good_dt[7]))
                    AddHitData(i + 1, theta, zhit, good_dt[7], offset, Esum_mean);
            }
            else
            {
                Double_t Esum_mean = Esum / nba;
                Double_t zhit = fZ0[i] + fZ1[i] * TMath::Sqrt(Esum_mean) + fZ2[i] * Esum_mean;
                // if (zhit > 0 && theta > -5000.)
                if (zhit > 0 && !TMath::IsNaN(good_dt[7]))
                    AddHitData(i + 1, theta, zhit, good_dt[7], offset, Esum_mean);
            }
        }
        // i = i + 2;
    }
#endif
    if (CalDat)
        delete[] CalDat;
    return;
}

// -----   Private method s467       --------------------------------------------
void R3BTwimCal2Hit::S467()
{
    // Reset entries in output arrays, local arrays
    Reset();

    Int_t nHits = fTwimCalDataCA->GetEntriesFast();
    if (nHits == 0)
        return;

    R3BTwimCalData** CalDat = new R3BTwimCalData*[nHits];

    Int_t secId, anodeId;
    Double_t energyperanode[fNumSec][fNumAnodes];
    Double_t dt[fNumSec][fNumAnodes];
    double good_dt[fNumAnodes];
    Int_t nbdet = 0;

    for (Int_t j = 0; j < fNumAnodes; j++)
    {
        good_dt[j] = std::nan("");
        for (Int_t i = 0; i < fNumSec; i++)
        {
            energyperanode[i][j] = 0;
            dt[i][j] = 0.;
        }
    }

    for (Int_t i = 0; i < nHits; i++)
    {
        CalDat[i] = dynamic_cast<R3BTwimCalData*>(fTwimCalDataCA->At(i));
        secId = CalDat[i]->GetSecID() - 1;
        anodeId = CalDat[i]->GetAnodeID() - 1;
        energyperanode[secId][anodeId] = CalDat[i]->GetEnergy();
        dt[secId][anodeId] = CalDat[i]->GetDTime();
    }

    // calculate truncated dE from 16 anodes, Twim-MUSIC
    for (Int_t i = 0; i < fNumSec; i++)
    {
        Double_t nba = 0., offset = 0., theta = -5000., Esum = 0.;
        fNumAnodesAngleFit = 0;
        for (Int_t j = 0; j < fNumAnodes; j++)
        {
            if (energyperanode[i][j] > 0 && energyperanode[i][j] < 8192 && StatusAnodes[i][j] == 1)
            {
                Esum = Esum + energyperanode[i][j];
                good_dt[fNumAnodesAngleFit] = dt[i][j];
                fPosAnodes[fNumAnodesAngleFit] = fCal_Par->GetAnodePos(j + 1);
                fNumAnodesAngleFit++;
                nba++;
            }
        }

        if (nba > 0 && (Esum / nba) > 0.)
        {
            if (fNumAnodesAngleFit > 4)
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
                offset = c_svd_r[0];
                theta = c_svd_r[1];
            }
            double zhit =
                fZ0[i] + fZ1[i] * TMath::Sqrt(Esum / nba) + fZ2[i] * TMath::Sqrt(Esum / nba) * TMath::Sqrt(Esum / nba);
            if (zhit > 0 && theta > -5000. && good_dt[7] != std::nan("") && nba > 0) // NOLINT
                AddHitData(i + 1, theta, zhit, good_dt[7], offset, Esum / nba);
        }
    }

    if (CalDat)
        delete[] CalDat;
    return;
}

// -----   Public method Reset   ------------------------------------------------
void R3BTwimCal2Hit::Reset()
{
    R3BLOG(debug1, "Clearing TwimHitData Structure");
    if (fTwimHitDataCA)
        fTwimHitDataCA->Clear();
}

// -----   Private method AddHitData  --------------------------------------------
R3BTwimHitData* R3BTwimCal2Hit::AddHitData(UShort_t secid,
                                           Double_t theta,
                                           Double_t charge_z,
                                           Double_t xpos,
                                           Double_t offset,
                                           Double_t ene_ave)
{
    // It fills the R3BTwimHitData
    TClonesArray& clref = *fTwimHitDataCA;
    Int_t size = clref.GetEntriesFast();
    return new (clref[size]) R3BTwimHitData(secid, theta, charge_z, xpos, offset, ene_ave);
}
void R3BTwimCal2Hit::FinishTask()
{
    if (fDebug)
    {
        for (Int_t s = 0; s < fNumSec; s++)
        {
            dx_canvas[s]->Write();
            dx_all_canvas[s]->Write();
            dx_vs_pos_canvas[s]->Write();
        }
    }
}

ClassImp(R3BTwimCal2Hit)
