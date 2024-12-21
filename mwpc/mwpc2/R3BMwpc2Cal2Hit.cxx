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

// --------------------------------------------------------------
// -----         R3BMwpc2Cal2Hit source file                -----
// -----   Created 11/10/19 by J.L. Rodriguez-Sanchez       -----
// -----      s455 method 17/06/22 by Antia GG              -----
// --------------------------------------------------------------

// ROOT headers
#include "TClonesArray.h"
#include "TMath.h"

// Fair headers
#include "FairLogger.h"
#include "FairRootManager.h"
#include "FairRunAna.h"
#include "FairRuntimeDb.h"

// MWPC headers
#include "R3BEventHeader.h"
#include "R3BLogger.h"
#include "R3BMwpc2Cal2Hit.h"
#include "R3BMwpcCalData.h"
#include "R3BMwpcHitData.h"

// R3BMwpc2Cal2Hit: Default Constructor --------------------------
R3BMwpc2Cal2Hit::R3BMwpc2Cal2Hit()
    : R3BMwpc2Cal2Hit("R3BMwpc2Cal2Hit", 1)
{
}

// R3BMwpc2Cal2Hit: Standard Constructor --------------------------
R3BMwpc2Cal2Hit::R3BMwpc2Cal2Hit(const char* name, Int_t iVerbose)
    : FairTask(name, iVerbose)
    , fMwpcCalDataCA(NULL)
    , fMwpcHitDataCA(NULL)
    , fwx(3.125)   // in mm
    , fwy(5.000)   // in mm
    , fSize(200.0) // in mm
    , fOnline(kFALSE)
    , fExpId(0)
{
}

// Virtual R3BMwpc2Cal2Hit: Destructor
R3BMwpc2Cal2Hit::~R3BMwpc2Cal2Hit()
{
    R3BLOG(debug1, "Destructor");
    if (fMwpcCalDataCA)
        delete fMwpcCalDataCA;
    if (fMwpcHitDataCA)
        delete fMwpcHitDataCA;
}

// -----   Public method Init   --------------------------------------------
InitStatus R3BMwpc2Cal2Hit::Init()
{
    R3BLOG(info, "");

    // INPUT DATA
    FairRootManager* rootManager = FairRootManager::Instance();
    if (!rootManager)
    {
        R3BLOG(fatal, "FairRootManager not found");
        return kFATAL;
    }

    header = dynamic_cast<R3BEventHeader*>(rootManager->GetObject("EventHeader."));
    if (!header)
        header = dynamic_cast<R3BEventHeader*>(rootManager->GetObject("R3BEventHeader"));

    fMwpcCalDataCA = dynamic_cast<TClonesArray*>(rootManager->GetObject("Mwpc2CalData"));
    if (!fMwpcCalDataCA)
    {
        R3BLOG(fatal, "Mwpc2CalData not found");
        return kFATAL;
    }

    // OUTPUT DATA
    // Hit data
    fMwpcHitDataCA = new TClonesArray("R3BMwpcHitData");
    rootManager->Register("Mwpc2HitData", "MWPC2 Hit", fMwpcHitDataCA, !fOnline);
    fMwpcHitDataCA->Clear();
    return kSUCCESS;
}

// -----   Public method ReInit   ----------------------------------------------
InitStatus R3BMwpc2Cal2Hit::ReInit() { return kSUCCESS; }

// ------ Some functions ---------------

// bool R3BMwpc2Cal2Hit::sortPairsmwpc2(const pair<Double_t, Int_t> &x, const pair<Double_t, Int_t> &y)
bool sortPairsmwpc2(const pair<Double_t, Int_t>& x, const pair<Double_t, Int_t>& y) { return x.first > y.first; }

Double_t R3BMwpc2Cal2Hit::GetPositionX(Double_t qmax, Int_t padmax, Double_t qleft, Double_t qright)
{
    Double_t a3 = TMath::Pi() * fwx / (TMath::ACosH(0.5 * (TMath::Sqrt(qmax / qleft) + TMath::Sqrt(qmax / qright))));
    // Double_t a2 = gRandom->Uniform(-fwx / 2,fwx / 2);
    Double_t a2 = (a3 / TMath::Pi()) * TMath::ATanH((TMath::Sqrt(qmax / qleft) - TMath::Sqrt(qmax / qright)) /
                                                    (2 * TMath::SinH(TMath::Pi() * fwx / a3)));

    return (-1. * padmax * fwx + (fSize / 2) - (fwx / 2) - a2); // Left is positive and right negative
}

/* ----   Protected method to obtain the position Y ---- */
Double_t R3BMwpc2Cal2Hit::GetPositionY(Double_t qmax, Int_t padmax, Double_t qdown, Double_t qup)
{
    Double_t a2 = 0;
    if (qdown != 0 && qup != 0)
    {
        Double_t a3 = TMath::Pi() * fwy / (TMath::ACosH(0.5 * (TMath::Sqrt(qmax / qdown) + TMath::Sqrt(qmax / qup))));
        // Double_t a2 = gRandom->Uniform(-fwy / 2, fwy / 2);
        a2 = (a3 / TMath::Pi()) * TMath::ATanH((TMath::Sqrt(qmax / qdown) - TMath::Sqrt(qmax / qup)) /
                                               (2 * TMath::SinH(TMath::Pi() * fwy / a3)));
    }

    return (padmax * fwy - (fSize / 2) + (fwy / 2) + a2);
}

/* ----   Public method Execution   ---- */
void R3BMwpc2Cal2Hit::Exec(Option_t* option)
{
    // At the moment we will use the expid to select the reconstruction
    // this should be changed in the future because expid is not necessary
    int expid = fExpId != 0 ? fExpId : header->GetExpId();
    if (expid == 455)
    {
        S455();
    }
    else
    {
        S467();
    }

    return;
}

// -----   Private method s467       --------------------------------------------
void R3BMwpc2Cal2Hit::S467()
{
    // Reset entries in output arrays, local arrays
    Reset();

    // Reading the Input -- Cal Data --
    Int_t nHits = fMwpcCalDataCA->GetEntriesFast();
    if (nHits == 0)
        return;

    // Data from cal level
    R3BMwpcCalData** calData;
    calData = new R3BMwpcCalData*[nHits];
    Int_t planeId;
    Int_t padId;
    Int_t padmx = -1, padmxu = -1, padmxd = -1, padmy = -1;
    Double_t q = 0., qmxu = 0., qmxd = 0., qleft = 0., qright = 0.;
    Double_t qmx = 0., qmy = 0., qdown = 0., qup = 0.;
    Double_t x = NAN, y = NAN;

    for (Int_t i = 0; i < Mw2PadsX; i++)
        fx[i] = 0;
    for (Int_t i = 0; i < Mw2PadsY; i++)
        fy[i] = 0;

    for (Int_t i = 0; i < nHits; i++)
    {
        calData[i] = dynamic_cast<R3BMwpcCalData*>(fMwpcCalDataCA->At(i));
        planeId = calData[i]->GetPlane();
        padId = calData[i]->GetPad() - 1; // From 0 to 63 for X down and up
        q = calData[i]->GetQ();

        // FIXME: in November this should be OK!
        if (planeId == 1 || planeId == 2)
            fx[padId] += q; // Xup+Xdown
        else
            fy[padId] = q;

        if (q > qmx && (planeId == 1 || planeId == 2))
        {
            qmx = q;
            padmx = padId;
        }
        if (q > qmy && planeId == 3)
        {
            qmy = q;
            padmy = padId;
        }
    }
    // Add Hit data ----
    if (padmx > 1 && padmy > 1 && padmx + 1 < Mw2PadsX && padmy + 1 < Mw2PadsY && qmx > 0 && qmy > 0)
    {
        // FIXME: in November this should be OK!
        // Obtain position X ----
        qleft = (Double_t)fx[padmx - 1];
        qright = (Double_t)fx[padmx + 1];
        if (qleft > 0 && qright > 0)
            x = GetPositionX(qmx, padmx, qleft, qright);
        // Obtain position Y ----
        qdown = fy[padmy - 1];
        qup = fy[padmy + 1];
        if (qdown > 0 && qup > 0)
            y = GetPositionY(qmy, padmy, qdown, qup);

        AddHitData(x, y);
    }

    if (calData)
        delete[] calData;
    return;
}

// -----   Private method s455       --------------------------------------------
void R3BMwpc2Cal2Hit::S455()
{
    // Reset entries in output arrays, local arrays
    Reset();
    // Reading the Input -- Cal Data --
    Int_t nHits = fMwpcCalDataCA->GetEntriesFast();
    R3BMwpcCalData** calData = new R3BMwpcCalData*[nHits];

    Int_t planeId = 0;
    Int_t padId = 0;
    Double_t q = 0;
    Int_t planex1 = 0, planex2 = 0;
    Int_t padmx1 = -1, padmy1 = -1, padmx2 = -1, padmy2 = -1;
    Int_t padmx_p1 = -1, padmx_p2 = -1;
    Int_t qmx_p1 = 0., qmx_p2 = 0.;
    Double_t qmx1 = 0., qmy1 = 0., qleft1 = 0., qright1 = 0., qdown1 = 0., qup1 = 0.;
    Double_t qmx2 = 0., qmy2 = 0., qleft2 = 0., qright2 = 0., qdown2 = 0., qup2 = 0.;
    Double_t x1 = 0., y1 = 0., x2 = 0., y2 = 0.;
    Int_t nx_p1 = 0, nx_p2 = 0, ny = 0;

    vector<pair<Double_t, Int_t>> QpadX_p1, QpadX_p2, QpadY;
    QpadX_p1.clear();
    QpadX_p2.clear();
    QpadY.clear();
    Bool_t xexists_p1 = false;
    Bool_t xexists_p2 = false;
    Bool_t yexists = false;
    // Double_t fx[Mw2PadsX], fy[Mw2PadsY];
    // cout << "NUEVO EVENTO" << endl;
    for (Int_t i = 0; i < Mw2PadsX; i++)
        fx_p1[i] = 0;
    for (Int_t i = 0; i < Mw2PadsX; i++)
        fx_p2[i] = 0;
    for (Int_t i = 0; i < Mw2PadsY; i++)
        fy[i] = 0;
    for (Int_t i = 0; i < nHits; i++)
    {
        calData[i] = dynamic_cast<R3BMwpcCalData*>(fMwpcCalDataCA->At(i));
        planeId = calData[i]->GetPlane();
        padId = calData[i]->GetPad() - 1;
        q = calData[i]->GetQ();
        // cout << "q = " << q << ", padId = " << padId << endl;// << ", planeId = " << planeId << endl;
        pair<Double_t, Int_t> hit_pair = make_pair(q, padId);
        if (planeId == 1)
        {
            fx_p1[padId] = q;
            QpadX_p1.push_back(hit_pair);
            // xexists_p1 = true;
            nx_p1 = nx_p1 + 1;
        }
        if (planeId == 2)
        {
            fx_p2[padId] = q;
            QpadX_p2.push_back(hit_pair);
            // xexists_p2 = true;
            nx_p2 = nx_p2 + 1;
        }
        if (planeId == 3)
        {
            fy[padId] = q;
            QpadY.push_back(hit_pair);
            // yexists = true;
            ny = ny + 1;
        }
    }

    if ((nx_p1 > 0 || nx_p2 > 0) && ny > 0)
    {

        if (nx_p1 > 0)
        {
            sort(QpadX_p1.begin(),
                 QpadX_p1.end(),
                 sortPairsmwpc2); // el vector se ordena por el primer elemento de cada par
            qmx_p1 = QpadX_p1[0].first;
            padmx_p1 = QpadX_p1[0].second;
            // cout << "qmx1 = " << qmx1 << ", padmx1 = " << padmx1 << endl;
            /*for (Int_t i=0; i<nx_p1; i++){
              cout << "Loop i = " << i << " q = " << QpadX_p1[i].first << endl;
            }*/
        }
        if (nx_p2 > 0)
        {
            sort(QpadX_p2.begin(),
                 QpadX_p2.end(),
                 sortPairsmwpc2); // el vector se ordena por el primer elemento de cada par
            qmx_p2 = QpadX_p2[0].first;
            padmx_p2 = QpadX_p2[0].second;
            // cout << "qmx1 = " << qmx1 << ", padmx1 = " << padmx1 << endl;
        }
        if (qmx_p1 > qmx_p2 && padmx_p1 + 1 < Mw2PadsX && padmx_p2 + 1 < Mw2PadsX)
        {
            qmx1 = qmx_p1;
            padmx1 = padmx_p1;
            qleft1 = fx_p1[padmx1 - 1];
            qright1 = fx_p1[padmx1 + 1];
            planex1 = 1;
            /*cout << "1>2" << endl;
            cout << "qmx1 = " << qmx1 << ", padmx1 = " << padmx1 << endl;
            cout << "qmx2 = " << qmx2 << ", padmx2 = " << padmx2 << endl;*/
            if (qmx1 > 10 && qleft1 > 0 && qright1 > 0)
            {
                x1 = GetPositionX(qmx1, padmx1, qleft1, qright1);
            }
            else
                x1 = -1000.;

            for (Int_t i = 0; i < nx_p1; i++)
            {
                if (QpadX_p1[i].second == padmx_p1 - 3 || QpadX_p1[i].second == padmx_p1 - 2 ||
                    QpadX_p1[i].second == padmx_p1 - 1 || QpadX_p1[i].second == padmx_p1 ||
                    QpadX_p1[i].second == padmx_p1 + 1 || QpadX_p1[i].second == padmx_p1 + 2 ||
                    QpadX_p1[i].second == padmx_p1 + 3)
                {
                    QpadX_p1[i].first = 0.;
                }
            }
            sort(QpadX_p1.begin(), QpadX_p1.end(), sortPairsmwpc2);
            qmx_p1 = QpadX_p1[0].first;
            padmx_p1 = QpadX_p1[0].second;

            if (qmx_p1 > qmx_p2 && padmx_p1 + 1 < Mw2PadsX)
            {
                qmx2 = qmx_p1;
                padmx2 = padmx_p1;
                qleft2 = fx_p1[padmx2 - 1];
                qright2 = fx_p1[padmx2 + 1];
                planex2 = 1;
                if (qmx2 > 10 && qleft2 > 0 && qright2 > 0)
                {
                    x2 = GetPositionX(qmx2, padmx2, qleft2, qright2);
                }
                else
                    x2 = -1000.;
            }
            else if (qmx_p2 > qmx_p1 && padmx_p2 + 1 < Mw2PadsX)
            {
                qmx2 = qmx_p2;
                padmx2 = padmx_p2;
                qleft2 = fx_p2[padmx2 - 1];
                qright2 = fx_p2[padmx2 + 1];
                planex2 = 2;
                if (qmx2 > 10 && qleft2 > 0 && qright2 > 0)
                {
                    x2 = GetPositionX(qmx2, padmx2, qleft2, qright2);
                }
                else
                    x2 = -1000.;
            }
        }

        else if (qmx_p2 > qmx_p1 && padmx_p1 + 1 < Mw2PadsX && padmx_p2 + 1 < Mw2PadsX)
        {

            qmx1 = qmx_p2;
            padmx1 = padmx_p2;
            qleft1 = fx_p2[padmx1 - 1];
            qright1 = fx_p2[padmx1 + 1];
            planex1 = 2;
            /*cout << "2>1" << endl;
            cout << "qmx1 = " << qmx1 << ", padmx1 = " << padmx1 << endl;
            cout << "qmx2 = " << qmx2 << ", padmx2 = " << padmx2 << endl;*/
            if (qmx1 > 10 && qleft1 > 0 && qright1 > 0)
            {
                x1 = GetPositionX(qmx1, padmx1, qleft1, qright1);
            }
            else
                x1 = -1000.;

            for (Int_t i = 0; i < nx_p2; i++)
            {
                if (QpadX_p2[i].second == padmx_p2 - 3 || QpadX_p2[i].second == padmx_p2 - 2 ||
                    QpadX_p2[i].second == padmx_p2 - 1 || QpadX_p2[i].second == padmx_p2 ||
                    QpadX_p2[i].second == padmx_p2 + 1 || QpadX_p2[i].second == padmx_p2 + 2 ||
                    QpadX_p2[i].second == padmx_p2 + 3)
                {
                    QpadX_p2[i].first = 0.;
                }
            }
            sort(QpadX_p2.begin(), QpadX_p2.end(), sortPairsmwpc2);
            qmx_p2 = QpadX_p2[0].first;
            padmx_p2 = QpadX_p2[0].second;

            if (qmx_p2 > qmx_p1 && padmx_p2 + 1 < Mw2PadsX)
            {
                qmx2 = qmx_p2;
                padmx2 = padmx_p2;
                qleft2 = fx_p2[padmx2 - 1];
                qright2 = fx_p2[padmx2 + 1];
                planex2 = 2;
                if (qmx2 > 10 && qleft2 > 0 && qright2 > 0)
                {
                    x2 = GetPositionX(qmx2, padmx2, qleft2, qright2);
                }
                else
                    x2 = -1000.;
            }

            else if (qmx_p1 > qmx_p2 && padmx_p1 + 1 < Mw2PadsX)
            {
                qmx2 = qmx_p1;
                padmx2 = padmx_p1;
                qleft2 = fx_p1[padmx2 - 1];
                qright2 = fx_p1[padmx2 + 1];
                planex2 = 1;
                if (qmx2 > 10 && qleft2 > 0 && qright2 > 0)
                {
                    x2 = GetPositionX(qmx2, padmx2, qleft2, qright2);
                }
                else
                    x2 = -1000.;
            }
        }

        if (ny > 0)
        {

            sort(QpadY.begin(), QpadY.end(), sortPairsmwpc2);
            qmy1 = QpadY[0].first;
            padmy1 = QpadY[0].second;
            /*for (Int_t i =0; i<ny; i++){
                cout << "1  i = " << i << ", q = " << QpadY[i].first << ", pad = " << QpadY[i].second << endl;
            }
            cout << "qmy1 = " << qmy1 << ", padmy1 = " << padmy1 << endl;*/

            if (padmy1 + 1 < Mw2PadsY)
            {
                // Obtain position Y for 1st charge ----
                qdown1 = fy[padmy1 - 1];
                qup1 = fy[padmy1 + 1];
                if (qmy1 > 10)
                { //&& qdown1 > 0 && qup1 > 0){
                    y1 = GetPositionY(qmy1, padmy1, qdown1, qup1);
                }
                else
                    y1 = -1000.;
            }

            for (Int_t i = 0; i < ny; i++)
            {
                if (QpadY[i].second == padmy1 - 3 || QpadY[i].second == padmy1 - 2 || QpadY[i].second == padmy1 - 1 ||
                    QpadY[i].second == padmy1 || QpadY[i].second == padmy1 + 1 || QpadY[i].second == padmy1 + 2 ||
                    QpadY[i].second == padmy1 + 3)
                {
                    QpadY[i].first = 0.;
                }
            }

            sort(QpadY.begin(), QpadY.end(), sortPairsmwpc2);
            qmy2 = QpadY[0].first;
            padmy2 = QpadY[0].second;
            /*for (Int_t i =0; i<ny; i++){
                cout << "2  i = " << i << ", q = " << QpadY[i].first << ", pad = " << QpadY[i].second << endl;
            }
            cout << "qmy2 = " << qmy2 << ", padmy2 = " << padmy2 << endl;*/

            if (padmy2 + 1 < Mw2PadsY)
            {
                // Obtain position Y for 2nd charge ----
                qdown2 = fy[padmy2 - 1];
                qup2 = fy[padmy2 + 1];
                if (qmy2 > 10)
                { // && qdown2 > 0 && qup2 > 0){
                    y2 = GetPositionY(qmy2, padmy2, qdown2, qup2);
                }
                else
                    y2 = -1000.;
            }
        } // if y

        // cout << "x1 = " << x1 << ", y1 = " << y1 << " ,planex1 = " << planex1 << ", x2 = " << x2 << ", y2 = " << y2
        // << " ,planex2 = " << planex2 << endl;
        AddHitData(x1, y1, planex1);
        AddHitData(x2, y2, planex2);
    } // if all

    if (calData)
        delete[] calData;
    return;
}

// -----   Public method Reset   ------------------------------------------------
void R3BMwpc2Cal2Hit::Reset()
{
    R3BLOG(debug1, "Clearing Mwpc2HitData Structure");
    if (fMwpcHitDataCA)
        fMwpcHitDataCA->Clear();
}

// -----   Private method AddHitData  --------------------------------------------
R3BMwpcHitData* R3BMwpc2Cal2Hit::AddHitData(Double_t x, Double_t y, Int_t plane)
{
    // It fills the R3BMwpcHitData
    TClonesArray& clref = *fMwpcHitDataCA;
    Int_t size = clref.GetEntriesFast();
    return new (clref[size]) R3BMwpcHitData(x, y, plane);
}

ClassImp(R3BMwpc2Cal2Hit);
