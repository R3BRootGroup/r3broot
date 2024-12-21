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
// -----                        R3BMwpc3Cal2Hit                     -----
// -----             Created 14/10/19 by G. García Jiménez          -----
// -----             by modifying J.L classes for MWPC0             -----
// -----               s455 method 17/06/22 by Antia GG             -----
// ----------------------------------------------------------------------

#include "TClonesArray.h"
#include "TF1.h"
#include "TGraph.h"
#include "TMath.h"

// Fair headers
#include "FairLogger.h"
#include "FairRootManager.h"
#include "FairRunAna.h"
#include "FairRuntimeDb.h"

// MWPC headers
#include "R3BEventHeader.h"
#include "R3BLogger.h"
#include "R3BMwpc3Cal2Hit.h"
#include "R3BMwpcCalData.h"
#include "R3BMwpcHitData.h"

bool compare(pair<int, int> p1, pair<int, int> p2) { return p1.second < p2.second; }

/* ---- R3BMwpc3Cal2Hit: Default Constructor ---- */

R3BMwpc3Cal2Hit::R3BMwpc3Cal2Hit()
    : R3BMwpc3Cal2Hit("R3BMwpc3Cal2Hit", 1)
{
}

/* ---- R3BMwpc3Cal2Hit: Standard Constructor ---- */
R3BMwpc3Cal2Hit::R3BMwpc3Cal2Hit(const char* name, Int_t iVerbose)
    : FairTask(name, iVerbose)
    , fMwpcCalDataCA(NULL)
    , fMwpcHitDataCA(NULL)
    , fTofWallHitDataCA(NULL)
    , fwy(5.000)
    , fwx(3.125)
    , fSizeX(900.0)
    , fSizeY(600.0) // in mm
    , fOnline(kFALSE)
    , fTofWallMatching(kFALSE)
    , fExpId(0)
{
}

/* ---- Virtual R3BMwpc3Cal2Hit: Destructor ---- */
R3BMwpc3Cal2Hit::~R3BMwpc3Cal2Hit()
{
    R3BLOG(debug1, "Destructor");
    if (fMwpcCalDataCA)
        delete fMwpcCalDataCA;
    if (fMwpcHitDataCA)
        delete fMwpcHitDataCA;
    if (fTofWallHitDataCA)
        delete fTofWallHitDataCA;
}

/* ---- Public method Init   ---- */
InitStatus R3BMwpc3Cal2Hit::Init()
{
    R3BLOG(info, "");
    // INPUT DATA
    FairRootManager* rootManager = FairRootManager::Instance();
    if (!rootManager)
    {
        R3BLOG(fatal, "FairRootManager not found");
        return kFATAL;
    }

    header = (R3BEventHeader*)rootManager->GetObject("EventHeader.");
    if (!header)
        header = (R3BEventHeader*)rootManager->GetObject("R3BEventHeader");

    fMwpcCalDataCA = (TClonesArray*)rootManager->GetObject("Mwpc3CalData");
    if (!fMwpcCalDataCA)
    {
        R3BLOG(fatal, "Mwpc1CalData not found");
        return kFATAL;
    }

    fTofWallHitDataCA = (TClonesArray*)rootManager->GetObject("TofWHitData");
    if (!fTofWallHitDataCA)
    {
        fTofWallMatching = kFALSE;
    }
    else
        fTofWallMatching = kTRUE;

    // OUTPUT DATA
    // Hit data
    fMwpcHitDataCA = new TClonesArray("R3BMwpcHitData");
    rootManager->Register("Mwpc3HitData", "MWPC3 Hit", fMwpcHitDataCA, !fOnline);
    fMwpcHitDataCA->Clear();
    return kSUCCESS;
}

/* ----   Public method ReInit   ---- */
InitStatus R3BMwpc3Cal2Hit::ReInit() { return kSUCCESS; }

bool sortPairs(const pair<Double_t, Int_t>& x, const pair<Double_t, Int_t>& y) { return x.first > y.first; }

/* ----   Public method Execution   ---- */
void R3BMwpc3Cal2Hit::Exec(Option_t* option)
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
void R3BMwpc3Cal2Hit::S467()
{
    // Reset entries in output arrays, local arrays
    Reset();

    if (fTofWallMatching)
        ReconstructHitWithTofWallMatching();

    else
    {
        // Reading the Input -- Cal Data --
        Int_t nHits = fMwpcCalDataCA->GetEntriesFast();
        if (nHits == 0)
            return;

        // Data from cal level
        R3BMwpcCalData** calData = new R3BMwpcCalData*[nHits];
        Int_t planeId = 0;
        Int_t padId = 0;
        Int_t padmx = -1, padmy = -1;
        Double_t q = 0., qmx = 0., qmy = 0., qleft = 0., qright = 0., qdown = 0., qup = 0.;
        Double_t x = NAN, y = NAN;

        vector<double> vQX, vQY;
        vector<int> vStripX, vStripY;
        vQX.clear();
        vQY.clear();
        vStripX.clear();
        vStripY.clear();

        for (Int_t i = 0; i < Mw3PadsX; i++)
            fx[i] = 0;
        for (Int_t i = 0; i < Mw3PadsY; i++)
            fy[i] = 0;

        for (Int_t i = 0; i < nHits; i++)
        {
            calData[i] = (R3BMwpcCalData*)(fMwpcCalDataCA->At(i));
            planeId = calData[i]->GetPlane();
            padId = calData[i]->GetPad() - 1;
            q = calData[i]->GetQ();
            if (planeId == 1)
            {
                fx[padId] = q;
                vQX.push_back(q);
                vStripX.push_back(padId);
            }
            else
            {
                fy[padId] = q;
                vQY.push_back(q);
                vStripY.push_back(padId);
            }

            if (q > qmx && planeId == 1)
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
        if (padmx > 0 && padmy > 0 && padmx + 1 < Mw3PadsX && padmy + 1 < Mw3PadsY && qmx > 0 && qmy > 0)
        {
            // Obtain position X ----
            qleft = (Double_t)fx[padmx - 1];
            qright = (Double_t)fx[padmx + 1];
            // std::cout<<qleft<<" "<<qright<<std::endl;
            if (qleft > 0 && qright > 0)
                x = GetPositionX(qmx, padmx, qleft, qright);
            // x = FittedHyperbolicSecant("X",vQX,vStripX,qmx,padmx);

            // Obtain position Y ----
            qdown = fy[padmy - 1];
            qup = fy[padmy + 1];
            // if(padmy==64) std::cout << padmy << " " << qmy << " " << qdown << " " << qup << std::endl;
            // if(padmy==63) std::cout << padmy << " " << qmy << " " << qdown << " " << qup << std::endl;
            if (qdown > 0 && qup > 0)
                y = GetPositionY(qmy, padmy, qdown, qup);
            // y = FittedHyperbolicSecant("Y",vQY,vStripY,qmy,padmy);

            AddHitData(x, y);
        }

        if (calData)
            delete[] calData;
    }
    return;
}

// -----   Private method s455       --------------------------------------------
void R3BMwpc3Cal2Hit::S455()
{
    // Reset entries in output arrays, local arrays
    Reset();

    if (fTofWallMatching)
        ReconstructHitWithTofWallMatching();

    else
    {
        // Reading the Input -- Cal Data --
        Int_t nHits = fMwpcCalDataCA->GetEntriesFast();
        R3BMwpcCalData** calData = new R3BMwpcCalData*[nHits];

        Int_t planeId = 0;
        Int_t padId = 0;
        Double_t q = 0;
        Int_t padmx1 = -1, padmy1 = -1, padmx2 = -1, padmy2 = -1;
        Double_t qmx1 = 0., qmy1 = 0., qleft1 = 0., qright1 = 0., qdown1 = 0., qup1 = 0.;
        Double_t qmx2 = 0., qmy2 = 0., qleft2 = 0., qright2 = 0., qdown2 = 0., qup2 = 0.;
        Double_t x1 = 0., y1 = 0., x2 = 0., y2 = 0.;
        Int_t nx = 0, ny = 0;

        vector<pair<Double_t, Int_t>> QpadX, QpadY;
        QpadX.clear();
        QpadY.clear();
        Bool_t xexists = false;
        Bool_t yexists = false;
        // Double_t fx[Mw3PadsX], fy[Mw3PadsY];

        for (Int_t i = 0; i < Mw3PadsX; i++)
            fx[i] = 0;
        for (Int_t i = 0; i < Mw3PadsY; i++)
            fy[i] = 0;
        for (Int_t i = 0; i < nHits; i++)
        {
            calData[i] = (R3BMwpcCalData*)(fMwpcCalDataCA->At(i));
            planeId = calData[i]->GetPlane();
            padId = calData[i]->GetPad() - 1;
            q = calData[i]->GetQ();
            // cout << "i = " << i << ", q = " << q << ", padId = " << padId << ", planeId = " << planeId << endl;
            pair<Double_t, Int_t> hit_pair = make_pair(q, padId);
            if (planeId == 1)
            {
                fx[padId] = q;
                QpadX.push_back(hit_pair);
                xexists = true;
                nx = nx + 1;
            }
            else
            {
                fy[padId] = q;
                QpadY.push_back(hit_pair);
                yexists = true;
                ny = ny + 1;
            }
        }

        if (xexists == true && yexists == true)
        {
            sort(QpadX.begin(), QpadX.end(), sortPairs); // el vector se ordena por el primer elemento de cada par
            qmx1 = QpadX[0].first;
            padmx1 = QpadX[0].second;
            // cout << "qmx1 = " << qmx1 << ", padmx1 = " << padmx1 << endl;

            sort(QpadY.begin(), QpadY.end(), sortPairs);
            qmy1 = QpadY[0].first;
            padmy1 = QpadY[0].second;
            // cout << "qmy1 = " << qmy1 << ", padmy1 = " << padmy1 << endl;

            if (padmx1 + 1 < Mw3PadsX && padmy1 + 1 < Mw3PadsY)
            {
                // Obtain position X for 1st charge----
                qleft1 = fx[padmx1 - 1];
                qright1 = fx[padmx1 + 1];
                if (qmx1 > 10 && qleft1 > 0 && qright1 > 0)
                {
                    x1 = GetPositionX(qmx1, padmx1, qleft1, qright1);
                }
                else
                    x1 = -1000.;

                // Obtain position Y for 1st charge ----
                qdown1 = fy[padmy1 - 1];
                qup1 = fy[padmy1 + 1];
                if (qmy1 > 10 && qdown1 > 0 && qup1 > 0)
                {
                    y1 = GetPositionY(qmy1, padmy1, qdown1, qup1);
                }
                else
                    y1 = -1000.;

                // if (x1!=0 && y1!=0){
                AddHitData(x1, y1);
                //}
            }

            for (Int_t i = 0; i < nx; i++)
            {
                if (QpadX[i].second == padmx1 - 3 || QpadX[i].second == padmx1 - 2 || QpadX[i].second == padmx1 - 1 ||
                    QpadX[i].second == padmx1 || QpadX[i].second == padmx1 + 1 || QpadX[i].second == padmx1 + 2 ||
                    QpadX[i].second == padmx1 + 3)
                {
                    QpadX[i].first = 0.;
                }
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
            sort(QpadX.begin(), QpadX.end(), sortPairs);
            qmx2 = QpadX[0].first;
            padmx2 = QpadX[0].second;
            // cout << "qmx2 = " << qmx2 << ", padmx2 = " << padmx2 << endl;

            sort(QpadY.begin(), QpadY.end(), sortPairs);
            qmy2 = QpadY[0].first;
            padmy2 = QpadY[0].second;
            // cout << "qmy2 = " << qmy2 << ", padmy2 = " << padmy2 << endl;

            if (padmx2 + 1 < Mw3PadsX && padmy2 + 1 < Mw3PadsY)
            {
                // Obtain position X for 2nd charge----
                qleft2 = fx[padmx2 - 1];
                qright2 = fx[padmx2 + 1];
                if (qmx2 > 10 && qleft2 > 0 && qright2 > 0)
                {
                    x2 = GetPositionX(qmx2, padmx2, qleft2, qright2);
                }
                else
                    x2 = -1000;

                // Obtain position Y for 2nd charge ----
                qdown2 = fy[padmy2 - 1];
                qup2 = fy[padmy2 + 1];
                if (qmy2 > 10 && qdown2 > 0 && qup2 > 0)
                {
                    y2 = GetPositionY(qmy2, padmy2, qdown2, qup2);
                }
                else
                    y2 = -1000.;
                // if (x2!=0 && y2!=0){
                AddHitData(x2, y2);
                //}
            }

        } // Condition if (xexists == true && yexists == true)

        if (calData)
            delete[] calData;
    }

    return;
}

/* ---- Fitted Hyperbolic Secant ---- */
Double_t R3BMwpc3Cal2Hit::FittedHyperbolicSecant(string XorY,
                                                 vector<double> vQ,
                                                 vector<int> vStrip,
                                                 int Qmax,
                                                 int StripMax)
{
    Double_t pos = -1000;
    static TF1* f = new TF1("sechs", "[0]/(cosh(TMath::Pi()*(x-[1])/[2])*cosh(TMath::Pi()*(x-[1])/[2]))", 1, 300);

    double StartingPoint = StripMax;

    f->SetParameters(Qmax, StartingPoint, 2.5);
    f->SetParLimits(0, 0.2 * Qmax, 1.2 * Qmax);

    static vector<double> y;
    static vector<double> q;
    y.clear();
    q.clear();
    double final_size = 0;
    unsigned int sizeQ = vQ.size();

    for (unsigned int i = 0; i < sizeQ; i++)
    {
        if (vQ[i] > Qmax * 0.05)
        {
            q.push_back(vQ[i]);
            y.push_back(vStrip[i]);
            final_size++;
        }
    }

    if (final_size < 3)
        return -1000;

    TGraph* g = new TGraph(q.size(), &y[0], &q[0]);
    g->Fit(f, "QN0");
    delete g;
    Double_t pospad = f->GetParameter(1);

    if (XorY == "Y")
        pos = pospad * fwy - fSizeY / 2 + fwy / 2;
    if (XorY == "X")
        pos = -(pospad * fwx - fSizeX / 2 + fwx / 2);

    return pos;
}

/* ---- Indentifying cluster and matching good hits using the tofwall information ---- */
void R3BMwpc3Cal2Hit::ReconstructHitWithTofWallMatching()
{
    // Getting Position information from tof wall
    Int_t twHits = fTofWallHitDataCA->GetEntriesFast();
    if (!twHits)
        return;
    R3BSofTofWHitData** twHitData;
    twHitData = new R3BSofTofWHitData*[twHits];
    vector<double> twX;
    vector<double> twY;
    twX.clear();
    twY.clear();
    for (Int_t i = 0; i < twHits; i++)
    {
        twHitData[i] = (R3BSofTofWHitData*)(fTofWallHitDataCA->At(i));
        twX.push_back(twHitData[i]->GetX());
        twY.push_back(twHitData[i]->GetY());
    }

    if (twHitData)
        delete twHitData;

    // Getting Position information from mwpc3
    Int_t nHits = fMwpcCalDataCA->GetEntriesFast();
    if (nHits == 0)
        return;
    // Data from cal level
    R3BMwpcCalData** calData;
    calData = new R3BMwpcCalData*[nHits];
    Int_t planeId;
    Int_t padId;
    Int_t q;
    /*fQX.clear();
    fQY.clear();
    fPadX.clear();
    fPadY.clear();*/
    fPairX.clear();
    fPairY.clear();
    fPairX.reserve(100);
    fPairY.reserve(100);
    fClusterX.clear();
    fClusterY.clear();

    for (Int_t i = 0; i < nHits; i++)
    {
        calData[i] = (R3BMwpcCalData*)(fMwpcCalDataCA->At(i));
        planeId = calData[i]->GetPlane();
        padId = calData[i]->GetPad() - 1;
        q = calData[i]->GetQ();
        pair<Int_t, Int_t> hit_pair = make_pair(padId, q);
        if (planeId == 1)
        {
            fPairX.push_back(hit_pair);
            // fQX.push_back(q);
            // fPadX.push_back(padId);
        }
        if (planeId == 3)
        {
            fPairY.push_back(hit_pair);
            // fQY.push_back(q);
            // fPadY.push_back(padId);
        }
    }
    if (fPairX.size() == 0 || fPairY.size() == 0)
    {
        if (calData)
            delete[] calData;
        return;
    }

    // sorting fPair with increasing pad number;
    sort(fPairX.begin(), fPairX.end());
    sort(fPairY.begin(), fPairY.end());

    fThresholdX = max(0.2 * GetChargeMax(fPairX), 500.);
    fThresholdY = max(0.2 * GetChargeMax(fPairY), 500.);

    while (GetChargeMax(fPairX) > fThresholdX && fPairX.size() > 3)
    {
        fClusterX.push_back(FindCluster(fPairX));
    }
    while (GetChargeMax(fPairY) > fThresholdY && fPairY.size() > 3)
    {
        fClusterY.push_back(FindCluster(fPairY));
    }

    if (fClusterX.size() == 0 || fClusterY.size() == 0)
    {
        if (calData)
            delete[] calData;
        return;
    }

    vector<double> Xpos;
    for (unsigned int i = 0; i < fClusterX.size(); i++)
    {
        vector<pair<int, int>> hit;
        hit = fClusterX[i];
        double x;
        int qleft = hit[0].second;
        int qmax = hit[1].second;
        int padmax = hit[1].first;
        int qright = hit[2].second;
        if (padmax > 0 && padmax + 1 < Mw3PadsX && qmax > 0 && qleft > 0 && qright > 0)
        {
            x = GetPositionX(qmax, padmax, qleft, qright);
            Xpos.push_back(x);
        }
    }

    vector<double> Ypos;
    for (unsigned int i = 0; i < fClusterY.size(); i++)
    {
        vector<pair<int, int>> hit;
        hit = fClusterY[i];
        double y;
        int qdown = hit[0].second;
        int qmax = hit[1].second;
        int padmax = hit[1].first;
        int qup = hit[2].second;
        if (padmax > 0 && padmax + 1 < Mw3PadsY && qmax > 0 && qdown > 0 && qup > 0)
        {
            y = GetPositionY(qmax, padmax, qdown, qup);
            Ypos.push_back(y);
        }
    }

    if (Xpos.size() == 0 || Ypos.size() == 0)
    {
        if (calData)
            delete[] calData;
        return;
    }

    if (twX.size() == twY.size())
    {
        TofWallMatching(twX, twY, Xpos, Ypos);
    }
    if (calData)
        delete[] calData;
    return;
}

/* ---- Matching position with the tofwall ---- */
void R3BMwpc3Cal2Hit::TofWallMatching(vector<double> twx,
                                      vector<double> twy,
                                      vector<double> mwpcx,
                                      vector<double> mwpcy)
{
    vector<double> good_posx;
    vector<double> good_posy;
    for (unsigned int i = 0; i < twx.size(); i++)
    {
        double tofx = twx[i];
        double tofy = twy[i];
        for (unsigned int k = 0; k < mwpcx.size(); k++)
        {
            double posx = mwpcx[k];
            if (abs(posx - tofx) < 12.5)
            {
                good_posx.push_back(posx);
            }
        }
        for (unsigned int p = 0; p < mwpcy.size(); p++)
        {
            double posy = mwpcy[p];
            if (abs(posy - tofy) < 10)
            {
                good_posy.push_back(posy);
            }
        }
    }
    if (good_posx.size() == 2 && good_posy.size() == 2)
    {
        for (unsigned int i = 0; i < good_posx.size(); i++)
        {
            AddHitData(good_posx[i], good_posy[i]);
        }
    }

    return;
}

/* ---- Getting a cluster of three pairs (pad,charge) around the maximum charge found ---- */
/* ---- The three points are then removed from the original fPair to find the next cluster ---- */
vector<pair<int, int>> R3BMwpc3Cal2Hit::FindCluster(vector<pair<int, int>>& p1)
{
    vector<pair<int, int>> hit;
    hit.clear();
    pair<int, int> couple;

    const auto pelt = max_element(p1.begin(), p1.end(), compare);
    int padmax = pelt->first;
    int qmax = pelt->second;
    int imax = distance(p1.begin(), pelt);

    couple = { p1[imax - 1].first, p1[imax - 1].second };
    // cout << "pad= " << p1[imax-1].first << " / Q= " << p1[imax-1].second << endl;
    hit.push_back(couple);
    couple = { p1[imax].first, p1[imax].second };
    // cout << "pad= " << p1[imax].first << " / Q= " << p1[imax].second << endl;
    hit.push_back(couple);
    couple = { p1[imax + 1].first, p1[imax + 1].second };
    // cout << "pad= " << p1[imax+1].first << " / Q= " << p1[imax+1].second << endl;
    hit.push_back(couple);

    // Removing those three points from fPair
    p1.erase(p1.begin() + imax - 1, p1.begin() + imax);
    p1.erase(p1.begin() + imax - 1, p1.begin() + imax);
    p1.erase(p1.begin() + imax - 1, p1.begin() + imax);

    return hit;
}
/* ---- Getting the maximum charger of the distribution ---- */
double R3BMwpc3Cal2Hit::GetChargeMax(vector<pair<int, int>> p1)
{
    const auto comp = max_element(p1.begin(), p1.end(), compare);
    double Qmax = comp->second;
    return Qmax;
}

/* ----   Protected method to obtain the position X ---- */
Double_t R3BMwpc3Cal2Hit::GetPositionX(Double_t qmax, Int_t padmax, Double_t qleft, Double_t qright)
{
    Double_t a3 = TMath::Pi() * fwx / (TMath::ACosH(0.5 * (TMath::Sqrt(qmax / qleft) + TMath::Sqrt(qmax / qright))));
    // Double_t a2 = gRandom->Uniform(-fwx / 2,fwx / 2);
    Double_t a2 = (a3 / TMath::Pi()) * TMath::ATanH((TMath::Sqrt(qmax / qleft) - TMath::Sqrt(qmax / qright)) /
                                                    (2 * TMath::SinH(TMath::Pi() * fwx / a3)));

    return (-1. * padmax * fwx + (fSizeX / 2) - (fwx / 2) - a2); // Left is positive and right negative
}

/* ----   Protected method to obtain the position Y ---- */
Double_t R3BMwpc3Cal2Hit::GetPositionY(Double_t qmax, Int_t padmax, Double_t qdown, Double_t qup)
{
    Double_t a3 = TMath::Pi() * fwy / (TMath::ACosH(0.5 * (TMath::Sqrt(qmax / qdown) + TMath::Sqrt(qmax / qup))));
    // Double_t a2 = gRandom->Uniform(-fwy / 2, fwy / 2);
    Double_t a2 = (a3 / TMath::Pi()) * TMath::ATanH((TMath::Sqrt(qmax / qdown) - TMath::Sqrt(qmax / qup)) /
                                                    (2 * TMath::SinH(TMath::Pi() * fwy / a3)));

    return (padmax * fwy - (fSizeY / 2) + (fwy / 2) + a2);
}

/* ----   Public method Reset  ---- */
void R3BMwpc3Cal2Hit::Reset()
{
    R3BLOG(debug1, "Clearing Mwpc3HitData Structure");
    if (fMwpcHitDataCA)
        fMwpcHitDataCA->Clear();
}

/* ----   Private method AddHitData  ---- */
R3BMwpcHitData* R3BMwpc3Cal2Hit::AddHitData(Double_t x, Double_t y)
{
    // It fills the R3BMwpcHitData
    TClonesArray& clref = *fMwpcHitDataCA;
    Int_t size = clref.GetEntriesFast();
    return new (clref[size]) R3BMwpcHitData(x, y);
}

ClassImp(R3BMwpc3Cal2Hit);
