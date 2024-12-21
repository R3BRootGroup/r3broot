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

// ------------------------------------------------------------------------
// -----            R3BFootStripCal2Hit source file                   -----
// -----       Created 05/11/21 by J.L. Rodriguez-Sanchez             -----
// ------------------------------------------------------------------------

// ROOT headers
#include <TClonesArray.h>
#include <TF1.h>
#include <TH1F.h>
#include <TMath.h>
#include <TRandom3.h>
#include <TSpectrum.h>
#include <iomanip>
#include <vector>

// FAIR headers
#include <FairLogger.h>
#include <FairRootManager.h>
#include <FairRunAna.h>
#include <FairRuntimeDb.h>

// FOOT headers
#include "R3BFootCalData.h"
#include "R3BFootHitData.h"
#include "R3BFootMappingPar.h"
#include "R3BFootStripCal2Hit.h"
#include "R3BLogger.h"

// R3BFootStripCal2Hit::Default Constructor -------------------------------------
R3BFootStripCal2Hit::R3BFootStripCal2Hit()
    : R3BFootStripCal2Hit("R3BFootStripCal2Hit", 1)
{
}

// R3BFootStripCal2HitPar::Standard Constructor ---------------------------------
R3BFootStripCal2Hit::R3BFootStripCal2Hit(const TString& name, Int_t iVerbose)
    : FairTask(name, iVerbose)
{
}

// Virtual R3BFootStripCal2Hit::Destructor --------------------------------------
R3BFootStripCal2Hit::~R3BFootStripCal2Hit()
{
    R3BLOG(debug1, "");
    if (fFootHitData)
        delete fFootHitData;
}

void R3BFootStripCal2Hit::SetParContainers()
{
    // Parameter Container
    // Reading footMappingPar from FairRuntimeDb
    FairRuntimeDb* rtdb = FairRuntimeDb::instance();
    R3BLOG_IF(fatal, !rtdb, "FairRuntimeDb not found");

    fMap_Par = dynamic_cast<R3BFootMappingPar*>(rtdb->getContainer("footMappingPar"));
    if (!fMap_Par)
    {
        R3BLOG(error, "Couldn't get handle on footMappingPar container");
    }
    else
    {
        R3BLOG(info, "footMappingPar found");
    }
}

void R3BFootStripCal2Hit::SetParameter()
{
    if (!fMap_Par)
    {
        R3BLOG(warn, "Container footMappingPar not found");
        return;
    }
    //--- Parameter Container ---
    fMaxNumDet = fMap_Par->GetNumDets(); // Number of foot detectors

    R3BLOG(info, "NumDet from mapping " << fMaxNumDet);
    for (int i = 0; i < fMaxNumDet; i++)
    {
        fDistTarget.push_back(fMap_Par->GetDist2target(i + 1));
        fAngleTheta.push_back(fMap_Par->GetAngleTheta(i + 1));
        fAnglePhi.push_back(fMap_Par->GetAnglePhi(i + 1));
        fOffsetX.push_back(fMap_Par->GetOffsetX(i + 1));
        fOffsetY.push_back(fMap_Par->GetOffsetY(i + 1));
    }
    fMap_Par->printParams();
}

// -----   Public method Init   -------------------------------------------------
InitStatus R3BFootStripCal2Hit::Init()
{
    R3BLOG(info, "");

    FairRootManager* rootManager = FairRootManager::Instance();
    R3BLOG_IF(fatal, rootManager == nullptr, "FairRootManager not found");

    // INPUT DATA
    fFootCalData = dynamic_cast<TClonesArray*>(rootManager->GetObject("FootCalData"));
    if (!fFootCalData)
    {
        R3BLOG(fatal, "FootCalData not found");
        return kFATAL;
    }

    // Set container with mapping parameters
    SetParameter();

    // Output data
    fFootHitData = new TClonesArray("R3BFootHitData");
    rootManager->Register("FootHitData", "FOOT Hit", fFootHitData, !fOnline);
    fFootHitData->Clear();

    char Name[255];
    for (Int_t i = 0; i < fMaxNumDet; i++)
    {
        sprintf(Name, "hssd_%d", i + 1);
        hssd[i] = new TH1F(Name, "", 640, -0.5, 640.5);
    }

    return kSUCCESS;
}

// -----   Public method ReInit   -----------------------------------------------
InitStatus R3BFootStripCal2Hit::ReInit()
{
    SetParContainers();
    SetParameter();
    return kSUCCESS;
}

// -----   Public method Execution   --------------------------------------------
void R3BFootStripCal2Hit::Exec(Option_t* /*option*/)
{
    // Reset entries in output arrays, local arrays
    Reset();

    // Reading the Input Cal Data
    int nHits = fFootCalData->GetEntriesFast();
    if (nHits == 0)
        return;

    // Data from cal level
    R3BFootCalData** calData = new R3BFootCalData*[nHits];
    int detId;
    int stripId;
    double energy;
    double x = 0., y = 0., z = 0.;
    std::vector<std::vector<int>> StripI;
    std::vector<std::vector<double>> StripE;

    // Clustering algorithm - A. Revel
    StripI.resize(fMaxNumDet);
    StripE.resize(fMaxNumDet);

    ClusterPos = std::vector<std::vector<double>>(fMaxNumDet, std::vector<double>(200, std::nan("")));

    ClusterESum = std::vector<std::vector<double>>(fMaxNumDet, std::vector<double>(200, std::nan("")));

    ClusterNStrip = std::vector<std::vector<int>>(fMaxNumDet, std::vector<int>(200, 0));

    ClusterI = std::vector<std::vector<std::vector<double>>>(
        fMaxNumDet, std::vector<std::vector<double>>(200, std::vector<double>(0)));

    ClusterE = std::vector<std::vector<std::vector<double>>>(
        fMaxNumDet, std::vector<std::vector<double>>(200, std::vector<double>(0)));

    // Filling vectors
    for (int i = 0; i < nHits; i++)
    {
        calData[i] = (R3BFootCalData*)(fFootCalData->At(i));
        detId = calData[i]->GetDetId() - 1;
        stripId = calData[i]->GetStripId() - 1;
        energy = calData[i]->GetEnergy();

        StripI[detId].push_back(stripId);
        StripE[detId].push_back(energy);
    }
    if (calData)
        delete[] calData;

    // Sort (should be good by default but just in case)
    for (int i = 0; i < fMaxNumDet; i++)
    {
        if (StripI[i].size() > 1)
        {
            for (int j = 0; j < StripI[i].size() - 1; j++)
            {
                for (int k = j + 1; k < StripI[i].size(); k++)
                {
                    if (StripI[i][j] > StripI[i][k])
                    {
                        auto tempI = StripI[i][j];
                        StripI[i][j] = StripI[i][k];
                        StripI[i][k] = tempI;
                        auto tempD = StripE[i][j];
                        StripE[i][j] = StripE[i][k];
                        StripE[i][k] = tempD;
                    }
                }
            }
        }
    }

    // Clustering
    for (int i = 0; i < fMaxNumDet; i++)
    {

        int ClusterCount = 0;
        auto TempI = 0;

        int j = 0;

        while (j < StripI[i].size())
        {

            if (j == 0)
            {
                TempI = StripI[i][j];
            }

            ClusterNStrip[i][ClusterCount]++;
            ClusterI[i][ClusterCount].push_back(StripI[i][j]);
            ClusterE[i][ClusterCount].push_back(StripE[i][j]);

            int k = j + 1;

            while (k < StripI[i].size() && (StripI[i][k] - TempI) == 1)
            {
                ClusterNStrip[i][ClusterCount]++;
                ClusterI[i][ClusterCount].push_back(StripI[i][k]);
                ClusterE[i][ClusterCount].push_back(StripE[i][k]);
                TempI = StripI[i][k];
                StripI[i].erase(StripI[i].begin() + k);
                StripE[i].erase(StripE[i].begin() + k);
            }

            StripI[i].erase(StripI[i].begin() + j);
            StripE[i].erase(StripE[i].begin() + j);

            ClusterCount++;
        }
        ClusterMult[i] = ClusterCount;
    }

    // Compute Sum Energy, Position and Eta
    for (int i = 0; i < fMaxNumDet; i++)
    {
        for (int j = 0; j < ClusterMult[i]; j++)
        {
            for (int k = 0; k < ClusterNStrip[i][j]; k++)
            {
                ClusterESum[i][j] += ClusterE[i][j][k];
                ClusterPos[i][j] += ClusterE[i][j][k] * ClusterI[i][j][k];
            }
            ClusterPos[i][j] = ClusterPos[i][j] / ClusterESum[i][j];
            // Eta[i][j] = ClusterPos[i][j] - (double)ClusterPos[i][j];
        }
    }

    // Sort Cluster from Higher to Lower Energy
    for (int i = 0; i < fMaxNumDet; i++)
    {
        for (int j = 0; j < ClusterMult[i] - 1; j++)
        {
            for (int k = j + 1; k < ClusterMult[i]; k++)
            {
                if (ClusterESum[i][j] < ClusterESum[i][k])
                {
                    auto tempI = ClusterNStrip[i][j];
                    ClusterNStrip[i][j] = ClusterNStrip[i][k];
                    ClusterNStrip[i][k] = tempI;

                    double tempD = ClusterESum[i][j];
                    ClusterESum[i][j] = ClusterESum[i][k];
                    ClusterESum[i][k] = tempD;

                    tempD = ClusterPos[i][j];
                    ClusterPos[i][j] = ClusterPos[i][k];
                    ClusterPos[i][k] = tempD;

                    // tempD = Eta[i][j];
                    // Eta[i][j] = Eta[i][k];
                    // Eta[i][k] = tempD;
                }
            }
        }
    }

    // Filling HitData
    for (uint8_t i = 0; i < fMaxNumDet; i++)
    {
        for (int j = 0; j < ClusterMult[i]; j++)
        {
            double pos = 100. * ClusterPos[i][j] / 640. - fMiddle;

            if (fAnglePhi[i] == 0.)
            { // X-Foot (StripId numbered from left to right)
                x = pos * TMath::Cos(fAngleTheta[i] * TMath::DegToRad()) + fOffsetX[i];
                y = fOffsetY[i];
                z = pos * TMath::Sin(fAngleTheta[i] * TMath::DegToRad()) + fDistTarget[i];
            }
            else if (fAnglePhi[i] == 90.)
            { // Y-Foot (StripId numbered from bottom to top)
                x = fOffsetX[i];
                y = pos + fOffsetY[i];
                z = fDistTarget[i];
            }
            else if (fAnglePhi[i] == 180.)
            { // X-Foot (StripId numbered from right to left)
                x = -pos * TMath::Cos(fAngleTheta[i] * TMath::DegToRad()) + fOffsetX[i];
                y = fOffsetY[i];
                z = -pos * TMath::Sin(fAngleTheta[i] * TMath::DegToRad()) + fDistTarget[i];
            }
            else if (fAnglePhi[i] == 270.)
            { // Y-Foot (StripId numbered from top to bottom)
                x = fOffsetX[i];
                y = -pos + fOffsetY[i];
                z = fDistTarget[i];
            }
            else
            {
                LOG(info) << "R3BFootStripCal2Hit::AnglePhi is Wrong !";
            }

            TVector3 master(x, y, z);
            // TODO: Eta correction not used at the moment

            if (ClusterESum[i][j] > fThSum && j < fMaxNumClusters)
            {
                AddHitData(i + 1, ClusterNStrip[i][j], pos, master, ClusterESum[i][j], ClusterMult[i]);
            }
        }
    }

    for (int i = 0; i < fMaxNumDet; i++)
        hssd[i]->Reset();

    StripI.clear();
    StripE.clear();
    ClusterMult.clear();
    ClusterPos.clear();
    ClusterESum.clear();
    ClusterNStrip.clear();
    // Eta.clear();
    ClusterI.clear();
    ClusterE.clear();

    return;
}

// -----   Public method Reset   ------------------------------------------------
void R3BFootStripCal2Hit::Reset()
{
    R3BLOG(debug, "Clearing HitData Structure");
    if (fFootHitData)
        fFootHitData->Clear();
}

// -----   Private method AddHitData  -------------------------------------------
R3BFootHitData* R3BFootStripCal2Hit::AddHitData(uint8_t detid,
                                                uint16_t numhit,
                                                double pos,
                                                TVector3 master,
                                                double energy_s,
                                                uint16_t mulS)
{
    // It fills the R3BFootHitData
    TClonesArray& clref = *fFootHitData;
    int size = clref.GetEntriesFast();
    return new (clref[size]) R3BFootHitData(detid, numhit, pos, master, energy_s, mulS);
}

ClassImp(R3BFootStripCal2Hit)
