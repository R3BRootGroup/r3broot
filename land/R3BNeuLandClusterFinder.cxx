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

// -----------------------------------------------------------------------------
// -----                   R3BNeuLandClusterFinder                         -----
// -----                Created 21-02-2012 by D.Kresan                     -----
// -----                The algorithm written by M.Heil                    -----
// -----------------------------------------------------------------------------

#include <algorithm>
#include <iostream>

#include "TClonesArray.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TMath.h"

#include "FairLogger.h"
#include "FairRootManager.h"

#include "R3BLandDigi.h"
#include "R3BNeuLandCluster.h"
#include "R3BNeuLandClusterFinder.h"

using std::cout;
using std::endl;

// Auxilliary function for sorting digis in time
// to be used with std::sort
bool SortDigiByTime(const R3BLandDigi* digi1, const R3BLandDigi* digi2);

// -----------------------------------------------------------------------------
R3BNeuLandClusterFinder::R3BNeuLandClusterFinder()
    : FairTask("R3BNeuLandClusterFinder")
{
    // Default constructor

    // Set default verbosity level
    fVerbose = 1;

    // Create array for storing the clusters
    fArrayCluster = new TClonesArray("R3BNeuLandCluster");
}
// -----------------------------------------------------------------------------

// -----------------------------------------------------------------------------
R3BNeuLandClusterFinder::R3BNeuLandClusterFinder(Int_t verbose)
    : FairTask("R3BNeuLandClusterFinder", verbose)
{
    // Standard constructor

    // Set verbosity level
    fVerbose = verbose;

    // Create array for storing the clusters
    fArrayCluster = new TClonesArray("R3BNeuLandCluster");
}
// -----------------------------------------------------------------------------

// -----------------------------------------------------------------------------
R3BNeuLandClusterFinder::~R3BNeuLandClusterFinder()
{
    // Destructor

    // Delete the array of clusters (first clear the array)
    if (fArrayCluster)
    {
        fArrayCluster->Clear();
        delete fArrayCluster;
    }
}
// -----------------------------------------------------------------------------

// -----------------------------------------------------------------------------
InitStatus R3BNeuLandClusterFinder::Init()
{
    // Task initialization

    // Get pointer to I/O Manager
    FairRootManager* rootMgr = FairRootManager::Instance();
    // Check pointer
    if (NULL == rootMgr)
    {
        cout << "-E- R3BNeuLandClusterFinder::Init : illegal pointer to root manager." << endl;
        // Stop here
        return kFATAL;
    }

    // Get pointer to LAND digi array
    fArrayDigi = dynamic_cast<TClonesArray*>(rootMgr->GetObject("LandDigi"));
    // Check pointer
    if (NULL == fArrayDigi)
    {
        cout << "-E- R3BNeuLandClusterFinder::Init : no LandDigi array." << endl;
        // Stop here
        return kFATAL;
    }

    // Register the output array (will be written to the output file)
    rootMgr->Register("NeuLandCluster", "Digi clusters in NeuLAND", fArrayCluster, kTRUE);

    // Initialize control histograms
    CreateHistograms();

    return kSUCCESS;
}
// -----------------------------------------------------------------------------

// -----------------------------------------------------------------------------
void R3BNeuLandClusterFinder::Exec(Option_t* option)
{
    // Reset the previous event
    Reset();

    // Prepare digis for sorting (fill the stl vector)
    Int_t nDigis = fArrayDigi->GetEntries();
    R3BLandDigi* digi1;
    for (Int_t i = 0; i < nDigis; i++)
    {
        digi1 = dynamic_cast<R3BLandDigi*>(fArrayDigi->At(i));
        fVectorDigi.push_back(digi1);
    } // loop over entries

    // Sort digis in time
    std::sort(fVectorDigi.begin(), fVectorDigi.end(), SortDigiByTime);
    nDigis = fVectorDigi.size();

    // find clusters and mark the position of the cluster by
    // time of first hit, position of first hit, and total energy
    Double_t /*distance,*/ posx, posy, posz, post, delt, delx, dely, delz;
    Int_t nClusters = 0;
    fhHits->Fill(nDigis);

    // Map from digi index to cluster index
    Int_t belongsToCluster[3000];
    for (Int_t i = 0; i < 3000; i++)
    {
        belongsToCluster[i] = -1;
    }

    // Declare variables outside of loop
    R3BLandDigi* digi2;
    R3BNeuLandCluster* cluster;
    Int_t oldPaddle;
    Int_t clusNo;

    // Check number of digis
    if (nDigis > 0 && nDigis <= 3000)
    {

        // Loop over all sorted digis
        for (Int_t i = 0; i < nDigis; i++)
        {
            // Get pointer to 1-st digi
            digi1 = fVectorDigi.at(i);
            // Store position, time and paddle number
            posx = digi1->GetXX();
            posy = digi1->GetYY();
            posz = digi1->GetZZ();
            post = digi1->GetTdc();
            oldPaddle = (Int_t)(digi1->GetPaddleNr() - 1);

            // find neighbour
            for (Int_t k = i; k < nDigis; k++)
            {
                // Get pointer to 2-nd digi
                digi2 = fVectorDigi.at(k);
                // Calculate differences in space-time coordinates
                delx = digi2->GetXX() - posx;
                dely = digi2->GetYY() - posy;
                delz = digi2->GetZZ() - posz;
                delt = digi2->GetTdc() - post;

                // Check if properly sorted
                if (delt < 0.)
                {
                    // 2-nd digi has earlier time. Stop here.
                    LOG(fatal) << "Exec(): Sorting of digis failed.";
                }

                // Check if neighbour
                if (TMath::Abs(delx) < 7.5 && TMath::Abs(dely) < 7.5 && TMath::Abs(delz) < 7.5 && delt >= 0.0 &&
                    delt < 1.0)
                {
                    // This is a neighbour

                    // Check if this cluster already exists
                    if (belongsToCluster[i] > -1)
                    {
                        if (i != k)
                        {
                            // Inside cluster
                            clusNo = belongsToCluster[i];
                            belongsToCluster[k] = nClusters - 1;
                            // Get pointer to the current cluster
                            cluster = dynamic_cast<R3BNeuLandCluster*>(fArrayCluster->At(clusNo));
                            // Update information of the cluster
                            cluster->SetE(cluster->GetE() + digi2->GetQdc());
                            cluster->SetSize(cluster->GetSize() + 1);
                            if (digi2->GetTdc() > cluster->GetStopT())
                            {
                                cluster->SetStopT(digi2->GetTdc());
                                cluster->SetStopPosition(digi2->GetXX(), digi2->GetYY(), digi2->GetZZ());
                            }
                        }
                    }
                    else
                    {
                        // Create new cluster
                        new ((*fArrayCluster)[nClusters]) R3BNeuLandCluster(digi2->GetXX(),
                                                                            digi2->GetYY(),
                                                                            digi2->GetZZ(),
                                                                            digi2->GetTdc(),
                                                                            digi2->GetXX(),
                                                                            digi2->GetYY(),
                                                                            digi2->GetZZ(),
                                                                            digi2->GetTdc(),
                                                                            digi2->GetQdc(),
                                                                            1);
                        belongsToCluster[k] = nClusters;
                        // Increment number of clusters
                        nClusters += 1;
                    }
                }
            }
        }

        // Fill control histograms
        Int_t totalSize = 0;
        for (Int_t i = 0; i < nClusters; i++)
        {
            cluster = dynamic_cast<R3BNeuLandCluster*>(fArrayCluster->At(i));
            fhClusterSize->Fill(cluster->GetSize());
            fhClusterEnergy->Fill(cluster->GetE());
            totalSize += cluster->GetSize();
        }
        fhClusters->Fill(nClusters);
        fhClusterNo_vs_Size->Fill(totalSize, nClusters);
    }

    // Control output
    if (fVerbose)
    {
        cout << "-I- R3BNeuLandClusterFinder : found " << nClusters << " clusters." << endl;
    }
}
// -----------------------------------------------------------------------------

// -----------------------------------------------------------------------------
void R3BNeuLandClusterFinder::CreateHistograms()
{
    // Initialize control histograms

    fhClusterSize = new TH1F("hClusterSize", "number of paddles in a cluster", 500, -0.5, 499.5);
    fhClusterEnergy = new TH1F("hClusterEnergy", "Cluster energy", 5000, 0., 500.);
    fhHits = new TH1F("hHits", "Number of hits in one event", 500, -0.5, 499.5);
    fhClusters = new TH1F("hClusters", "Number of clusters in one event", 100, -0.5, 99.5);
    fhClusterNo_vs_Size =
        new TH2F("hClusterNo_vs_Size", "Number of Clusters vs. Total Size", 100, -0.5, 99.5, 100, -0.5, 99.5);
}
// -----------------------------------------------------------------------------

// -----------------------------------------------------------------------------
void R3BNeuLandClusterFinder::Finish()
{
    // Finish task execution

    // Write histograms to the output file
    fhClusterSize->Write();
    fhClusterEnergy->Write();
    fhHits->Write();
    fhClusters->Write();
    fhClusterNo_vs_Size->Write();
}
// -----------------------------------------------------------------------------

// -----------------------------------------------------------------------------
void R3BNeuLandClusterFinder::Reset()
{
    // Reset an event

    // Clear sorted vector
    fVectorDigi.clear();

    // Clear output array
    if (fArrayCluster)
    {
        fArrayCluster->Clear();
    }
}
// -----------------------------------------------------------------------------

// -----------------------------------------------------------------------------
bool SortDigiByTime(const R3BLandDigi* digi1, const R3BLandDigi* digi2)
{
    // Auxilliary function for sorting digis in time.
    // To be used with std::sort() from <algorithm>

    if (digi1->GetTdc() < digi2->GetTdc())
    {
        return true;
    }
    else
    {
        return false;
    }
}
// -----------------------------------------------------------------------------

ClassImp(R3BNeuLandClusterFinder);
