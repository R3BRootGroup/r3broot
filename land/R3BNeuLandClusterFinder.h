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
// -----                    R3BNeuLandClusterFinder                        -----
// -----                Created 21-02-2012 by D.Kresan                     -----
// -----------------------------------------------------------------------------

#ifndef R3BNEULANDCLUSTERFINDER_H
#define R3BNEULANDCLUSTERFINDER_H

#include <vector>

#include "FairTask.h"

class TClonesArray;
class TH1F;
class TH2F;
class R3BLandDigi;

class R3BNeuLandClusterFinder : public FairTask
{

  public:
    // Default constructor
    R3BNeuLandClusterFinder();

    // Standard constructor
    R3BNeuLandClusterFinder(Int_t verbose);

    // Destructor
    virtual ~R3BNeuLandClusterFinder();

    // Task initialization
    virtual InitStatus Init();

    // Process an event
    virtual void Exec(Option_t* option);

    // Finish of task
    virtual void Finish();

    // Reset the containers
    void Reset();

  private:
    TClonesArray* fArrayDigi;              // Array of digis - input
    TClonesArray* fArrayCluster;           // Array of clusters - output
    std::vector<R3BLandDigi*> fVectorDigi; // Vector of digis (for sorting in time)

    // Control histograms
    TH1F* fhClusterSize;
    TH1F* fhClusterEnergy;
    TH1F* fhClusters;
    TH1F* fhHits;
    TH2F* fhClusterNo_vs_Size;

    // Initialization of control histograms
    void CreateHistograms();

    ClassDef(R3BNeuLandClusterFinder, 1);
};

#endif
