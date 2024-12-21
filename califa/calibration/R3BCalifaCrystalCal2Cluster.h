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

#pragma once

#ifndef R3BCALIFACRYSTALCAL2CLUSTER_H
#define R3BCALIFACRYSTALCAL2CLUSTER_H 1

#include <FairTask.h>
#include <R3BCalifaClusterData.h>
#include <R3BCalifaGeometry.h>
#include <Rtypes.h>

#include <TH2F.h>
#include <TVector3.h>

class TClonesArray;
class R3BTGeoPar;
class R3BCalifaMappingPar;

class R3BCalifaCrystalCal2Cluster : public FairTask
{
  public:
    /** Default constructor
     **/
    R3BCalifaCrystalCal2Cluster();

    /** Destructor **/
    virtual ~R3BCalifaCrystalCal2Cluster();

    /** Virtual method Exec **/
    void Exec(Option_t* /*opt*/) override;

    /** Virtual method Reset **/
    virtual void Reset();

    /** Public method SelectGeometryVersion **/
    void SelectGeometryVersion(Int_t version);

    void SetCrystalThreshold(Double_t cryThresh) { fCrystalThreshold = cryThresh; }

    void SetGammaClusterThreshold(Double_t clusterThresh) { fGammaClusterThreshold = clusterThresh; }

    void SetProtonClusterThreshold(Double_t clusterThresh) { fProtonClusterThreshold = clusterThresh; }

    /** Virtual method SetParContainers **/
    void SetParContainers() override;

    /** Accessor to select online mode **/
    void SetOnline(Bool_t option) { fOnline = option; }

    /** Virtual method Init **/
    InitStatus Init() override;

    /** Virtual method ReInit **/
    InitStatus ReInit() override;

    void SetRandomization(Bool_t rand) { fRand = rand; }

    void SetRandomizationFile(TString file)
    {
        fRandFile = file;
        fHistoFile = new TFile(fRandFile);
    }

    void IsSimulation(Bool_t simu) { fSimulation = simu; }

    void SetTotalCrystals(Int_t total) { fTotalCrystals = total; }

    bool InsideClusterWindow(TVector3 mother, TVector3 crystal);

    void SetRectangularWindow(Float_t thetaL, Float_t phiL);

    void SetRoundWindow(Double_t window);

  private:
    TClonesArray* fCrystalCalData;
    TClonesArray* fCalifaClusterData;

    TVector3 fTargetPos = { 0, 0, 0 };
    TVector3 fCalifaPos = { 0, 0, 0 };

    R3BTGeoPar* fTargetGeoPar;
    R3BTGeoPar* fCalifaGeoPar;

    Bool_t fOnline;         // Selector for online data storage
    Int_t fGeometryVersion; // Selecting the geometry of the CALIFA calorimeter
    Int_t fTotalCrystals;

    Double_t fCrystalThreshold; // Minimum energy requested in a crystal to be included in a cluster
    Double_t fProtonClusterThreshold =
        50.;                         // Minimum energy in a crystal to be considered as a proton cluster candidate
    Double_t fGammaClusterThreshold; // Minimum energy in a crystal to be considered as a gamma cluster candidate
    Double_t fProtonThreshold;       // Defines the cut energy between proton and gamma readout

    Double_t fRoundWindow; // Cluster window
    Bool_t fSimulation;    // Simulation flag

    // R3BCalifaGeometry* fCalifaGeo = nullptr;
    Bool_t fRand;      // Flag to set randomization procedure
    TString fRandFile; // File with angular coverages for each crystal
    TFile* fHistoFile = NULL;
    TH2F** fAngularDistributions;
    TString fWindowAlg;
    Float_t fThetaLimit;
    Float_t fPhiLimit;
    /** Private method AddCluster
    **
    ** Adds a CalifaCluster to the ClusterCollection
    **/
    R3BCalifaClusterData* AddCluster(std::vector<uint16_t> crystalList,
                                     Double_t ene,
                                     Double_t Nf,
                                     Double_t Ns,
                                     Double_t pAngle,
                                     Double_t aAngle,
                                     ULong64_t time,
                                     uint8_t clusterType);

  public:
    ClassDefOverride(R3BCalifaCrystalCal2Cluster, 3);
};

#endif /* R3BCalifaCrystalCal2Cluster_H */
