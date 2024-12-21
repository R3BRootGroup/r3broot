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

#ifndef R3BNEULANDCLUSTERMON_H
#define R3BNEULANDCLUSTERMON_H

/** Neuland Cluster Monitor
 * @author Jan Mayer
 * @since  10.06.2015
 *
 *  Input:  Clusters. Should work both on Clusters from Monte Carlo simulations
 *          as well as experimental data on digi level.
 *  Output: Soon: Various diagrams.
 *          Currently: A 3D scatter plot for each event.
 */

#include "FairTask.h"
#include "Filterable.h"
#include "R3BNeulandCluster.h"
#include "TCAConnector.h"
#include <vector>

class TH1D;
class TH2D;
class TH3D;

class R3BNeulandClusterMon : public FairTask
{
  public:
    R3BNeulandClusterMon(TString input = "NeulandClusters",
                         TString output = "NeulandClusterMon",
                         const Option_t* option = "");

    ~R3BNeulandClusterMon() override = default;

    // No copy and no move is allowed (Rule of three/five)
    R3BNeulandClusterMon(const R3BNeulandClusterMon&) = delete;            // copy constructor
    R3BNeulandClusterMon(R3BNeulandClusterMon&&) = delete;                 // move constructor
    R3BNeulandClusterMon& operator=(const R3BNeulandClusterMon&) = delete; // copy assignment
    R3BNeulandClusterMon& operator=(R3BNeulandClusterMon&&) = delete;      // move assignment

  protected:
    InitStatus Init() override;
    void Finish() override;

  public:
    void Exec(Option_t*) override;
    void SetBeta(double b) { fBeta = b; }
    void AddFilter(const Filterable<R3BNeulandCluster*>::Filter& f) { fClusterFilters.Add(f); }

  private:
    TCAInputConnector<R3BNeulandCluster> fNeulandClusters;
    TString fOutput;

    double fBeta; // just for ClusterRValue

    Bool_t fIs3DTrackEnabled;
    TH3D* fh3;

    TH1D* fhClusters;
    TH1D* fhClusterSize;
    TH1D* fhClusterEnergy;
    TH1D* fhClusterRValue;
    TH2D* fhClusterEnergyVSSize;
    TH1D* fhClusterEToF;
    TH1D* fhClusterTime;
    TH2D* fhClusterNumberVSEnergy;
    TH2D* fhClusterEToFVSEnergy;
    TH2D* fhClusterEToFVSTime;
    TH2D* fhClusterEVSTime;

    TH2D* fhClusterEnergyVSEToF;
    TH2D* fhClusterSizeVSEToF;
    TH3D* fhClusterEnergyVSSizeVSEToF;

    TH2D* fhElasticTargetMass;

    TH2D* fhClusterForemostMinusCentroidVSEnergy;
    TH2D* fhClusterForemostMinusMaxEnergyDigiPosVSEnergy;

    TH2D* fhClusterCentroidMinusFirstDigiPosVSEnergy;
    TH2D* fhClusterMaxEnergyDigiMinusFirstDigiPosVSEnergy;
    TH2D* fhClusterMaxEnergyDigiMinusCentroidVSEnergy;
    TH2D* fhClusterEnergyMomentVSEnergy;
    TH2D* fhClusterEnergyMomentVSClusterSize;
    TH1D* fhClusterEnergyMoment;
    TH1D* fhClusterMaxEnergyDigiMinusFirstDigiMag;
    TH2D* fhClusterLastMinusFirstDigiMagVSEnergy;

    TH2D* fhENFromScatterVSEToF;

    TH2D* fhEToFVSEelastic;
    TH2D* fhScatteredNEnergyVSAngle;
    TH2D* fhScatteredNEnergyVSEdep;
    TH2D* fhScatterAngleVSRecoilAngle;
    TH2D* fhSumAngleVSRatioErecoEtof;
    TH2D* fhClusterEnergyVSScatteredRecoilAngle;
    TH2D* fhClusterEnergyVSScatteredNeutronAngle;

    TH1D* fhZ;
    TH2D* fhZVSEToF;
    TH1D* fhDistFromCenter;
    TH2D* fhDistFromCenterVSEToF;
    TH1D* fhDeltaT;
    TH1D* fhForemostMinusFirstDigiTime;

    TH2D* fhThetaEDigi;
    TH2D* fhThetaEDigiCosTheta;

    TH1D* hT;
    TH1D* hTNeigh;

    Filterable<R3BNeulandCluster*> fClusterFilters;

  public:
    ClassDefOverride(R3BNeulandClusterMon, 0);
};

#endif // R3BNEULANDCLUSTERMON_H
