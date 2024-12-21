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

// ---------------------------------------------------------------------
// -----            R3BFootStripCal2Hit source file                -----
// -----       Created 05/11/21 by J.L. Rodriguez-Sanchez          -----
// ---------------------------------------------------------------------

#pragma once

#include "R3BFootHitData.h"

#include <FairTask.h>

#include <Rtypes.h>
#include <TVector3.h>
#include <vector>

class TClonesArray;
class TH1F;
class R3BFootMappingPar;

class R3BFootStripCal2Hit : public FairTask
{

  public:
    /** Default constructor **/
    R3BFootStripCal2Hit();

    /** Standard constructor **/
    R3BFootStripCal2Hit(const TString& name, Int_t iVerbose = 1);

    /** Destructor **/
    virtual ~R3BFootStripCal2Hit();

    /** Virtual method Exec **/
    void Exec(Option_t* /*option*/) override;

    /** Virtual method Reset **/
    virtual void Reset();

    // Fair specific
    /** Virtual method Init **/
    InitStatus Init() override;

    /** Virtual method ReInit **/
    InitStatus ReInit() override;

    /** Virtual method SetParContainers **/
    void SetParContainers() override;

    /** Accessor for selecting online mode **/
    inline void SetOnline(Bool_t option) { fOnline = option; }

    /** Accessor for selecting max. number of clusters per ams detector **/
    inline void SetMaxNumClusters(int max) { fMaxNumClusters = max; }

    /** Accessor to set up the threshold for the cluster energy sum **/
    inline void SetClusterEnergy(double thsum) { fThSum = thsum; }

  private:
    void SetParameter();

    double fPitch = 157.7;
    double fMiddle = 50.;
    double fThSum = 20.;
    int fMaxNumDet = 16;
    int fMaxNumClusters = 10;
    std::vector<double> fDistTarget;
    std::vector<double> fAngleTheta;
    std::vector<double> fAnglePhi;
    std::vector<double> fOffsetX;
    std::vector<double> fOffsetY;
    std::vector<TH1F*> hssd;

    std::vector<int> ClusterMult;                 // Cluster multiplicity
    std::vector<std::vector<double>> ClusterPos;  // Position of Cluster from Weighted Average
    std::vector<std::vector<double>> ClusterESum; // Sum of Energies in the Cluster
    // std::vector<std::vector<double>> Nu;       // Nu for Energy/Position correction
    std::vector<std::vector<int>> ClusterNStrip;            // Number of Strips in Cluster
    std::vector<std::vector<std::vector<double>>> ClusterI; // Id Strip in Cluster
    std::vector<std::vector<std::vector<double>>> ClusterE; // Energy of Strip in Cluster

    R3BFootMappingPar* fMap_Par = nullptr; // Parameter container with mapping
    TClonesArray* fFootCalData = nullptr;  // Array with FOOT Cal-input data
    TClonesArray* fFootHitData = nullptr;  // Array with FOOT Hit-output data

    bool fOnline = false; // Don't store data for online
    Double_t* fChannelPeaks;

    // Private method AddHitData
    R3BFootHitData* AddHitData(uint8_t detid,
                               uint16_t numhit,
                               double s,
                               TVector3 master,
                               double energy_s,
                               uint16_t mulS);

  public:
    // Class definition
    ClassDefOverride(R3BFootStripCal2Hit, 1)
};
