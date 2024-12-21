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

#ifndef R3BNEULANDMCMON_H
#define R3BNEULANDMCMON_H

/** Neuland Monte Carlo Monitor
 * @author Jan Mayer
 * @since  07.07.2015
 *
 *  Input:  Monte Carlo Tracks "MCTrack".
 *  Output: Several histograms:
 *          - Energy of primary neutrons
 *          - Number of particles by PDG code created by the first neutron interaction
 *          - Energy of non-neutron tracks in kLAND created by primary neutron interaction(s), by PID
 *          - Energy of neutron tracks created by primary neutron interaction
 *          - Total energy of non-neutron tracks in kLAND created by primary neutron interaction(s)
 *          - Total energy of non-neutron tracks in kLAND created by primary neutron interaction(s), by PID
 *          - IDs of tracks with a primary mother
 *          - Distribution of track mother IDs
 *          - ...
 *          - Option FULLSIMANA: Create histograms for reactions products of primary neutron interactions
 *            Grouped by Pdg, Reaction and Reaction&Pdg (several 1000 histograms)
 *          - Option 3DTRACK: Prepare 3D Histogram for each event.
 *  Usage:  run->AddTask(new R3BNeulandMCMon("FULLSIMANA"));
 *          Requires MCTracks, NeulandPoints and NPNIPS from R3BNeuland simulation.
 */

#include "FairTask.h"
#include "R3BMCTrack.h"
#include "R3BNeulandPoint.h"
#include "TCAConnector.h"
#include <map>

class TH1D;
class TH2D;
class TH3D;

class R3BNeulandMCMon : public FairTask
{
  public:
    R3BNeulandMCMon(const Option_t* option = "");

    ~R3BNeulandMCMon() override = default;

    // No copy and no move is allowed (Rule of three/five)
    R3BNeulandMCMon(const R3BNeulandMCMon&) = delete;            // copy constructor
    R3BNeulandMCMon(R3BNeulandMCMon&&) = delete;                 // move constructor
    R3BNeulandMCMon& operator=(const R3BNeulandMCMon&) = delete; // copy assignment
    R3BNeulandMCMon& operator=(R3BNeulandMCMon&&) = delete;      // move assignment

  protected:
    InitStatus Init() override;
    void Finish() override;

  public:
    void Exec(Option_t*) override;

  private:
    Bool_t fIs3DTrackEnabled;
    Bool_t fIsFullSimAnaEnabled;

    TCAInputConnector<R3BNeulandPoint> fPrimaryNeutronInteractionPoints;
    TCAInputConnector<R3BMCTrack> fMCTracks;
    TCAInputConnector<R3BNeulandPoint> fNeulandPoints;

    int nEvents;

    TH1D* fhPDG;
    TH1D* fhEPrimarys;
    TH1D* fhEPrimarys2;
    TH1D* fhEPrimaryNeutrons;
    TH1D* fhErelMC;
    TH1D* fhEtot;
    TH1D* fhEtotPrim;
    TH1D* fhESecondaryNeutrons;
    TH1D* fhMotherIDs;
    TH1D* fhPrimaryDaughterIDs;
    TH1D* fhMCToF;
    TH1D* fhnNPNIPs;
    TH2D* fhNPNIPsEToFVSTime;
    TH2D* fhNPNIPSrvsz;
    TH2D* fhNPNIPSxy;
    TH2D* fhThetaLight;
    TH2D* fhElossVSLight;
    TH2D* fhElossVSLightLog;
    std::map<Int_t, TH1D*> fhmEPdg;
    std::map<Int_t, TH1D*> fhmEtotPdg;
    std::map<Int_t, TH1D*> fhmEtotPdgRel;

    std::map<Int_t, TH1D*> fhmEnergyByProductPdg;
    std::map<Int_t, TH1D*> fhmCountByProductPdg;
    std::map<std::string, TH1D*> fhmEnergyByProductPdgReduced;
    std::map<std::string, TH1D*> fhmCountByProductPdgReduced;
    std::map<Int_t, TH2D*> fhmElossVSLightLogPdg;
    std::map<TString, TH1D*> fhmEnergyByReaction;
    std::map<TString, std::map<Int_t, TH1D*>> fhmmEnergyByReactionByProductPdg;
    TH1D* fhnProducts;
    TH1D* fhnProductsCharged;
    TH1D* fhSumProductEnergy;
    TH1D* fhnSecondaryProtons;
    TH1D* fhnSecondaryNeutrons;

    TH3D* fh3;
    TH3D* fh3PNIP;

    ClassDefOverride(R3BNeulandMCMon, 0)
};

#endif // R3BNEULANDMCMON_H
