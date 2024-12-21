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

#include <FairTask.h>
#include <TCanvas.h>
#include <TMath.h>
#include <array>
#include <boost/multi_array.hpp>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <sstream>
#include <vector>

#include "R3BShared.h"

using namespace boost;

constexpr const int Nb_Sides = 2;
constexpr const int Nb_Rings = 6;
constexpr const int Nb_Preamps = 16;
constexpr const int Nb_PreampCh = 16;
constexpr const int Nb_SlotandModule = 4; // Febex slot and module info: 0 slot and 1 module, (PR) 2 slot and 3 module
constexpr const int BinsChannelFebex = 5000; // Number of Bins per Febex channel
constexpr const int MaxBinChannelFebex = 65535;
constexpr const int MaxNbCrystals = 5088;  // gamma + proton range channels
constexpr const int BarrelCrystals = 3904; // gamma + proton range channels
constexpr const int iPhosCrystals = 4864;  // gamma + proton range channels

class TClonesArray;
class R3BCalifaMappingPar;
class TH1F;
class TH1I;
class TH2F;
class R3BEventHeader;

/**
 * This taks reads CALIFA data and plots online histograms
 */
class R3BCalifaOnlineSpectra : public FairTask
{
  public:
    /**
     * Default constructor.
     * Creates an instance of the task with default parameters.
     */
    R3BCalifaOnlineSpectra();

    /**
     * Standard constructor.
     * Creates an instance of the task.
     * @param name a name of the task.
     * @param iVerbose a verbosity level.
     */
    explicit R3BCalifaOnlineSpectra(const TString& name, Int_t iVerbose = 1);

    /**
     * Destructor.
     * Frees the memory used by the object.
     */
    ~R3BCalifaOnlineSpectra() = default;

    /** Virtual method SetParContainers **/
    void SetParContainers() override;

    /**
     * Method for task initialization.
     * This function is called by the framework before
     * the event loop.
     * @return Initialization status. kSUCCESS, kERROR or kFATAL.
     */
    InitStatus Init() override;

    /** Virtual method ReInit **/
    InitStatus ReInit() override;

    /**
     * Method for event loop implementation.
     * Is called by the framework every time a new event is read.
     * @param option an execution option.
     */
    void Exec(Option_t* /*option*/) override;

    /**
     * A method for finish of processing of an event.
     * Is called by the framework for each event after executing
     * the tasks.
     */
    void FinishEvent() override;

    /**
     * Method for finish of the task execution.
     * Is called by the framework after processing the event loop.
     */
    void FinishTask() override;

    /**
     * Method for setting number of rings
     */
    inline void SetRings(Int_t rings) { fNumRings = rings; }

    /**
     * Method for setting max. energy per crystal in MeV for Barrel histograms at CAL level
     */
    inline void SetMaxEnergyperCrystalBarrel(Int_t maxenergy) { fMaxEnergyBarrel = maxenergy; }

    /**
     * Method for setting max. energy per crystal in MeV for Iphos histograms at CAL level
     */
    inline void SetMaxEnergyperCrystalIphos(Int_t maxenergy) { fMaxEnergyIphos = maxenergy; }

    /**
     * Method for setting the configuration parameters file
     */
    inline void SetCalifaConfigFile(TString file) { fCalifaFile = file; }

    /**
     * Method to select binning and max range
     */
    inline void SetRange_bins(Int_t Histos_bins) { fMapHistos_bins = Histos_bins; }
    inline void SetRange_max(Int_t Histos_max) { fMapHistos_max = Histos_max; }

    /**
     * Method for setting the number of bins of Febex histograms
     */
    inline void SetBinChannelFebex(Int_t bin) { fBinsChannelFebex = bin; }

    /**
     * Method for setting max range of Febex histograms
     */
    inline void SetMaxBinFebex(Int_t max) { fMaxBinChannelFebex = max; }

    /**
     * Method for setting min proton energy (in keV) for opening angle histogram
     */
    inline void SetMinProtonEnergyForOpening(Float_t min) { fMinProtonE = min; }

    /**
     * Method for setting Tot histograms
     */
    inline void SetTotHist(Bool_t opt) { fTotHist = opt; }

    /**
     * Method to reset histograms
     */
    void Reset_CALIFA_Histo();

    /**
     * Method to change histogram scales
     */
    void Log_CALIFA_Histo();

    /**
     * Method for setting histogram sequence (Febex or Preamp. channels)
     */
    void Febex2Preamp_CALIFA_Histo();

    /**
     * Method for setting the trigger
     */
    inline void SetTrigger(Int_t trigger) { fTrigger = trigger; }

    /**
     * Method for selecting tpat values.
     */
    inline void SetTpat(Int_t tpat1, Int_t tpat2)
    {
        fTpat1 = tpat1;
        fTpat2 = tpat2;
    }

  private:
    void SetParameter();

    int fMapHistos_max = 4000;
    int fMapHistos_bins = 500;
    int fTpat1 = -1;
    int fTpat2 = -1;

    R3BCalifaMappingPar* fMap_Par = nullptr;    /**< Container with mapping parameters. >*/
    TClonesArray* fMappedItemsCalifa = nullptr; /**< Array with mapped items.    */
    TClonesArray* fTrigMappedItemsCalifa = nullptr;
    TClonesArray* fCalItemsCalifa = nullptr; /**< Array with cal items.       */
    TClonesArray* fHitItemsCalifa = nullptr; /**< Array with hit items.       */
    TClonesArray* fWRItemsCalifa = nullptr;  /**< Array with WR-Califa items. */
    TClonesArray* fWRItemsMaster = nullptr;  /**< Array with WR-Master items. */

    // Check for trigger should be done globablly (somewhere else)
    R3BEventHeader* header = nullptr; /**< Event header.  */
    int fNEvents = 0;                 /**< Event counter. */
    int fTrigger = -1;

    int fNbCalifaCrystals = MaxNbCrystals;        /**< Number of Crystals in Califa. */
    int fNumSides = Nb_Sides;                     /**< Number of Sides, left and right.   */
    int fNumRings = Nb_Rings;                     /**< Number of Rings.   */
    int fNumPreamps = Nb_Preamps;                 /**< Number of Preamps per ring.   */
    int fNumCrystalPreamp = Nb_PreampCh;          /**< Number of Crystals/Channels per Preamp. */
    int fBinsChannelFebex = BinsChannelFebex;     /**< Number of Bins per Febex channel. */
    int fMaxBinChannelFebex = MaxBinChannelFebex; /**< Maximum bin for Febex histograms. */

    // Selector for febex or preamp sequence
    const std::array<int, Nb_PreampCh> fOrderFebexPreamp{ 6, 5, 4, 3, 2, 1, 0, 7, 8, 15, 14, 13, 12, 11, 10, 9 };

    float fMinProtonE = 50000.; /**< Min proton energy (in keV) to calculate the opening angle */

    TString fCalifaFile;       /**< Config file name. */
    int fMaxEnergyBarrel = 10; /**< Max. energy for Barrel histograms at CAL level. */
    int fMaxEnergyIphos = 30;  /**< Max. energy for Iphos histograms at CAL level. */
    bool fLogScale = true;     /**< Selecting scale. */
    bool fRaw2Cal = false;     /**< Mapped or Cal selector. */
    bool fFebex2Preamp = true; /**< Febex or Preamp selector. */
    bool fTotHist = false;     /**< Tot histograms selector. */
    multi_array<int, 4> fFebexInfo;

    // Canvas
    TCanvas* cCalifaMult;
    TCanvas* cCalifa_cry_energy;
    TCanvas* cCalifa_cry_energy_cal;
    TCanvas* cMap_RingR[Nb_Rings];
    TCanvas* cMap_RingL[Nb_Rings];
    TCanvas* cMapCry[Nb_Sides][Nb_Rings][Nb_Preamps];
    TCanvas* cMapCryTot[Nb_Sides][Nb_Rings][Nb_Preamps];
    TCanvas* cMapCryCal[Nb_Sides][Nb_Rings][Nb_Preamps];
    TCanvas* cMapCryP[Nb_Sides][Nb_Rings][Nb_Preamps];
    TCanvas* cMapCryPTot[Nb_Sides][Nb_Rings][Nb_Preamps];
    TCanvas* cMapCryPCal[Nb_Sides][Nb_Rings][Nb_Preamps];
    TCanvas* cCalifaCoinE;
    TCanvas* cCalifaCoinPhi;
    TCanvas* cCalifaCoinTheta;
    TCanvas* cCalifa_angles;
    TCanvas* cCalifa_theta_energy;
    TCanvas* cCalifa_hitenergy;
    TCanvas* cCalifa_opening;
    TCanvas* cCalifa_NsNf;
    TCanvas* cCalifaTriggers;
    TCanvas* cCalifa_opening_tpat;

    // WR data
    TCanvas* cCalifa_wr;
    TH1I* fh1_Califa_wr;
    TCanvas* cWrs;
    TH1I* fh1_wrs[2];
    TCanvas* cCalifa_sync;
    TH1F* fh1_Califa_sync[3];
    TCanvas* cCalifa_wr_energy;

    // Histograms
    TH1F* fh1_Califa_Mult;
    TH1F* fh1_Califa_MultHit;
    TH2F* fh2_Califa_cryId_energy;
    TH2F* fh2_Preamp_vs_ch_R[Nb_Rings];
    TH2F* fh2_Preamp_vs_ch_L[Nb_Rings];
    TH1F* fh1_crystals[Nb_Sides][Nb_Rings][Nb_Preamps][Nb_PreampCh];
    TH2F* fh2_crystalsETot[Nb_Sides][Nb_Rings][Nb_Preamps][Nb_PreampCh];
    TH1F* fh1_crystals_p[Nb_Sides][Nb_Rings][Nb_Preamps][Nb_PreampCh];
    TH2F* fh2_crystalsETot_p[Nb_Sides][Nb_Rings][Nb_Preamps][Nb_PreampCh];
    TH2F* fh2_Califa_cryId_energy_cal;
    TH1F* fh1_crystals_cal[Nb_Sides][Nb_Rings][Nb_Preamps][Nb_PreampCh];
    TH1F* fh1_crystals_p_cal[Nb_Sides][Nb_Rings][Nb_Preamps][Nb_PreampCh];
    TH2F* fh2_Califa_coinE;
    TH2F* fh2_Califa_coinE_p2p;
    TH2F* fh2_Califa_coinTheta;
    TH2F* fh2_Califa_coinTheta_cutOPA;
    TH2F* fh2_Califa_coinPhi;
    TH2F* fh2_Califa_theta_phi;
    TH2F* fh2_Califa_theta_energy;
    TH1F* fh1_Califa_total_energy;
    TH1F* fh1_openangle;
    TH2F* fh2_openangle_tpat;
    TH2F* fh2_Cal_wr_energy_l;
    TH2F* fh2_Cal_wr_energy_r;
    std::vector<TH2F*> fh2_Califa_NsNf;
    TH2F* fh2_Califa_EtrigCor[4];
    TH1F* fh1_Califa_Etrig[2];
    TH1F* fh1_CalifaTriggers;

  public:
    ClassDefOverride(R3BCalifaOnlineSpectra, 1)
};
