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

#ifndef R3BCalifavsFootOnlineSpectra_H
#define R3BCalifavsFootOnlineSpectra_H 1

#include "FairTask.h"
#include "TCanvas.h"
#include "TMath.h"
#include <array>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <sstream>

#define Nb_Sides 2
#define Nb_Rings 5
#define Nb_Preamps 16
#define Nb_PreampCh 16

class TClonesArray;
class TH1F;
class TH1I;
class TH2F;
class R3BEventHeader;

class R3BCalifavsFootOnlineSpectra : public FairTask
{
  public:
    /**
     * Default constructor.
     * Creates an instance of the task with default parameters.
     */
    R3BCalifavsFootOnlineSpectra();

    /**
     * Standard constructor.
     * Creates an instance of the task.
     * @param name a name of the task.
     * @param iVerbose a verbosity level.
     */
    R3BCalifavsFootOnlineSpectra(const TString& name, Int_t iVerbose = 1);

    /**
     * Destructor.
     * Frees the memory used by the object.
     */
    virtual ~R3BCalifavsFootOnlineSpectra();

    /** Virtual method SetParContainers **/
    virtual void SetParContainers();

    /**
     * Method for task initialization.
     * This function is called by the framework before
     * the event loop.
     * @return Initialization status. kSUCCESS, kERROR or kFATAL.
     */
    virtual InitStatus Init();

    /** Virtual method ReInit **/
    virtual InitStatus ReInit();

    /**
     * Method for event loop implementation.
     * Is called by the framework every time a new event is read.
     * @param option an execution option.
     */
    virtual void Exec(Option_t* option);

    /**
     * A method for finish of processing of an event.
     * Is called by the framework for each event after executing
     * the tasks.
     */
    virtual void FinishEvent();

    /**
     * Method for finish of the task execution.
     * Is called by the framework after processing the event loop.
     */
    virtual void FinishTask();

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
     * Method to select binning and max range
     */
    inline void SetRange_bins(Int_t Histos_bins) { fMapHistos_bins = Histos_bins; }
    inline void SetRange_max(Int_t Histos_max) { fMapHistos_max = Histos_max; }

    /**
     * Method for setting min proton energy (in keV) for opening angle histogram
     */
    inline void SetMinProtonEnergyForOpening(Float_t min) { fMinProtonE = min; }

    /**
     * Method to reset histograms
     */
    void Reset_CALIFAFOOT_Histo();

    /**
     * Method for setting the fTpat
     */
    void SetTpat(Int_t tpat) { fTpat = tpat; }

  private:
    // void SetParameter();

    Int_t fMapHistos_max;
    Int_t fMapHistos_bins;
    TClonesArray* fCalItemsCalifa; /**< Array with cal items.       */
    TClonesArray* fHitItemsCalifa; /**< Array with hit items.       */
    TClonesArray* fMappedItemsFoot;
    TClonesArray* fCalItemsFoot;
    TClonesArray* fHitItemsFoot;

    // Check for trigger should be done globablly (somewhere else)
    R3BEventHeader* header; /**< Event header.  */
    Int_t fNEvents;         /**< Event counter. */
    Int_t fTpat;

    Int_t fNbCalifaCrystals;   /**< Number of Crystals in Califa. */
    Int_t fNumSides;           /**< Number of Sides, left and right.   */
    Int_t fNumRings;           /**< Number of Rings.   */
    Int_t fNumPreamps;         /**< Number of Preamps per ring.   */
    Int_t fNumCrystalPreamp;   /**< Number of Crystals/Channels per Preamp. */
    Int_t fBinsChannelFebex;   /**< Number of Bins per Febex channel. */
    Int_t fMaxBinChannelFebex; /**< Maximum bin for Febex histograms. */
    Float_t fMinProtonE;       /**< Min proton energy (in keV) to calculate the opening angle */

    Int_t fMaxEnergyBarrel; /**< Max. energy for Barrel histograms at CAL level. */
    Int_t fMaxEnergyIphos;  /**< Max. energy for Iphos histograms at CAL level. */
    Bool_t fRaw2Cal;        /**< Mapped or Cal selector. */

    TH2F* fh2_Califa_theta_phi[2]; // 0: all, 1: with foot
    TCanvas* cCalifa_angles;

  public:
    ClassDef(R3BCalifavsFootOnlineSpectra, 1)
};

#endif /* R3BCalifavsFootOnlineSpectra_H */
