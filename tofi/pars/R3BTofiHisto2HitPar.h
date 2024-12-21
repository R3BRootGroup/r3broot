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
// ----------------------------------------------------
// -----      Create hit parameters for TOFI      -----
// -----        Created Jul 2019 by L.Bott        -----
// ----------------------------------------------------

#ifndef R3BTOFIHISTO2HITPAR_H
#define R3BTOFIHISTO2HITPAR_H

#define N_TOFI_HIT_PLANE_MAX 1
#define N_TOFI_HIT_PADDLE_MAX 24

#include "FairTask.h"
#include "R3BTCalEngine.h"
#include "TObject.h"

class R3BTofiHitPar;
class TClonesArray;
class R3BEventHeader;
class TH1F;
class TH2F;

#ifdef __CINT__
#pragma link C++ class R3BTofiHitModulePar + ;
#endif
class R3BTofiHisto2HitPar : public FairTask
{
  public:
    /**
     * Default constructor.
     * Creates an instance of the task with default parameters.
     */
    R3BTofiHisto2HitPar();
    /**
     * Standard constructor.
     * Creates an instance of the task.
     * @param name a name of the task.
     * @param iVerbose a verbosity level.
     */
    R3BTofiHisto2HitPar(const char* name, Int_t iVerbose = 1);
    /**
     * Destructor.
     * Frees the memory used by the object.
     */
    virtual ~R3BTofiHisto2HitPar();
    /**
     * Method for task initialization.
     * This function is called by the framework before
     * the event loop.
     * @return Initialization status. kSUCCESS, kERROR or kFATAL.
     */
    virtual InitStatus Init();

    virtual void SetParContainers();

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
     * Method for setting the update rate for control histograms
     * @param rate an update rate value (events).
     */
    inline void SetUpdateRate(Int_t rate) { fUpdateRate = rate; }

    /**
     * Method for setting minimum required statistics per module.
     * Only detector modules with number of entries in TDC
     * distribution greater than minimum statistics will be
     * calibrated.
     * @param minStats a value of minimum statistics required.
     */
    inline void SetMinStats(Int_t minStats) { fMinStats = minStats; }

    /**
     * Method for selecting events with certain trigger value.
     * @param trigger 1 - onspill, 2 - offspill, -1 - all events.
     */
    inline void SetTrigger(Int_t trigger) { fTrigger = trigger; }

    /**
     * Method for selecting a certain paramter calculation.
     */
    inline void SetParameter(Int_t parameter) { fParameter = parameter; }

    /**
     * Method for setting number of LOS detectors and channels.
     * @param nDets number of detectors.
     * @param nCh number of channels per detector (4+master trigger?)
     */
    inline void SetNofModules(Int_t nDets, Int_t nCh)
    {
        fNofPlanes = nDets;
        fPaddlesPerPlane = nCh;
        fNofModules = nDets * nCh * 1; // 1 edge per channel
    }
    /**
     * Method for setting the y-position of a horizonzal sweep run for calibration of effective velocity of light
     */
    inline void SetTofiPos(Double_t Y) { fTofiY = Y; }
    /**
     * Method for setting the nuclear charge of main beam
     */
    inline void SetTofiQ(Double_t Q) { fTofiQ = Q; }
    /**
     * Method for setting the lower range of ToT for offset calibration
     */
    inline void SetTofiTotLow(Double_t TotLow) { fTofiTotLow = TotLow; }
    /**
     * Method for setting the upper range of ToT for offset calibration
     */
    inline void SetTofiTotHigh(Double_t TotHigh) { fTofiTotHigh = TotHigh; }
    /**
     * Method for using smiley or double exponential charge correction
     */
    inline void SetTofiSmiley(Bool_t Smiley) { fTofiSmiley = Smiley; }
    /**
     *
     */
    inline void ReadParaFile(TString file) { fParaFile = file; }
    /**
     *
     */
    inline void ReadHistoFile(TString file) { fHistoFile = file; }

    /**
     * Method for calculation of offset.
     */
    virtual void calcOffset();

    /**
     * Method for calculation of ToT offset.
     */
    virtual void calcToTOffset(Double_t TotLow, Double_t TotHigh);

    /**
     * Method for calculation of sync offset between paddles.
     */
    virtual void calcSync();

    /**
     * Method for calculation of effective speed of light.
     */
    virtual void calcVeff();

    /**
     * Method for calculation of light attenuation factor.
     */
    virtual void calcLambda(Double_t TotLow, Double_t TotHigh);

    /**
     * Method for calculation of position dependence charge of scintillator.
     */
    virtual void doubleExp(TH2F* histo, Double_t min, Double_t max, Double_t*);

    /**
     * Method for calculation of position dependent charge of scintillator.
     */
    virtual void smiley(TH2F* histo, Double_t min, Double_t max, Double_t*);

    /**
     * Method for calculation of z correction.
     */
    virtual void zcorr(TH2F* histo, Int_t min, Int_t max, Double_t*);

  private:
    Int_t fUpdateRate; /**< An update rate. */
    Int_t fMinStats;   /**< Minimum statistics required per module. */
    Int_t fTrigger;    /**< Trigger value. */
    Int_t fParameter;  /**< Parameter #. */

    UInt_t fNofPlanes;       /**< Number of planes. */
    UInt_t fPaddlesPerPlane; /**< Number of bars per plane. */
    UInt_t fNofModules;      /**< Total number of modules (=edges) to calibrate */

    Int_t fNEvents;             /**< Event counter. */
    R3BTofiHitPar* fCal_Par;    /**< Parameter container. */
    TClonesArray* fCalItemsLos; /**< Array with cal items. */
    TClonesArray* fCalData;     /**< Array with mapped data - input data. */
    R3BEventHeader* header;     /**< Event header - input data. */
    Double_t fClockFreq;        /**< Clock cycle in [ns]. */
    Double_t fTofiY;
    Double_t fTofiQ;
    Double_t fTofiTotLow;
    Double_t fTofiTotHigh;
    Int_t countskip = 0;
    Bool_t fTofiSmiley;
    TString fParaFile;
    TString fHistoFile;

  public:
    ClassDef(R3BTofiHisto2HitPar, 1)
};

#endif
