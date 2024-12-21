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

#ifndef R3BCALIFAMAPPED2CRYSTALCALPAR_H
#define R3BCALIFAMAPPED2CRYSTALCALPAR_H 1

#include <FairTask.h>
#include <TArrayF.h>
#include <TH1F.h>

class TClonesArray;
class R3BCalifaMappingPar;
class R3BCalifaCrystalCalPar;
class R3BEventHeader;

class R3BCalifaMapped2CrystalCalPar : public FairTask
{

  public:
    /** Default constructor **/
    R3BCalifaMapped2CrystalCalPar();

    /** Standard constructor **/
    R3BCalifaMapped2CrystalCalPar(const char* name, Int_t iVerbose = 1);

    /** Destructor **/
    virtual ~R3BCalifaMapped2CrystalCalPar();

    /** Virtual method Init **/
    virtual InitStatus Init();

    /** Virtual method Exec **/
    virtual void Exec(Option_t* opt);

    /** Virtual method FinishEvent **/
    virtual void FinishEvent();

    /** Virtual method FinishTask **/
    virtual void FinishTask();

    /** Virtual method Reset **/
    virtual void Reset();

    /** Virtual method ReInit **/
    virtual InitStatus ReInit();

    /** Virtual method Search peaks and calibrate **/
    virtual void SearchPeaks();

    /** Virtual method SetParContainers **/
    virtual void SetParContainers();

    /** Accessor functions **/
    const Int_t GetNumCrystals() { return fNumCrystals; }
    const Int_t GetCalRange_left() { return fMapHistos_left; }
    const Int_t GetCalRange_right() { return fMapHistos_right; }
    const Int_t GetCalRange_bins() { return fMapHistos_bins; }
    const Int_t GetNumPeaks() { return fNumPeaks; }
    const Double_t GetSigma() { return fSigma; }
    const Double_t GetThreshold() { return fThreshold; }
    const Int_t GetNumParameterFit() { return fNumParam; }
    const Int_t GetMinStadistics() { return fMinStadistics; }

    TArrayF* GetEnergyPeaks() { return fEnergyPeaks; }

    void SetNumCrystals(Int_t numberCry) { fNumCrystals = numberCry; }
    void SetCalRange_left(Int_t Histos_left) { fMapHistos_left = Histos_left; }
    void SetCalRange_right(Int_t Histos_right) { fMapHistos_right = Histos_right; }
    void SetCalRange_bins(Int_t Histos_bins) { fMapHistos_bins = Histos_bins; }
    void SetCalRangeP_left(Int_t Histos_left) { fMapHistos_leftp = Histos_left; }
    void SetCalRangeP_right(Int_t Histos_right) { fMapHistos_rightp = Histos_right; }
    void SetCalRangeP_bins(Int_t Histos_bins) { fMapHistos_binsp = Histos_bins; }
    void SetNumPeaks(Int_t numberpeaks) { fNumPeaks = numberpeaks; }
    void SetSigma(Double_t sigma) { fSigma = sigma; }
    void SetThreshold(Double_t threshold) { fThreshold = threshold; }
    void SetNumParameterFit(Int_t numberParFit) { fNumParam = numberParFit; }
    void SetMinStadistics(Int_t minstad) { fMinStadistics = minstad; }

    void SetDebugMode(Bool_t debug) { fDebugMode = debug; }

    void SetEnergyPeaks(TArrayF* thePeaks)
    {
        fEnergyPeaks = thePeaks;
        fNumPeaks = thePeaks->GetSize();
    }

  private:
    void SetParameter();
    Bool_t fDebugMode;
    Int_t fNumCrystals;
    Int_t fMapHistos_left; // gamma range
    Int_t fMapHistos_right;
    Int_t fMapHistos_bins;
    Int_t fMapHistos_leftp; // particle range
    Int_t fMapHistos_rightp;
    Int_t fMapHistos_binsp;

    Int_t fNumParam;
    Int_t fMinStadistics;

    Int_t fNumPeaks;
    Double_t fSigma;
    Double_t fThreshold;

    TArrayF* fEnergyPeaks;
    Double_t* fChannelPeaks;

    R3BCalifaMappingPar* fMap_Par;     /**< Parameter container with mapping. >*/
    R3BCalifaCrystalCalPar* fCal_Par;  /**< Container for Cal parameters. >*/
    TClonesArray* fCalifaMappedDataCA; /**< Array with CALIFA Mapped-input data. >*/

    TH1F** fh_Map_energy_crystal;

  public:
    ClassDef(R3BCalifaMapped2CrystalCalPar, 2);
};

#endif
