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

// -------------------------------------------------------------------------
// -----         R3BAmsMapped2StripCalPar source file                  -----
// -----             Created 29/05/18  by J.L. Rodriguez-Sanchez       -----
// -------------------------------------------------------------------------

#ifndef R3BAmsMapped2StripCalPar_H
#define R3BAmsMapped2StripCalPar_H

#include "FairTask.h"
#include "TH1F.h"

class TClonesArray;
class R3BAmsMappingPar;
class R3BAmsStripCalPar;
class R3BEventHeader;

class R3BAmsMapped2StripCalPar : public FairTask
{

  public:
    /** Default constructor **/
    R3BAmsMapped2StripCalPar();

    /** Standard constructor **/
    R3BAmsMapped2StripCalPar(const TString& name, Int_t iVerbose = 1);

    /** Destructor **/
    virtual ~R3BAmsMapped2StripCalPar();

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

    /** Virtual method Search pedestals **/
    virtual void SearchPedestals();

    /** Virtual method to print parameters **/
    virtual void PrintParamsDaq();

    /** Virtual method SetParContainers **/
    virtual void SetParContainers();

    void SetPrintParamsDaq(Bool_t print) { fPrint = print; }

    /** Accessor functions **/
    const Int_t GetNumDetectors() { return fNumDets; }
    const Int_t GetNumStrips() { return fNumStrips; }
    const Int_t GetNumStripsS() { return fNumStripsS; }
    const Int_t GetNumStripsK() { return fNumStripsK; }
    const Int_t GetCalRange_left() { return fMapHistos_left; }
    const Int_t GetCalRange_right() { return fMapHistos_right; }
    const Int_t GetCalRange_bins() { return fMapHistos_bins; }
    const Double_t GetSigma() { return fSigma; }
    const Double_t GetMean() { return fMean; }
    const Int_t GetMinStadistics() { return fMinStadistics; }
    const Int_t GetMaxSigma() { return fMaxSigma; }

    void SetNumDetectors(Int_t numberDet) { fNumDets = numberDet; }
    void SetNumStrips(Int_t numberStrip) { fNumStrips = numberStrip; }
    void SetNumStripsS(Int_t numberStrip) { fNumStripsS = numberStrip; }
    void SetNumStripsK(Int_t numberStrip) { fNumStripsK = numberStrip; }
    void SetCalRange_left(Int_t Histos_left) { fMapHistos_left = Histos_left; }
    void SetCalRange_right(Int_t Histos_right) { fMapHistos_right = Histos_right; }
    void SetCalRange_bins(Int_t Histos_bins) { fMapHistos_bins = Histos_bins; }
    void SetSigma(Double_t sigma) { fSigma = sigma; }
    void SetMean(Double_t mean) { fMean = mean; }
    void SetMinStadistics(Int_t minstad) { fMinStadistics = minstad; }
    void SetMaxSigma(Double_t sigma) { fMaxSigma = sigma; }

  protected:
    void SetParameter();
    // Number of histograms, limits and bining
    Int_t fNumDets;
    Int_t fNumStrips;
    Int_t fNumStripsS;
    Int_t fNumStripsK;
    Int_t fMapHistos_left;
    Int_t fMapHistos_right;
    Int_t fMapHistos_bins;
    Int_t fMaxSigma; // Maximum sigma to mark dead strips, pedestal -1
    Double_t parameters[8 * 1024][2];

    // Minimum stadistics and parameters
    Int_t fMinStadistics;
    Double_t fSigma;
    Double_t fMean;
    Bool_t fPrint;

    R3BAmsMappingPar* fMap_Par;     /**< Parameter container with mapping. >*/
    R3BAmsStripCalPar* fStrip_Par;  /**< Parameter container. >*/
    TClonesArray* fAmsMappedDataCA; /**< Array with AMS Mapped input data. >*/

    TH1F** fh_Map_energy_strip;

  public:
    ClassDef(R3BAmsMapped2StripCalPar, 1);
};

#endif
