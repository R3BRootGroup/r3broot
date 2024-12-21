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

#ifndef R3BCalifaCrystalCal2TotCalPar_H
#define R3BCalifaCrystalCal2TotCalPar_H 1

#include <FairTask.h>
#include <TGraph.h>
#include <TH1F.h>

class TClonesArray;
class R3BCalifaMappingPar;
class R3BCalifaTotCalPar;
class R3BEventHeader;

class R3BCalifaCrystalCal2TotCalPar : public FairTask
{

  public:
    /** Default constructor **/
    R3BCalifaCrystalCal2TotCalPar();

    /** Standard constructor **/
    R3BCalifaCrystalCal2TotCalPar(const char* name, Int_t iVerbose = 1);

    /** Destructor **/
    virtual ~R3BCalifaCrystalCal2TotCalPar();

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

    /** Virtual method SetParContainers **/
    virtual void SetParContainers();

    virtual void Search_TotParams();

    /** Accessor functions **/
    const Int_t GetNumCrystals() { return fNumCrystals; }
    const Double_t GetThreshold() { return fThreshold; }
    const Int_t GetNumParameterFit() { return fNumParam; }
    const Int_t GetMinStadistics() { return fMinStadistics; }

    void SetNumCrystals(Int_t numberCry) { fNumCrystals = numberCry; }
    void SetThreshold(Double_t threshold) { fThreshold = threshold; }
    void SetNumParameterFit(Int_t numberParFit) { fNumParam = numberParFit; }
    void SetMinStadistics(Int_t minstad) { fMinStadistics = minstad; }
    void SetDebugMode(Bool_t debug) { fDebugMode = debug; }
    void SetRange(Double_t H_left, Double_t H_right)
    {
        fLeft = H_left;
        fRight = H_right;
    }

  private:
    void SetParameter();
    Bool_t fDebugMode;
    Int_t fNumCrystals;
    Int_t fNumParam;
    Int_t fMinStadistics;
    Int_t* idx;

    Double_t fThreshold;
    Double_t fLeft;
    Double_t fRight;

    R3BCalifaMappingPar* fMap_Par;   /**< Parameter container with mapping. >*/
    R3BCalifaTotCalPar* fTotCal_Par; /**< Container for Tot-Cal parameters. >*/
    TClonesArray* fCrystalCalDataCA; /**< Array with CALIFA energy calibrated - input data. >*/

    TGraph** energy_vs_tot_crystal;

  public:
    ClassDef(R3BCalifaCrystalCal2TotCalPar, 1);
};

#endif
