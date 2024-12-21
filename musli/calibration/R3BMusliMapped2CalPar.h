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

#ifndef R3BMusliMapped2CalPar_H
#define R3BMusliMapped2CalPar_H 1

#include "FairTask.h"
#include "TGraph.h"
#include "TH1F.h"
#include "TH2D.h"

class TClonesArray;
class R3BMusliCalPar;
class R3BTGeoPar;

class R3BMusliMapped2CalPar : public FairTask
{
  public:
    /** Default constructor **/
    R3BMusliMapped2CalPar();

    /** Standard constructor **/
    R3BMusliMapped2CalPar(const TString& name,
                          Int_t iVerbose = 1,
                          const TString& namedeta = "Mwpc1",
                          const TString& namedetb = "Mwpc0");

    /** Destructor **/
    virtual ~R3BMusliMapped2CalPar();

    /** Virtual method Exec **/
    virtual void Exec(Option_t* option);

    /** Virtual method FinishTask **/
    virtual void FinishTask();

    virtual void SetParContainers();

    /** Virtual method Init **/
    virtual InitStatus Init();

    /** Virtual method ReInit **/
    virtual InitStatus ReInit();

    /** Method to set up the limits for fit **/
    void SetFitLimits(Int_t min, Int_t max)
    {
        fLim_DTrawMin = min;
        fLim_DTrawMax = max;
    }

    void SetNumSignals(Int_t n) { fNumSignals = n; }
    void SetNumGroupsAnodes(Int_t n) { fNumGroupsAnodes = n; }
    void SetMaxMult(Int_t n) { fMaxMult = n; }
    void SetMinStatistics(Int_t n) { fMinStatistics = n; }
    void SetNumParamsEneFit(Int_t n) { fNumParamsEneFit = n; }
    void SetNumParamsPosFit(Int_t n) { fNumParamsPosFit = n; }
    virtual void CalculatePosCalPar();
    // virtual void CalculateEneCalPar();

  private:
    Int_t fNumSignals;
    Int_t fNumGroupsAnodes;
    Int_t fMaxMult;
    Int_t fMinStatistics;
    Int_t fNumParamsEneFit;
    Int_t fNumParamsPosFit;
    Int_t fLim_DTrawMin;
    Int_t fLim_DTrawMax;
    TArrayF* EneCalParams; // for energy loss alignement
    TArrayF* PosCalParams; // for position calibration (interpolation from Mw)

    TString fNameDetA;
    TString fNameDetB;

    Int_t multMap[18];         // mapped data
    Double_t fEneMap[20][18];  // mapped data
    Double_t fTimeMap[20][18]; // mapped data

    Double_t fXA; // hit data of MwA
    Double_t fXB; // hit data of MwB

    R3BMusliCalPar* fCal_Par; /**< Parameter container. >*/
    R3BTGeoPar* fMusliGeo_Par;
    R3BTGeoPar* fMwAGeo_Par;
    R3BTGeoPar* fMwBGeo_Par;
    TClonesArray* fMusliMappedDataCA; /**< Array with Musli Mapped-input data. >*/
    TClonesArray* fMwAHitDataCA;      /**< Array with hit items. */
    TClonesArray* fMwBHitDataCA;      /**< Array with hit items. */

    TH2D** fh2_XvsDT;

  public:
    // Class definition
    ClassDef(R3BMusliMapped2CalPar, 1)
};

#endif
