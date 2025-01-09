/******************************************************************************
 *   Copyright (C) 2022 GSI Helmholtzzentrum für Schwerionenforschung GmbH    *
 *   Copyright (C) 2022-2025 Members of R3B Collaboration                     *
 *                                                                            *
 *             This software is distributed under the terms of the            *
 *                 GNU General Public Licence (GPL) version 3,                *
 *                    copied verbatim in the file "LICENSE".                  *
 *                                                                            *
 * In applying this license GSI does not waive the privileges and immunities  *
 * granted to it by virtue of its status as an Intergovernmental Organization *
 * or submit itself to any jurisdiction.                                      *
 ******************************************************************************/

// --------------------------------------------------------------
// -----          R3BAlpideDigitizer source file            -----
// -----    Created 12/10/22 by J.L. Rodriguez Sanchez      -----
// --------------------------------------------------------------

#pragma once

#include <FairTask.h>

#include "R3BAlpideHitData.h"

#include <Rtypes.h>
#include <TRotation.h>
#include <string>

constexpr double pixelSizeFactor = 150000.0;
constexpr double pixelSizeOffset = 0.8;

class TClonesArray;
class R3BTGeoPar;
class R3BAlpideGeometry;
class R3BAlpideMappingPar;

class R3BAlpideDigitizer : public FairTask
{
  public:
    /** Default constructor **/
    R3BAlpideDigitizer();

    /** Standard constructor **/
    R3BAlpideDigitizer(const TString& name, Int_t iVerbose = 1);

    /** Destructor **/
    ~R3BAlpideDigitizer() override;

    /** Virtual method Init **/
    InitStatus Init() override;

    /** Virtual method ReInit **/
    InitStatus ReInit() override;

    /** Virtual method Exec **/
    void Exec(Option_t*) override;

    // Fair specific
    void SetParContainers() override;

    /** Setter for sigma **/
    void SetSigma(Float_t sigma) { fsigma = sigma; }

    /** Setter for Lab. frame **/
    void SetLabframe() { fLabframe = true; }

  private:
    void SetParameter();
    void Reset();

    Int_t fGeoversion;
    TClonesArray* fMCTrack;
    TClonesArray* fAlpidePoints;
    TClonesArray* fAlpideHits;
    TString fName;
    Float_t fsigma;
    bool fLabframe;
    TVector3 fTrans;
    TRotation fRot;

    R3BAlpideMappingPar* fMappingPar;
    R3BAlpideGeometry* fAlpideGeo;

    /** Private method AddHitData **/
    R3BAlpideHitData* AddHitData(UInt_t sensorId, uint16_t clustersize, Double_t x, Double_t y, Double_t z = 0.0);

  public:
    // Class definition
    ClassDefOverride(R3BAlpideDigitizer, 1);
};
