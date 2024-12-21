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

#include <FairDetector.h>
#include <TGeoMatrix.h>
#include <TString.h>

class R3BDetector : public FairDetector
{
  public:
    R3BDetector();

    R3BDetector(const TString& name, const Int_t detId);

    R3BDetector(const TString& name,
                const Int_t detId,
                const TString& geoFile,
                const TGeoTranslation& trans,
                const TGeoRotation& rot = TGeoRotation());

    R3BDetector(const TString& name,
                const Int_t detId,
                const TString& geoFile,
                const TGeoCombiTrans& combi = TGeoCombiTrans());

    /** Construct geometry from root files */
    void ConstructGeometry() override;
    void ConstructRootGeometry(TGeoMatrix* geo = nullptr) override;

    /** Method to setup the position of the detector-origin. */
    virtual void SetPosition(const TGeoTranslation& trans) { fCombiTrans.SetTranslation(trans); }

    /** Method to setup the rotation of the detector */
    virtual void SetRotation(const TGeoRotation& rot) { fCombiTrans.SetRotation(rot); }

    /** Method to setup the translation + rotation of the detector */
    virtual void SetCombiTrans(const TGeoCombiTrans& combi) { fCombiTrans = combi; }

    virtual void SetEnergyCut(const Double_t cutE) { fCutE = cutE; }

    virtual void SetActive(const Bool_t flag) { fActive = flag; }

  protected:
    TGeoCombiTrans fCombiTrans;
    Double_t fCutE;

  public:
    ClassDefOverride(R3BDetector, 3)
};
