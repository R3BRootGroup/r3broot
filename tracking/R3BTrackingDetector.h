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

#ifndef R3B_TRACKING_DETECTOR
#define R3B_TRACKING_DETECTOR

#include <vector>

#include "FairTask.h"
#include "R3BTrackingParticle.h"
#include "TObject.h"
#include "TString.h"
#include "TVector3.h"

enum EDetectorType
{
    kBeforeTarget,
    kTarget,
    kTargetGlad,
    kAfterGlad,
    kTof
};

class R3BHit;
class TClonesArray;
class R3BTGeoPar;

/* Generic detector class that holds all infos needed for the track fitting
 */
class R3BTrackingDetector : public TObject
{
  public:
    R3BTrackingDetector(const char* detectorName,
                        EDetectorType type,
                        const char* geoParName,
                        const char* hitArray = NULL);
    virtual ~R3BTrackingDetector();

    void SetParContainers();

    InitStatus Init();

    void CopyHits();

    void Draw(Option_t* option = "");

    void GlobalToLocal(const TVector3& posGlobal, Double_t& x_local, Double_t& y_local);
    void LocalToGlobal(TVector3& posGlobal, Double_t x_local, Double_t y_local);

    const TString& GetDetectorName() const { return fDetectorName; }

    Double_t GetEnergyLoss(const R3BTrackingParticle* particle);

    inline R3BTGeoPar* GetGeoPar() { return fGeo; }

  public:
    TString fDetectorName;
    TString fGeoParName;
    TString fDataName;

    // coordinates of the plane
    TVector3 pos0;
    TVector3 pos1;
    TVector3 pos2;

    TVector3 norm;

    // which section of the setup: 0=before target 1=target-glad 2=after glad
    EDetectorType section;

    // material + thickness
    // ??
    R3BTGeoPar* fGeo;

    // resolutions (for the chi2)
    Double_t res_x;
    Double_t res_y;
    Double_t res_t;

    // TClonesArray of hits. Holding a detector dependent structure.
    // In the Exec() function the user has to copy the hit position
    // of interest into the hit_# structures.
    std::vector<R3BHit*> hits; // not used directly by the fitter
    TClonesArray* fArrayHits;

    // Double_t hit_x; // local coordinates (on the detector plane)
    // Double_t hit_y;
    // TVector3 hit_xyz; // same in global coordinates
    // Double_t hit_time;
    // we could store a time here as well. Maybe later.

    // Double_t track_x; // local coordinates (on the detector plane)
    // Double_t track_y;
    // TVector3 track_xyz; // same in global coordinates

  public:
    ClassDef(R3BTrackingDetector, 0)
};

#endif
