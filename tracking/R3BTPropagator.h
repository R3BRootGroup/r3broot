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

#ifndef R3B_T_PROPAGATOR
#define R3B_T_PROPAGATOR

#include "TCanvas.h"
#include "TObject.h"
#include "TVector3.h"

class R3BGladFieldMap;
class FairRKPropagator;
class R3BTGeoPar;
class FairField;
class R3BTrackingParticle;
class R3BTrackingDetector;

constexpr int N_PARTICLE_INFO = 7;

class R3BTPropagator : public TObject
{
  public:
    R3BTPropagator(R3BGladFieldMap* field, Bool_t vis = kFALSE);

    virtual ~R3BTPropagator();

    inline void SetmTofGeo(R3BTGeoPar* geo) { fmTofGeo = geo; }

    Bool_t PropagateToDetector(R3BTrackingParticle* particle, R3BTrackingDetector* detector);

    Bool_t PropagateToDetectorBackward(R3BTrackingParticle* particle, R3BTrackingDetector* detector);

    Bool_t PropagateToPlane(R3BTrackingParticle* particle, const TVector3& v1, const TVector3& v2, const TVector3& v3);

    Bool_t PropagateToPlaneBackward(R3BTrackingParticle* particle,
                                    const TVector3& v1,
                                    const TVector3& v2,
                                    const TVector3& v3);

    Bool_t PropagateToPlaneRK(R3BTrackingParticle* particle,
                              const TVector3& v1,
                              const TVector3& v2,
                              const TVector3& v3);

    Bool_t PropagateToPlaneRK_eloss(R3BTrackingParticle* particle,
                                    const TVector3& point1,
                                    const TVector3& point2,
                                    const TVector3& point3,
                                    double& Tot_Tof,
                                    double step,
                                    bool eloss);

    Bool_t LineIntersectPlane(const TVector3& pos,
                              const TVector3& mom,
                              const TVector3& v1,
                              const TVector3& normal,
                              TVector3& intersect);

    static Double_t Energy_loss_in_air(double beta, double step, double charge, double mass);

    static void update_particle(R3BTrackingParticle* particle,
                                Double_t length,
                                std::array<double, N_PARTICLE_INFO>& vecOut,
                                bool eloss);

    void SetVis(Bool_t vis = kTRUE) { fVis = vis; }

  private:
    FairRKPropagator* fFairProp;

    FairField* fField;

    R3BTGeoPar* fmTofGeo;

    Bool_t fVis;

    TVector3 fPlane1[3];
    TVector3 fPlane2[3];
    TVector3 fNorm1;
    TVector3 fNorm2;

    TCanvas* fc4;

    ClassDef(R3BTPropagator, 1)
};

#endif //! R3B_T_PROPAGATOR
