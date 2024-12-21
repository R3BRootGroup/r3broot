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

#ifndef R3BNEUTRONTRACKER_H
#define R3BNEUTRONTRACKER_H 1

#include "FairTask.h"
#include "R3BLandDigi.h"
#include "R3BLandDigiPar.h"
#include "R3BLandFirstHits.h"
#include "TLorentzVector.h"
#include <map>
#include <string>

class TClonesArray;
class TObjectArray;
class TH1F;
class TH2F;
class R3BNeutronTrack;

struct PM_HIT
{
    Double_t paddle, tdcL, tdcR, tdc, qdcL, qdcR, qdc, x, y, z, x0, y0, z0, t0;
};
struct PRIM_PART
{
    Double_t pdg, px, py, pz, pt, p, x, y, z, t, A, M;
};
struct NEUT1_HIT
{
    Double_t x, y, z, t, px, py, pz, p;
};
struct NEUT2_HIT
{
    Double_t x, y, z, t, px, py, pz, p;
};
struct CLUSTER
{
    Double_t xStart, yStart, zStart, xEnd, yEnd, zEnd, tStart, tStop, e, size;
};

class R3BNeutronTracker : public FairTask
{

  public:
    /** Default constructor **/
    R3BNeutronTracker();

    /** Destructor **/
    ~R3BNeutronTracker();

    /** Virtual method Init **/
    virtual InitStatus Init();

    /** Virtual method Exec **/
    virtual void Exec(Option_t* opt);

    virtual void Finish();
    virtual void Reset();

    void UseBeam(Double_t _beam_energy, Double_t _beam_beta);

  private:
    TLorentzVector fPosIn, fPosOut; //!  position
    TLorentzVector fMomIn, fMomOut; //!  momentum
    Double32_t fTime;               //!  time

    R3BNeutronTrack* AddHit(TVector3 pos_in, TVector3 pos_out, TVector3 momOut, Double_t time);

  protected:
    // Reused structure from previous
    // Algorithm
    TClonesArray* fLandPoints;
    TClonesArray* fLandMCTrack;
    TClonesArray* fLandDigi;
    TClonesArray* fLandFirstHits;

    // Parameter class
    R3BLandDigiPar* fLandDigiPar;

    // Created structure from this
    // Algorithm
    TClonesArray* fNeutronTracks;

    // Add here control hist for tracker
    TH1F* hNeutmult;
    TH1F* hErel1;
    TH1F* hErel2;
    TH1F* hErel3;
    TH1F* hErel4;

    TH1F* hMinv;
    TH1F* hMinv0;
    TH1F* hMinv1;
    TH1F* hMinv2;
    TH1F* hDeltaPx1;
    TH1F* hDeltaPy1;
    TH1F* hDeltaPz1;
    TH1F* hDeltaPx2;
    TH1F* hDeltaPy2;
    TH1F* hDeltaPz2;
    TH1F* hClusterSize;
    TH1F* hClusterEnergy;
    TH1F* hClusters;
    TH1F* hClusters1;
    TH1F* hClusters2;
    TH1F* hHits;
    TH2F* hClusterNo_vs_Size;

    TH1F* hDeltaX;
    TH1F* hDeltaY;
    TH1F* hDeltaZ;
    TH1F* hDeltaT;

    TH1F* hFirstHitZ;

    TH1F* hDeltaP1;
    TH1F* hDeltaP2;
    TH1F* hDeltaP3;
    TH1F* hDeltaP4;
    TH1F* hDeltaP5;
    TH1F* hDeltaP6;

    TH1F* hDelta;

    Int_t npaddles;
    Int_t nplanes;
    Int_t printing;

    Double_t beamEnergy;
    Double_t beamBeta;
    Double_t plength;  // half length of paddle
    Double_t att;      // light attenuation factor [1/cm]
    Double_t amu;      // atomic mass unit in MeV/c**2
    Double_t mNeutron; // mass of neutron in amu
    Double_t c;
    Double_t cMedia;    // speed of light in material in cm/ns
    Double_t calFactor; // calibration factor energy of LAND paddles
    Int_t eventNo;
    PM_HIT* PM_hit;
    PRIM_PART* PRIM_part;
    PRIM_PART* PRIM_frag;
    PRIM_PART* PRIM_gamma;
    NEUT1_HIT* NEUT1_hit;
    NEUT2_HIT* NEUT2_hit;
    CLUSTER* Cluster;

  private:
    virtual void SetParContainers();

    ClassDef(R3BNeutronTracker, 1);
};

#endif
