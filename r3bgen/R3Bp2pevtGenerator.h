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

/*
 * R3Bp2pevtGenerator.h
 *
 *  Created on: Jul 23, 2018
 *      Author: gu92joq
 */

#ifndef R3BP2PEVTGENERATOR_H_
#define R3BP2PEVTGENERATOR_H_

#include <map>
#include <utility>

#include <TRandom3.h>
#include <TVector3.h>

#include <FairGenerator.h>

/**
 * Event generator for p2p events in inverse and direct kinematics, alpha version
 * Supply energy, nucleus A and excitation energy by separate method SetValues
 */

struct cm_values
{
    // Internal cluster
    double e_clust;
    double p_clust;
    double theta_clust;
    // Scattered particle
    double e_scat;
    double p_scat;
    double theta_scat;
    // indicates satisfactory kinematics (i.e. energy & momentum conservation)
    bool good;
    double T;
};

class R3Bp2pevtGenerator : public FairGenerator
{
  public:
    // R3Bp2pevtGenerator() : FairGenerator("R3Bp2pevtGenerator", "R3Bp2pevtGenerator") {}
    // R3Bp2pevtGenerator(const char *name, const char *title) : FairGenerator(name, title) {}
    R3Bp2pevtGenerator();
    ~R3Bp2pevtGenerator(){};

    virtual Bool_t Init();
    virtual Bool_t ReadEvent(FairPrimaryGenerator* primGen);
    void SetHeavyNucleus(int A, double M_def, double M_res); // Set Mass of heavy nucleus and Mass Defect
    void SetInverse(bool inverse);                           // Set Inverse or direct kinematics
    void SetIsotropic(bool isotropic);          // Set Isotropic or parametric cross section for proton emission
    void SetMomDistrib(double mom);             // Set Sigma of momentum distribution
    void SetExcitation(double exe);             // Set Excitation energy
    void SetLightNucleus(double ma, double mi); // Set Incoming outgoing particle mass
    void SetBeamEnergy(double E);
    void Print();

  protected:
    /**
     * Set Values for target end beam nucleus, Cross section and kinematics
     */
    double ENERGY;
    int A_BEAM;
    double MA;
    double MB;
    double Ma;
    double Mi;
    double T_LIMIT;
    double MOM_SIGMA;
    bool ISOTROPIC;
    bool INVERSE;
    TRandom3 fRandom;
    double PI = 3.14159265358979323846;
    TVector3 DREHUNG(TVector3 v1, TVector3 v2);
    double CINEMA(double, double, double);
    cm_values CENMASS(double, double, double, double, bool);
    double momentum_CM(double, double, double);
    double get_T(double, double);
    std::pair<double, double> Lorentz(double, double, double, double);
    void SetValues(double E, int A, double MOM, double exe, bool invert, bool iso);

    ClassDef(R3Bp2pevtGenerator, 1);
};

#endif
