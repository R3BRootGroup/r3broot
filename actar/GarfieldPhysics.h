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

#ifndef GARFIELDMODELCONFIG_HH_
#define GARFIELDMODELCONFIG_HH_

//------------------------------------------------
// The Virtual Monte Carlo examples
// Copyright (C) 2007 - 2016 Ivana Hrivnacova
// All rights reserved.
//
// For the licensing terms see geant4_vmc/LICENSE.
// Contact: root-vmc@cern.ch
//-------------------------------------------------

/// \file  Garfield/include/GarfieldPhysics.h
/// \brief Definition of the GarfieldPhysics class
///
/// Garfield garfieldpp example adapted to Virtual Monte Carlo.
/// This class is imported from garfieldpp example with no modification
/// in the code.
///
/// \date 28/10/2015
/// \author D. Pheiffer, CERN

/*
 * GarfieldPhysics.hh
 *
 *  Created on: Jul 17, 2014
 *      Author: dpfeiffe
 */

#include <iostream>
#include <map>
#include <vector>

#include "AvalancheMC.hh"
#include "AvalancheMicroscopic.hh"
#include "ComponentAnalyticField.hh"
#include "GeometryRoot.hh"
#include "GeometrySimple.hh"
#include "MediumMagboltz.hh"
#include "Sensor.hh"
#include "SolidTube.hh"
#include "TrackHeed.hh"
#include "TrackSimple.hh"

typedef std::pair<double, double> EnergyRange_MeV;
typedef std::map<const std::string, EnergyRange_MeV> MapParticlesEnergy;

class GarfieldParticle
{
  public:
    GarfieldParticle(std::string particleName,
                     double ekin_eV,
                     double time,
                     double x_cm,
                     double y_cm,
                     double z_cm,
                     double dx,
                     double dy,
                     double dz)
        : fParticleName(particleName)
        , fEkin_MeV(ekin_eV / 1000000)
        , fTime(time)
        , fx_mm(10 * x_cm)
        , fy_mm(10 * y_cm)
        , fz_mm(10 * z_cm)
        , fdx(dx)
        , fdy(dy)
        , fdz(dz)
    {
    }
    ~GarfieldParticle(){};

    std::string getParticleName() { return fParticleName; }
    double getX_mm() { return fx_mm; }
    double getY_mm() { return fy_mm; }
    double getZ_mm() { return fz_mm; }
    double getEkin_MeV() { return fEkin_MeV; }
    double getTime() { return fTime; }
    double getDX() { return fdx; }
    double getDY() { return fdy; }
    double getDZ() { return fdz; }

  private:
    std::string fParticleName;
    double fEkin_MeV, fTime, fx_mm, fy_mm, fz_mm, fdx, fdy, fdz;
};

class GarfieldPhysics
{
  public:
    static GarfieldPhysics* GetInstance();
    static void Dispose();

    void InitializePhysics();
    void CreateGeometry();

    void DoIt(std::string particleName,
              double ekin_MeV,
              double time,
              double x_cm,
              double y_cm,
              double z_cm,
              double dx,
              double dy,
              double dz);

    void AddParticleName(const std::string particleName, double ekin_min_keV, double ekin_max_keV);
    bool FindParticleName(const std::string name);
    bool FindParticleNameEnergy(std::string name, double ekin_keV);
    void SetIonizationModel(std::string model, bool useDefaults = true);
    std::string GetIonizationModel();
    std::vector<GarfieldParticle*>* GetSecondaryParticles();
    void DeleteSecondaryParticles();
    inline void EnableCreateSecondariesInGeant4(bool flag) { createSecondariesInGeant4 = flag; };
    inline bool GetCreateSecondariesInGeant4() { return createSecondariesInGeant4; };
    inline double GetEnergyDeposit_MeV() { return fEnergyDeposit / 1000000; };
    inline double GetAvalancheSize() { return fAvalancheSize; };
    inline double GetGain() { return fGain; };
    inline void Clear()
    {
        fEnergyDeposit = 0;
        fAvalancheSize = 0;
        fGain = 0;
        nsum = 0;
    }

  private:
    GarfieldPhysics();
    ~GarfieldPhysics();

    std::string fIonizationModel;

    static GarfieldPhysics* fGarfieldPhysics;
    MapParticlesEnergy* fMapParticlesEnergy;
    TGeoManager* fGeoManager;
    Garfield::MediumMagboltz* fMediumMagboltz;
    Garfield::Sensor* fSensor;
    Garfield::AvalancheMC* fDrift;
    Garfield::AvalancheMicroscopic* fAvalanche;
    Garfield::TrackHeed* fTrackHeed;
    Garfield::GeometryRoot* fGeometryRoot;
    Garfield::GeometrySimple* fGeometrySimple;
    Garfield::ComponentAnalyticField* fComponentAnalyticField;
    Garfield::SolidTube* fTube;

    std::vector<GarfieldParticle*>* fSecondaryParticles;

    bool createSecondariesInGeant4;
    double fEnergyDeposit;
    double fAvalancheSize;
    double fGain;
    int nsum;

    ClassDef(GarfieldPhysics, 1)
};
#endif /* GARFIELDMODELCONFIG_HH_ */
