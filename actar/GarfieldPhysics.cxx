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

//------------------------------------------------
// The Virtual Monte Carlo examples
// Copyright (C) 2007 - 2016 Ivana Hrivnacova
// All rights reserved.
//
// For the licensing terms see geant4_vmc/LICENSE.
// Contact: root-vmc@cern.ch
//-------------------------------------------------

/// \file  Garfield/src/GarfieldPhysics.cxx
/// \brief Definition of the GarfieldPhysics class
///
/// Garfield garfieldpp example adapted to Virtual Monte Carlo.
/// This class is imported from garfieldpp example with no modification
/// in the code.
///
/// \date 28/10/2015
/// \author D. Pheiffer, CERN

#include "GarfieldPhysics.h"

#include "TGeoBBox.h"
#include "TGeoManager.h"
#include "TGeoTube.h"

#include "FairLogger.h"

double rTube = 0.;
double lTube = 0.;

// I.H.
// #include "TGDMLParse.h"
// #include "G4SystemOfUnits.hh"

GarfieldPhysics* GarfieldPhysics::fGarfieldPhysics = 0;
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

GarfieldPhysics* GarfieldPhysics::GetInstance()
{
    if (!fGarfieldPhysics)
    {
        fGarfieldPhysics = new GarfieldPhysics();
    }
    return fGarfieldPhysics;
}

void GarfieldPhysics::Dispose()
{
    delete fGarfieldPhysics;
    fGarfieldPhysics = 0;
}

GarfieldPhysics::GarfieldPhysics()
{
    fMapParticlesEnergy = new MapParticlesEnergy();
    fSecondaryParticles = new std::vector<GarfieldParticle*>();
    fMediumMagboltz = 0;
    fSensor = 0;
    fAvalanche = 0;
    fDrift = 0;
    fComponentAnalyticField = 0;
    fTrackHeed = 0;
    fGeoManager = 0;
    fGeometryRoot = 0;
    fGeometrySimple = 0;
    fTube = 0;
    createSecondariesInGeant4 = false;
    fIonizationModel = "Geant4";
    Clear();
}

GarfieldPhysics::~GarfieldPhysics()
{
    delete fMapParticlesEnergy;
    DeleteSecondaryParticles();
    delete fSecondaryParticles;
    delete fMediumMagboltz;
    delete fSensor;
    delete fAvalanche;
    delete fDrift;
    delete fComponentAnalyticField;
    delete fTrackHeed;
    delete fGeometryRoot;
    delete fGeometrySimple;

    std::cout << "Deconstructor GarfieldPhysics" << std::endl;
}

std::string GarfieldPhysics::GetIonizationModel() { return fIonizationModel; }

void GarfieldPhysics::SetIonizationModel(std::string model, bool useDefaults)
{
    if (model != "Geant4" && model != "Heed")
    {

        std::cout << "Unknown ionization model " << model << std::endl;
        std::cout << "Using Geant4 as default model!" << std::endl;
        model = "Geant4";
    }
    fIonizationModel = model;

    if (fIonizationModel == "Geant4")
    {
        if (useDefaults == true)
        {
            // Particle types and energies for which the G4FastSimulationModel with Garfield++ is valid
            this->AddParticleName("e-", 1e-6, 1e-3);
            this->AddParticleName("gamma", 1e-6, 1e+8);
        }
    }
    else if (fIonizationModel == "Heed")
    {
        if (useDefaults == true)
        {
            // Particle types and energies for which the G4FastSimulationModel with Garfield++ is valid
            this->AddParticleName("gamma", 1e-6, 1e+8);
            this->AddParticleName("e-", 6e-2, 1e+7);
            this->AddParticleName("e+", 6e-2, 1e+7);
            this->AddParticleName("mu-", 1e+1, 1e+8);
            this->AddParticleName("mu+", 1e+1, 1e+8);
            this->AddParticleName("pi-", 2e+1, 1e+8);
            this->AddParticleName("pi+", 2e+1, 1e+8);
            this->AddParticleName("kaon-", 1e+1, 1e+8);
            this->AddParticleName("kaon+", 1e+1, 1e+8);
            this->AddParticleName("proton", 9.e+1, 1e+8);
            this->AddParticleName("anti_proton", 9.e+1, 1e+8);
            this->AddParticleName("deuteron", 2.e+2, 1e+8);
            this->AddParticleName("alpha", 4.e+2, 1e+8);
        }
    }
}

void GarfieldPhysics::AddParticleName(const std::string particleName, double ekin_min_MeV, double ekin_max_MeV)
{
    if (ekin_min_MeV >= ekin_max_MeV)
    {
        std::cout << "Ekin_min=" << ekin_min_MeV << " keV is larger than Ekin_max=" << ekin_max_MeV << " keV"
                  << std::endl;
        return;
    }
    std::cout << "Garfield model is applicable for G4Particle " << particleName << " between " << ekin_min_MeV
              << " keV and " << ekin_max_MeV << " keV" << std::endl;
    fMapParticlesEnergy->insert(std::make_pair(particleName, std::make_pair(ekin_min_MeV, ekin_max_MeV)));
}

bool GarfieldPhysics::FindParticleName(std::string name)
{
    MapParticlesEnergy::iterator it;
    it = fMapParticlesEnergy->find(name);
    if (it != fMapParticlesEnergy->end())
    {
        return true;
    }
    return false;
}

bool GarfieldPhysics::FindParticleNameEnergy(std::string name, double ekin_MeV)
{
    MapParticlesEnergy::iterator it;
    it = fMapParticlesEnergy->find(name);
    if (it != fMapParticlesEnergy->end())
    {
        EnergyRange_MeV range = it->second;
        if (range.first <= ekin_MeV && range.second >= ekin_MeV)
        {
            return true;
        }
    }
    return false;
}

void GarfieldPhysics::InitializePhysics()
{

    fMediumMagboltz = new Garfield::MediumMagboltz();

    fMediumMagboltz->SetComposition("ar", 70., "co2", 30.);
    fMediumMagboltz->SetTemperature(293.15);
    fMediumMagboltz->SetPressure(760.);
    fMediumMagboltz->EnableDebugging();
    fMediumMagboltz->Initialise(true);
    fMediumMagboltz->DisableDebugging();
    // Set the Penning transfer efficiency.
    const double rPenning = 0.57;
    const double lambdaPenning = 0.;
    fMediumMagboltz->EnablePenningTransfer(rPenning, lambdaPenning, "ar");
    fMediumMagboltz->LoadGasFile("ar_70_co2_30_1000mbar.gas");

    fSensor = new Garfield::Sensor();
    fDrift = new Garfield::AvalancheMC();
    fAvalanche = new Garfield::AvalancheMicroscopic();
    fComponentAnalyticField = new Garfield::ComponentAnalyticField();

    // CreateGeometry();

    fDrift->SetSensor(fSensor);
    fAvalanche->SetSensor(fSensor);

    fTrackHeed = new Garfield::TrackHeed();
    fTrackHeed->EnableDebugging();
    fTrackHeed->SetSensor(fSensor);

    fTrackHeed->EnableDeltaElectronTransport();

    Garfield::GeometryRoot* geo = new Garfield::GeometryRoot();
    geo->SetGeometry(gGeoManager);
    geo->SetMedium("ArCO2", fMediumMagboltz);
    fComponentAnalyticField->SetGeometry(geo);

    TGeoTube* chamberInner = (TGeoTube*)gGeoManager->GetVolume("Actar_chamber_inner")->GetShape();
    rTube = chamberInner->GetRmax();
    lTube = chamberInner->GetDZ();
    //    fComponentAnalyticField->AddWire(0., 0., 2 * 25.e-4, 1000., "w");
    //    fComponentAnalyticField->AddTube(chamberInner->GetRmax(), 0., 0, "t");
    fComponentAnalyticField->AddPlaneX(-1. * chamberInner->GetDZ(), 0, "1");
    fComponentAnalyticField->AddPlaneX(chamberInner->GetDZ(), 1000., "2");

    fSensor->SetArea(-1. * chamberInner->GetDX(),
                     -1. * chamberInner->GetDY(),
                     -1. * chamberInner->GetDZ(),
                     chamberInner->GetDX(),
                     chamberInner->GetDY(),
                     chamberInner->GetDZ());
    fSensor->AddComponent(fComponentAnalyticField);
}

void GarfieldPhysics::CreateGeometry()
{
    /*	// Wire radius [cm]
        const double rWire = 25.e-4;
        // Outer radius of the tube [cm]
        const double rTube = 1.451;
        // Half-length of the tube [cm]
        const double lTube = 10.;

        fGeometrySimple = new Garfield::GeometrySimple();
        // Make a tube (centered at the origin, inner radius: 0, outer radius: rTube).
        fTube = new Garfield::SolidTube(0., 0., 0, rWire, rTube, lTube);
        // Add the solid to the geometry, together with the medium inside.
        fGeometrySimple->AddSolid(fTube, fMediumMagboltz);
        fComponentAnalyticField->SetGeometry(fGeometrySimple);

        // Voltages
        const double vWire = 1000.;
        const double vTube = 0.;
        // Add the wire in the center.
        fComponentAnalyticField->AddWire(0., 0., 2 * rWire, vWire, "w");
        // Add the tube.
        fComponentAnalyticField->AddTube(rTube, vTube, 0, "t");

        fSensor->AddComponent(fComponentAnalyticField);
    */
}

void GarfieldPhysics::DoIt(std::string particleName,
                           double ekin_MeV,
                           double time,
                           double x_cm,
                           double y_cm,
                           double z_cm,
                           double dx,
                           double dy,
                           double dz)
{

    DeleteSecondaryParticles();

    // Wire radius [cm]
    const double rWire = 25.e-4;
    // Outer radius of the tube [cm]
    // const double rTube = 1.45;
    // Half-length of the tube [cm]
    // const double lTube = 10.;

    double eKin_eV = ekin_MeV * 1e+6;

    double xc = 0., yc = 0., zc = 0., tc = 0.;
    // Number of electrons produced in a collision
    int nc = 0;
    // Energy loss in a collision
    double ec = 0.;
    // Dummy variable (not used at present)
    double extra = 0.;
    fEnergyDeposit = 0;
    // fAvalancheSize = 0;
    // fGain = 0;

    if (fIonizationModel != "Heed" || particleName == "gamma")
    {
        if (particleName == "gamma")
        {
            fTrackHeed->TransportPhoton(x_cm, y_cm, z_cm, time, eKin_eV, dx, dy, dz, nc);
        }
        else
        {
            fTrackHeed->TransportDeltaElectron(x_cm, y_cm, z_cm, time, eKin_eV, dx, dy, dz, nc);
            fEnergyDeposit = eKin_eV;
        }

        for (int cl = 0; cl < nc; cl++)
        {
            double xe, ye, ze, te;
            double ee, dxe, dye, dze;
            fTrackHeed->GetElectron(cl, xe, ye, ze, te, ee, dxe, dye, dze);
            if (ze < lTube && ze > -lTube && sqrt(xe * xe + ye * ye) < rTube)
            {
                nsum++;
                if (particleName == "gamma")
                {
                    fEnergyDeposit += fTrackHeed->GetW();
                }
                // TO DO: add this histo in MCApplication
                // analysisManager->FillH3(1, ze * 10, xe * 10, ye * 10);
                if (createSecondariesInGeant4)
                {
                    double newTime = te;
                    if (newTime < time)
                    {
                        newTime += time;
                    }
                    fSecondaryParticles->push_back(new GarfieldParticle("e-", ee, newTime, xe, ye, ze, dxe, dye, dze));
                }

                fDrift->DriftElectron(xe, ye, ze, te);

                double xe1, ye1, ze1, te1;
                double xe2, ye2, ze2, te2;

                int status;
                fDrift->GetElectronEndpoint(0, xe1, ye1, ze1, te1, xe2, ye2, ze2, te2, status);

                if (0 < xe2 && xe2 < rWire)
                {
                    xe2 += 2 * rWire;
                }
                if (0 > xe2 && xe2 > -rWire)
                {
                    xe2 += -2 * rWire;
                }
                if (0 < ye2 && ye2 < rWire)
                {
                    ye2 += 2 * rWire;
                }
                if (0 > ye2 && ye2 > -rWire)
                {
                    ye2 += -2 * rWire;
                }

                double e2 = 0.1;
                fAvalanche->AvalancheElectron(xe2, ye2, ze2, te2, e2, 0, 0, 0);

                int ne = 0, ni = 0;
                fAvalanche->GetAvalancheSize(ne, ni);
                fAvalancheSize += ne;
            }
        }
    }
    else
    {
        fTrackHeed->SetParticle(particleName);
        fTrackHeed->SetKineticEnergy(eKin_eV);
        fTrackHeed->NewTrack(x_cm, y_cm, z_cm, time, dx, dy, dz);

        while (fTrackHeed->GetCluster(xc, yc, zc, tc, nc, ec, extra))
        {
            if (zc < lTube && zc > -lTube && sqrt(xc * xc + yc * yc) < rTube)
            {
                nsum += nc;
                fEnergyDeposit += ec;
                for (int cl = 0; cl < nc; cl++)
                {
                    double xe, ye, ze, te;
                    double ee, dxe, dye, dze;
                    fTrackHeed->GetElectron(cl, xe, ye, ze, te, ee, dxe, dye, dze);
                    if (ze < lTube && ze > -lTube && sqrt(xe * xe + ye * ye) < rTube)
                    {
                        // TO DO: add this histo in MCApplication
                        // analysisManager->FillH3(1, ze * 10, xe * 10, ye * 10);
                        if (createSecondariesInGeant4)
                        {
                            double newTime = te;
                            if (newTime < time)
                            {
                                newTime += time;
                            }
                            fSecondaryParticles->push_back(
                                new GarfieldParticle("e-", ee, newTime, xe, ye, ze, dxe, dye, dze));
                        }

                        fDrift->DriftElectron(xe, ye, ze, te);

                        double xe1, ye1, ze1, te1;
                        double xe2, ye2, ze2, te2;

                        int status;
                        fDrift->GetElectronEndpoint(0, xe1, ye1, ze1, te1, xe2, ye2, ze2, te2, status);

                        if (0 < xe2 && xe2 < rWire)
                        {
                            xe2 += 2 * rWire;
                        }
                        if (0 > xe2 && xe2 > -rWire)
                        {
                            xe2 += -2 * rWire;
                        }
                        if (0 < ye2 && ye2 < rWire)
                        {
                            ye2 += 2 * rWire;
                        }
                        if (0 > ye2 && ye2 > -rWire)
                        {
                            ye2 += -2 * rWire;
                        }

                        double e2 = 0.1;
                        fAvalanche->AvalancheElectron(xe2, ye2, ze2, te2, e2, 0, 0, 0);

                        int ne = 0, ni = 0;
                        fAvalanche->GetAvalancheSize(ne, ni);
                        fAvalancheSize += ne;
                    }
                }
            }
        }
    }
    if (nsum > 0)
    {
        fGain = fAvalancheSize / nsum;
    }
}

std::vector<GarfieldParticle*>* GarfieldPhysics::GetSecondaryParticles() { return fSecondaryParticles; }

void GarfieldPhysics::DeleteSecondaryParticles()
{
    if (!fSecondaryParticles->empty())
    {
        fSecondaryParticles->erase(fSecondaryParticles->begin(), fSecondaryParticles->end());
    }
}

ClassImp(GarfieldPhysics)
