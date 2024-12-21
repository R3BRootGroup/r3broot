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

#ifndef GARFIELDG4FASTSIMULATIONMODEL_H_
#define GARFIELDG4FASTSIMULATIONMODEL_H_

//------------------------------------------------
// The Virtual Monte Carlo examples
// Copyright (C) 2007 - 2016 Ivana Hrivnacova
// All rights reserved.
//
// For the licensing terms see geant4_vmc/LICENSE.
// Contact: root-vmc@cern.ch
//-------------------------------------------------

/// \file  Garfield/geant4/include/GarfieldG4FastSimulationModel.h
/// \brief Definition of the GarfieldG4FastSimulationModel class
///
/// Garfield garfieldpp example adapted to Virtual Monte Carlo.
/// This class is imported from garfieldpp example.
/// There is disabled WriteGeometryToGDML function in order to avoid
/// requiring Geant4 optional library.
///
/// \date 28/10/2015
/// \author D. Pheiffer, CERN

/*
 * GarfieldModel.h
 *
 *  Created on: Apr 9, 2014
 *      Author: dpfeiffe
 */

#include "G4VFastSimulationModel.hh"
#include "GarfieldPhysics.h"

class G4VPhysicalVolume;

class GarfieldG4FastSimulationModel : public G4VFastSimulationModel
{
  public:
    //-------------------------
    // Constructor, destructor
    //-------------------------
    GarfieldG4FastSimulationModel(G4String, G4Region*);
    GarfieldG4FastSimulationModel(G4String);
    ~GarfieldG4FastSimulationModel();

    void SetPhysics(GarfieldPhysics* fGarfieldPhysics);
    // I.H. make this optional
    // void WriteGeometryToGDML(G4VPhysicalVolume* physicalVolume);

    virtual G4bool IsApplicable(const G4ParticleDefinition&);
    virtual G4bool ModelTrigger(const G4FastTrack&);
    virtual void DoIt(const G4FastTrack&, G4FastStep&);

  private:
    GarfieldPhysics* fGarfieldPhysics;
};

#endif /* GARFIELDMODEL_H_ */
