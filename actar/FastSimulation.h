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

#ifndef GARFIELD_FAST_SIMULATION_H
#define GARFIELD_FAST_SIMULATION_H

//------------------------------------------------
// The Virtual Monte Carlo examples
// Copyright (C) 2007 - 2016 Ivana Hrivnacova
// All rights reserved.
//
// For the licensing terms see geant4_vmc/LICENSE.
// Contact: root-vmc@cern.ch
//-------------------------------------------------

/// \file Garfield/geant4/include/FastSimulation.h
/// \brief Definition of the Garfield::FastSimulation class
///
/// Garfield garfieldpp example adapted to Virtual Monte Carlo.
///
/// \date 28/10/2015
/// \author I. Hrivnacova; IPN, Orsay

#include "TG4VUserFastSimulation.h"

class GarfieldMessenger;

namespace VMC
{
    namespace Garfield
    {

        /// \ingroup Garfield
        /// \brief Special class for definition of fast simulation models
        ///
        /// \author I. Hrivnacova; IPN, Orsay

        class FastSimulation : public TG4VUserFastSimulation
        {
          public:
            FastSimulation();
            virtual ~FastSimulation();

            // methods
            virtual void Construct();

          private:
            GarfieldMessenger* fMessenger;
        };

    } // namespace Garfield
} // namespace VMC

#endif // GFLASH_FAST_SIMULATION_H
