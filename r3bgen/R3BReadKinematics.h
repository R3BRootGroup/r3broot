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

/////////////////////////////////////////////////////////////////
//*-- AUTHOR : Oleg Kiselev		O.Kiselev@gsi.de
//*-- Date: 11/2005
//*-- Last Update: 20/07/05 Hector Alvarez
// --------------------------------------------------------------
// Description:
//   Event generator kinematics reader
//
// --------------------------------------------------------------
// Comments:
//
// --------------------------------------------------------------
//
/////////////////////////////////////////////////////////////////

#ifndef READKINEM_H
#define READKINEM_H

#include "TObject.h"
#include <fstream>
#include <iostream>
#include <string>

class R3BReadKinematics : public TObject
{
  public:
    R3BReadKinematics();
    ~R3BReadKinematics();

  public:
    //    File *pFile;
    Double_t LabAngle[1000];
    Double_t LabEnergy[1000];
    Double_t T_xsec[5000];
    Double_t CrossSection_t[5000];
    // Integrated matrix (cross section)
    Double_t Icross_section[5000];

    ClassDef(R3BReadKinematics, 1) // ROOT CINT
};
#endif
