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
//*-- AUTHOR : Saul Beceiro Novo        saul.beceiro@rai.usc.es
//             Hector Alvarez-Pol       hapol@fpddux.usc.es
//             Denis Bertini            D.Bertini@gsi.de
//*-- Date: 05/2009
//*-- Last Update: 05/08/09 by <D.Bertini@gsi.de>
// --------------------------------------------------------------
// Description:
//   Coulomb Disociation generators interface
//
// --------------------------------------------------------------
// Comments:
//     - 05/08/09  Adapting to R3BRoot
//     - 07/05/07 Adding setters for including emittance
//     - 03/05/07 Cleaning and moving to svn repository
//
// --------------------------------------------------------------
/////////////////////////////////////////////////////////////////

#ifndef R3BCDGenerator_h
#define R3BCDGenerator_h 1

#include "TObject.h"

#include <fstream>
#include <iostream>
#include <string>

using namespace std;

class R3BCDGenerator : public TObject
{

  private:
    Double_t px1; // momentum components for particle 1
    Double_t py1;
    Double_t pz1;
    Double_t E1;  // energy for particle 1
    Double_t px2; // momentum components for particle 2
    Double_t py2;
    Double_t pz2;
    Double_t E2; // energy for particle 1

    ifstream* inputFile; //!

  public:
    R3BCDGenerator();
    R3BCDGenerator(const R3BCDGenerator&);
    ~R3BCDGenerator();

    R3BCDGenerator& operator=(const R3BCDGenerator&) { return *this; }

    void ReadNewLine();

    Double_t GetPx1() { return px1; }
    Double_t GetPy1() { return py1; }
    Double_t GetPz1() { return pz1; }
    Double_t GetE1() { return E1; }
    Double_t GetPx2() { return px2; }
    Double_t GetPy2() { return py2; }
    Double_t GetPz2() { return pz2; }
    Double_t GetE2() { return E2; }

    void SetPx1(Double_t px) { px1 = px; }
    void SetPy1(Double_t py) { py1 = py; }
    void SetPz1(Double_t pz) { pz1 = pz; }
    void SetE1(Double_t E) { E1 = E; }
    void SetPx2(Double_t px) { px2 = px; }
    void SetPy2(Double_t py) { py2 = py; }
    void SetPz2(Double_t pz) { pz2 = pz; }
    void SetE2(Double_t E) { E2 = E; }

    ClassDef(R3BCDGenerator, 0);
};
#endif
