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

// -------------------------------------------------------------------------
// -----          based on FairAsciiGenerator header file              -----
// -----          Created 30/09/10  by M. Labiche                      -----
// -------------------------------------------------------------------------

/**  R3Bp2pGenerator.h
 Derived from FairGenerator.
**/

#ifndef FAIR_ASCIIGENERATOR_H
#define FAIR_ASCIIGENERATOR_H

#include "FairGenerator.h"

#include <fstream>

using namespace std;

class TDatabasePDG;
class FairPrimaryGenerator;

class R3Bp2pGenerator : public FairGenerator
{

  public:
    /** Default constructor without arguments should not be used. **/
    R3Bp2pGenerator();

    /** Standard constructor.
     ** @param fileName The input file name
     **/
    R3Bp2pGenerator(const char* fileName);

    R3Bp2pGenerator(const R3Bp2pGenerator&);

    R3Bp2pGenerator& operator=(const R3Bp2pGenerator&) { return *this; }

    /** Destructor. **/
    virtual ~R3Bp2pGenerator();

    /** Reads on event from the input file and pushes the tracks onto
     ** the stack. Abstract method in base class.
     ** @param primGen  pointer to the FairPrimaryGenerator
     **/
    virtual Bool_t ReadEvent(FairPrimaryGenerator* primGen);

  private:
    ifstream* fInputFile;    //! Input file stream
    const Char_t* fFileName; //! Input file Name

    ofstream* fOutputFile;

    /** Private method CloseInput. Just for convenience. Closes the
     ** input file properly. Called from destructor and from ReadEvent. **/
    void CloseInput();

    /** PDG data base */

    //  TDatabasePDG *fPDG; //!

    ClassDef(R3Bp2pGenerator, 1);
};

#endif
