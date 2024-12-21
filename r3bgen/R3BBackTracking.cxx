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

#include "R3BBackTracking.h"
#include "TString.h"
#include "TSystem.h"
#include <iostream>

using namespace std;

R3BBackTracking::R3BBackTracking()
    : px1(0.)
    , py1(0.)
    , pz1(0.)
    , pprim(0.)
    , x1(0.)
    , y1(0.)
    , z1(0.)
    , x0(0.)
    , y0(0.)
    , z0(0.)
    , inputFile(NULL)
{
    //
    // Constructor
    //
    TString WorkDir = gSystem->Getenv("VMCWORKDIR");
    TString kine = WorkDir + "/input/kine/reconstructionSim.dat";
    inputFile = new ifstream(kine.Data());

    // cout << "Entrada en fichero"<< *inputFile.tellg() << endl;
    // while(!eof) *inputFile.peek();
    // cout << "Salida de fichero"<<*inputFile.tellg()<< endl;
}

R3BBackTracking::R3BBackTracking(const R3BBackTracking& right)
    : px1(right.px1)
    , py1(right.py1)
    , pz1(right.pz1)
    , pprim(right.pprim)
    , x1(right.x1)
    , y1(right.y1)
    , z1(right.z1)
    , x0(right.x0)
    , y0(right.y0)
    , z0(right.z0)
    , inputFile(right.inputFile)
{
}

R3BBackTracking::~R3BBackTracking()
{
    //
    // Destructor
    //

    delete inputFile;
}

void R3BBackTracking::ReadNewLine()
{
    //
    // Reads a new line from the CD input file
    //

    Char_t* myChar = new Char_t[256];

    //  cout << inputFile << "" << *inputFile << endl;
    // cout << "Estoy en el caracter del fichero: "<< (*inputFile).tellg() << endl;
    //(*inputFile).seekg((int)(12100000*0.00001));
    // cout << "y me muevo a: "<<(*inputFile).tellg() << endl;

    // cout << "mirame, tonto!" << endl;
    // cout << &myChar << "          "<< myChar << "             " <<*myChar << endl;

    (*inputFile).getline(myChar, 256);

    // cout << myChar << " " <<*myChar << endl;

    *inputFile >> px1 >> py1 >> pz1 >> x1 >> y1 >> z1 >> x0 >> y0 >> z0 >> pprim;
    // cout << "cuchufletas las justas!" << endl;
}

ClassImp(R3BBackTracking)
