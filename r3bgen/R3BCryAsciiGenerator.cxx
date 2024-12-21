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

// ------------------------------------------------------------------------
//                 R3BCryAsciiGenerator source file                   -----
// ------------------------------------------------------------------------
// ----------  Use this generator with an appropiate ascii file  ----------
// ----------     generated with CRY with the "pdg" option       ----------
// ------------------------------------------------------------------------

#include "R3BCryAsciiGenerator.h"

#include "FairLogger.h"
#include "FairPrimaryGenerator.h"
#include "FairRunSim.h"

#include "TDatabasePDG.h"
#include "TRandom.h"
#include "TString.h"

using namespace std;

R3BCryAsciiGenerator::R3BCryAsciiGenerator()
    : fFileName()
    , fTopDist(0.0)
{
}

R3BCryAsciiGenerator::R3BCryAsciiGenerator(std::string fileName)
    : fFileName(std::move(fileName))
    , fTopDist(0.0)
{
    infile.open(fileName);
    if (!infile.is_open())
        LOG(error) << "R3BCryAsciiGenerator: Cannot open input file.";
}

R3BCryAsciiGenerator::R3BCryAsciiGenerator(const TString& fileName)
    : R3BCryAsciiGenerator(std::string(fileName))
{
}

R3BCryAsciiGenerator::R3BCryAsciiGenerator(const char* fileName)
    : R3BCryAsciiGenerator(std::string(fileName))
{
}

R3BCryAsciiGenerator::~R3BCryAsciiGenerator() { CloseInput(); }

bool R3BCryAsciiGenerator::ReadEvent(FairPrimaryGenerator* primGen)
{

    streampos prevpos;

    // Define event variable to be read from file
    Int_t eventId = 0; // event id
    Int_t nTracks = 0; // number of tracks in the event
    Int_t fPrevEventId = 0;

    // Define track variables
    Int_t pdgType = 0;  // PDG code
    Double_t KEn = 0.;  // kinetic energy read in MeV
    Double_t cosA = 0.; // direction cosines
    Double_t cosB = 0.;
    Double_t cosC = 0.;
    Double_t vx = 0.; // vertex componets
    Double_t vy = 0.;
    Double_t vz = 0.;

    // Intermediate varibles
    Double_t p = 0.;  // total momentum
    Double_t px = 0.; // momentum componets
    Double_t py = 0.;
    Double_t pz = 0.;
    Double_t mass = 0.; // particle mass read in GeV

    Bool_t read = true;

    infile >> eventId >> nTracks >> pdgType >> KEn >> vx >> vy >> vz >> cosA >> cosB >> cosC;

    if (infile.eof())
    {
        LOG(info) << "R3BCryAsciiGenerator: End of input file reached ";
        CloseInput();
        return kFALSE;
    }

    while (read && !infile.eof())
    {
        // We extract the particle mass from its object
        TParticlePDG* particle = TDatabasePDG::Instance()->GetParticle(pdgType);
        mass = particle->Mass();

        // Calculations
        vx = vx * 1.; // input should be in cm
        vy = fTopDist;
        vz = vy * 1.;
        KEn = KEn / 1000.0; // energy should be in GeV
        p = TMath::Sqrt((KEn + mass) * (KEn + mass) - mass * mass);
        px = p * cosA;
        py = p * cosB;
        pz = p * cosC;

        // We create the track
        primGen->AddTrack(pdgType, px, py, pz, vx, 10., vy);

        fPrevEventId = eventId;

        // Read the next line in the file, check if it belongs to the same event, and go back to previous line
        prevpos = infile.tellg();
        infile >> eventId >> nTracks >> pdgType >> KEn >> vx >> vy >> vz >> cosA >> cosB >> cosC;
        if (fPrevEventId != eventId)
        {
            read = false;
            infile.seekg(prevpos);
        }
    }

    return true;
}

// -----   Private method CloseInput   ------------------------------------
void R3BCryAsciiGenerator::CloseInput()
{
    if (infile)
    {
        if (infile.is_open())
        {
            LOG(info) << "R3BCryAsciiGenerator: Closing input file ";
            infile.close();
        }
    }
}

ClassImp(R3BCryAsciiGenerator)
