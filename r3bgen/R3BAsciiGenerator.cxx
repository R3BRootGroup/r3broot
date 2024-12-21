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

#include "R3BAsciiGenerator.h"
#include "FairIon.h"
#include "FairLogger.h"
#include "FairPrimaryGenerator.h"
#include "FairRunSim.h"
#include "G4NistManager.hh"
#include "TRandom.h"
#include <boost/filesystem.hpp>
// #include <boost/iostreams/filter/bzip2.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <map>

int GetIonPdg(int z, int a) { return 1000000000 + 10 * 1000 * z + 10 * a; }

R3BAsciiGenerator::R3BAsciiGenerator()
    : fFileName()
    , fFile()
    , fBuf()
    , fInput(&fBuf)
    , fX(0.)
    , fY(0.)
    , fZ(0.)
    , fPointVtxIsSet(false)
    , fDX(0.)
    , fDY(0.)
    , fDZ(0.)
    , fBoxVtxIsSet(false)
{
}

R3BAsciiGenerator::R3BAsciiGenerator(std::string fileName)
    : fFileName(std::move(fileName))
    , fFile()
    , fBuf()
    , fInput(&fBuf)
    , fX(0.)
    , fY(0.)
    , fZ(0.)
    , fPointVtxIsSet(false)
    , fDX(0.)
    , fDY(0.)
    , fDZ(0.)
    , fBoxVtxIsSet(false)
{
    RegisterIons();
}

R3BAsciiGenerator::R3BAsciiGenerator(const TString& fileName)
    : R3BAsciiGenerator(std::string(fileName))
{
}

R3BAsciiGenerator::R3BAsciiGenerator(const char* fileName)
    : R3BAsciiGenerator(std::string(fileName))
{
}

R3BAsciiGenerator::~R3BAsciiGenerator() {}

bool R3BAsciiGenerator::ReadEvent(FairPrimaryGenerator* primGen)
{
    // Event variable to be read from file
    int eventId = -1;
    int nTracks = -1;
    // Track variables to be read from file
    int iPid = -1;
    int iZ = 0;
    int iA = 0;
    double px = 0.;
    double py = 0.;
    double pz = 0.;
    double vx = 0.;
    double vy = 0.;
    double vz = 0.;
    double ivx = 0.;
    double ivy = 0.;
    double ivz = 0.;

    OpenOrRewindFile();

    // Read event number and number of primary particles.
    if (!(fInput >> eventId >> nTracks))
    {
        // That might fail for the newline at the end of the file - that's ok - just reopen & try again
        if (fInput.eof())
        {
            return ReadEvent(primGen);
        }
        // Throw error if its somewhere in the file
        LOG(fatal) << "R3BAsciiGenerator: Could not read event header " << eventId << "\t" << nTracks;
    }
    // Ignore the other stuff that might still be on that line
    fInput.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
    LOG(debug) << "R3BAsciiGenerator: Event " << eventId << " nTracks " << nTracks;

    if (fPointVtxIsSet)
    {
        if (fBoxVtxIsSet)
        {
            vx = gRandom->Gaus(fX, fDX);
            vy = gRandom->Gaus(fY, fDY);
            vz = gRandom->Uniform(fZ - fDZ / 2.0, fZ + fDZ / 2.0);
        }
        else
        {
            vx = fX;
            vy = fY;
            vz = fZ;
        }
    }

    // Loop over tracks in the current event
    for (int itrack = 0; itrack < nTracks; itrack++)
    {
        if (!(fInput >> iPid >> iZ >> iA >> px >> py >> pz >> ivx >> ivy >> ivz))
        {
            LOG(fatal) << "R3BAsciiGenerator: Error while reading particles for event" << eventId;
        }
        // Ignore the other stuff that might still be on that line
        fInput.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

        // Ions: -1, Particles +1
        int pdg = iPid == -1 ? GetIonPdg(iZ, iA) : iPid;

        if (!fPointVtxIsSet)
        {
            vx = ivx;
            vy = ivy;
            vz = ivz;
        }

        LOG(debug) << "R3BAsciiGenerator: Adding track " << iPid << "\t" << iZ << "\t" << iA << "\t" << px << "\t" << py
                   << "" << pz << "\t" << vx << "\t" << vy << "" << vz;
        primGen->AddTrack(pdg, px, py, pz, vx, vy, vz);

    } //! tracks

    return true;
}

void R3BAsciiGenerator::RegisterIons()
{
    LOG(info) << "R3BAsciiGenerator: Looking for ions ...";

    //  Event variable to be read from file
    int eventId = -1;
    int nTracks = -1;
    // Track variables to be read from file
    int iPid = -1;
    int iZ = 0;
    int iA = 0;
    // Keep a list of ions to register
    std::map<int, FairIon*> ions;

    OpenOrRewindFile();
    while (!fInput.eof())
    {
        // Read event number and number of primary particles.
        if (!(fInput >> eventId >> nTracks))
        {
            // That might fail for the newline at the end of the file - that's ok - we're done here
            if (fInput.eof())
            {
                break;
            }
            // Throw error if its somewhere in the file
            LOG(fatal) << "R3BAsciiGenerator: Could not read event header " << eventId << "\t" << nTracks;
        }
        // Ignore the other stuff that might still be on that line
        fInput.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
        LOG(debug) << "R3BAsciiGenerator: Event " << eventId << " nTracks " << nTracks;

        for (int iTrack = 0; iTrack < nTracks; iTrack++)
        {
            if (!(fInput >> iPid >> iZ >> iA))
            {
                LOG(fatal) << "R3BAsciiGenerator: Error while reading particles for event" << eventId;
            }
            // Ignore the other stuff that might still be on that line
            fInput.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

            if (iPid == -1)
            {
                const int pdg = GetIonPdg(iZ, iA);
                if (ions.find(pdg) == ions.end())
                {
                    const double mass = G4NistManager::Instance()->GetIsotopeMass(iZ, iA) / CLHEP::GeV;
                    LOG(debug) << "R3BAsciiGenerator: New ion " << iZ << "\t" << iA << "\t" << mass;
                    ions[pdg] = new FairIon(TString::Format("Ion_%d_%d", iA, iZ), iZ, iA, iZ, 0., mass);
                }
            }
        }
    }

    for (const auto& kv : ions)
    {
        // Note: FairRoot will not register ions known to TDatabasePDG (e.g. alphas)
        FairRunSim::Instance()->AddNewIon(kv.second);
    }

    LOG(info) << "R3BAsciiGenerator: " << ions.size() << " ions registered.";
}

void R3BAsciiGenerator::SetXYZ(Double32_t x, Double32_t y, Double32_t z)
{
    fX = x;
    fY = y;
    fZ = z;
    fPointVtxIsSet = kTRUE;
}

void R3BAsciiGenerator::SetDxDyDz(Double32_t sx, Double32_t sy, Double32_t sz)
{
    fDX = sx;
    fDY = sy;
    fDZ = sz;
    fBoxVtxIsSet = kTRUE;
}

void R3BAsciiGenerator::OpenOrRewindFile()
{
    // Rewind
    if (fInput.eof())
    {
        fFile.close();
        fBuf.reset();
        fInput.clear();
    }

    if (!fFile.is_open())
    {
        LOG(info) << "R3BAsciiGenerator: Reading input file " << fFileName;
        fFile.open(fFileName);
        const std::string extension = boost::filesystem::path(fFileName).extension().string();
        if (extension == ".gz")
        {
            LOG(debug) << "R3BAsciiGenerator: Using gzip";
            fBuf.push(boost::iostreams::gzip_decompressor());
        }
        // if (extension == ".bz2")
        //{
        //    LOG(debug) << "R3BAsciiGenerator: Using bzip2";
        //    fBuf.push(boost::iostreams::bzip2_decompressor());
        //}
        fBuf.push(fFile);
    }

    // Check for input file
    if (!fFile.is_open())
    {
        LOG(fatal) << "R3BAsciiGenerator: Input file not open!";
    }
}

ClassImp(R3BAsciiGenerator);
