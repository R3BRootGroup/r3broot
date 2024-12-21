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

#pragma once

#include <FairGenerator.h>
#include <TFile.h>
#include <TString.h>
#include <TTree.h>
#include <fstream>
#include <iostream>
#include <map>
#include <sstream>
#include <string>

class FairPrimaryGenerator;

class R3BINCLRootGenerator : public FairGenerator
{
  public:
    /** Default constructor without arguments should not be used. **/
    R3BINCLRootGenerator();

    /** Standard constructor.
     ** @param fileName The input file name
     **/
    explicit R3BINCLRootGenerator(std::string fileName);
    explicit R3BINCLRootGenerator(const TString& fileName);
    explicit R3BINCLRootGenerator(const char* fileName); // for old macros

    /** Destructor. **/
    ~R3BINCLRootGenerator() override;

    /** Reads on event from the input file and pushes the tracks onto
     ** the stack. Abstract method in base class.
     ** @param primGen  pointer to the R3BPrimaryGenerator
     **/
    bool ReadEvent(FairPrimaryGenerator* primGen) override;

    /**
     ** Set vertex point
     **/
    void SetXYZ(Double32_t x = 0, Double32_t y = 0, Double32_t z = 0);

    /**
     ** Set dispersion for the vertex point
     ** Gaussian distribution for x and y
     ** Uniform distribution for z, between z - sz and z + sz
     **/
    void SetDxDyDz(Double32_t sx = 0, Double32_t sy = 0, Double32_t sz = 0);

    /**
     ** Method to simulate only fission events
     **/
    void SetOnlyFission(Bool_t Opt = true) { fOnlyFission = Opt; }

    /**
     **  Method to simulate only spallation events
     **/
    void SetOnlySpallation(Bool_t Opt = true) { fOnlySpallation = Opt; }

    /**
     **  Method to simulate only p2p-fission events
     **/
    void SetOnlyP2pFission(Bool_t Opt = true) { fOnlyP2pFission = Opt; }

  private:
    TString fFileName; // Input file name
    TFile* fInput;
    Bool_t fOnlyFission;    // True if we want to simulate only fission events
    Bool_t fOnlyP2pFission; // True if we want to simulate only p2p-fission events
    Bool_t fOnlySpallation; // True if we want to simulate only spallation events

    /** Private method RegisterIons. Goes through the input file and registers
     ** any ion needed. TODO: Should not be needed by FairRoot. **/
    void RegisterIons();

    inline Int_t GetIonPdgId(int z, int a) { return 1000000000 + 10000 * z + 10 * a; }

    Int_t fEvt;
    TTree* Tree;
    Int_t fEvtRoot;
    Short_t fParticles;
    Short_t fOrigin[160];
    Short_t fMass[160];
    Short_t fCharge[160];
    Int_t fPdgCode[160];
    Float_t fPzPrime[160];
    Float_t fThetaPrime[160];
    Float_t fPhi[160];

    Double32_t fX, fY, fZ;    // Point vertex coordinates [cm]
    Bool_t fPointVtxIsSet;    // True if point vertex is set
    Double32_t fDX, fDY, fDZ; // Point vertex coordinates [cm]
    Bool_t fBoxVtxIsSet;      // True if point vertex is set

    ClassDefOverride(R3BINCLRootGenerator, 0);
};
