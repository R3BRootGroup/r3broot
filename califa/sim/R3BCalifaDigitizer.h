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

#ifndef R3BCALIFADIGITIZER_H
#define R3BCALIFADIGITIZER_H 1

#include <FairTask.h>
#include <R3BCalifa.h>
#include <R3BCalifaCrystalCalData.h>
#include <R3BCalifaCrystalPars4Sim.h>
#include <R3BCalifaPoint.h>
#include <R3BIOConnector.h>
#include <Rtypes.h>
#include <string>

class TClonesArray;

constexpr int barrelCrystals = 1952;

class R3BCalifaDigitizer : public FairTask
{
  public:
    /** Standard contructor **/
    R3BCalifaDigitizer();

    ~R3BCalifaDigitizer();

    /** Virtual method Init **/
    InitStatus Init() override;

    /** Virtual method Exec **/
    void Exec(Option_t*) override;

    void SetParContainers() override;

    /** Public method SetRealConfig
     **
     ** Defines the REAL experimental resolution and Thresholds of the CsI(Tl)
     *Crystals
     *@param isRealSet  Bool parameter used to set the experimental Resolution and
     *Thresholds
     **/
    void SetRealConfig(bool isRealSet = false);

    /** Public method SetExpEnergyRes
     **
     ** Sets the experimental energy resolution of the CsI crystals
     **/
    void SetExpGammaEnergyRes(double crysRes);
    void SetExpProtonEnergyRes(double crysRes);

    /** Public method SetFullProtonRange
     **
     ** Set up full CALIFA in proton range
     **/
    void SetFullProtonRange() { fFullProtonRange = true; }

    /** Public method SetComponentRes
     **
     ** Defines the experimental resolution of the CsI(Tl) components.
     *@param componentRes  Double parameter used to set the experimental
     *resolution in MeV
     **/
    void SetComponentRes(double componentRes);

    /** Public method SetDetectionThreshold
     **
     ** Defines the minimum energy requested in a crystal to be included as a
     *CrystalCal
     *@param thresholdEne  Double parameter used to set the threshold (in GeV)
     **/
    void SetDetectionThreshold(double thresholdEne);

    /** Public method SetNonUniformity
     **
     ** Defines the fNonUniformity parameter in % deviation from the central value
     *@param nonU  Double parameter setting the maximum non-uniformity allowed
     **/
    void SetNonUniformity(double nonU);

  private:
    void SetParameter();
    void FillRealConfig(int nbcry);
    void FillIdealConfig(int nbcry);
    void FillIdealConfigWithProtonRange(int nbcry);

    // Input array
    TClonesArray* fCalifaPointDataCA = nullptr;

    // Output array to new data level
    // R3B::OutputConnector<std::vector<R3BCalifaCrystalCalData>> output_data{ "CalifaCrystalCalData" };

    //  Array to store output data
    TClonesArray* fCalifaCryCalDataCA = nullptr;

    double fNonUniformity = 1.;      // Experimental non-uniformity parameter
    double fResolution = 0.;         // Experimental resolution
    double fGammaResolution = 0.;    // Experimental resolution @ 1 MeV
    double fProtonResolution = 0.;   // Experimental resolution for protons
    double fComponentRes = 0.;       // Experimental resolution for Nf and Ns
    double fThreshold = 0.;          // Minimum energy requested to create a Cal
    double fGammaSaturation = 31.;   // MeV
    double fProtonSaturation = 322.; // MeV
    bool fRealConfig = false;        // Real Configuration in CALIFA
    bool fFullProtonRange = false;   // Full proton range configuration
    int fNumberOfParams = 0;
    int fNumCrystals = 2544; // Real number of crystals, since Feb. 2024

    R3BCalifaCrystalPars4Sim* fSim_Par = nullptr; // Parameter Container for a Realistic Simulation

    /** Private method NUSmearing
     **
     ** Smears the energy according to some non-uniformity distribution
     ** Very simple scheme where the NU is introduced as a flat random
     ** distribution with limits fNonUniformity (%) of the energy value.
     **/
    double NUSmearing(double inputEnergy);

    /** Private method ExpResSmearing
     **
     ** Smears the energy according to some Experimental Resolution distribution
     ** Very simple scheme where the Experimental Resolution
     ** is introduced as a gaus random distribution with a width given by the
     ** parameter fCrystalResolution(in % @ MeV). Scales according to 1/sqrt(E)
     **/
    double ExpResSmearing(double inputEnergy);

    /** Private method CompSmearing
     **
     ** Smears the CsI(Tl) components Ns and Nf
     **/
    double CompSmearing(double inputComponent);

    /**
     ** Private method AddCrystalCal
     **/
    R3BCalifaCrystalCalData* AddCrystalCal(int crysid,
                                           double energy,
                                           double Nf,
                                           double Ns,
                                           uint64_t time,
                                           double tot_energy = 0.);

  public:
    ClassDefOverride(R3BCalifaDigitizer, 1);
};

#endif /* R3BCALIFADIGITIZER_H */
