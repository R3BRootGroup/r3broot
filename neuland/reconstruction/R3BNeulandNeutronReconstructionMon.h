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

#ifndef R3BNEULANDNEUTRONRECONSTRUCTIONMON_H
#define R3BNEULANDNEUTRONRECONSTRUCTIONMON_H

/** Neuland Neutron Reconstruction Evaluation
 * @author Jan Mayer
 *
 * Work in progress
 */

#include "FairTask.h"
#include <map>

class TClonesArray;
class TH1D;
class TH2D;

class R3BNeulandNeutronReconstructionMon : public FairTask
{
  public:
    R3BNeulandNeutronReconstructionMon(const TString input = "NeulandNeutrons",
                                       const TString output = "NeulandNeutronReconstructionMon");
    ~R3BNeulandNeutronReconstructionMon();

  private:
    // No copy and no move is allowed (Rule of three/five)
    R3BNeulandNeutronReconstructionMon(const R3BNeulandNeutronReconstructionMon&);            // copy constructor
    R3BNeulandNeutronReconstructionMon(R3BNeulandNeutronReconstructionMon&&);                 // move constructor
    R3BNeulandNeutronReconstructionMon& operator=(const R3BNeulandNeutronReconstructionMon&); // copy assignment
    R3BNeulandNeutronReconstructionMon& operator=(R3BNeulandNeutronReconstructionMon&&);      // move assignment

  protected:
    InitStatus Init() override;
    void Finish() override;

  public:
    void Exec(Option_t*) override;

  private:
    TString fInput;
    TString fOutput;

    TH1D* fhCountN;
    TH1D* fhCountNdiff;
    TH1D* fhScore;
    TH1D* fhEdiff;
    TH1D* fhErel;
    TH1D* fhErelMC;
    TH2D* fhErelVSnNreco;
    TH2D* fhErelVSnNrecoNPNIPs;
    TH2D* fhNreacNreco;
    std::map<Int_t, TH1D*> fhmErelnReco;

    TClonesArray* fPrimaryNeutronInteractionPoints;
    TClonesArray* fReconstructedNeutrons;
    TClonesArray* fMCTracks;

  public:
    ClassDefOverride(R3BNeulandNeutronReconstructionMon, 0);
};

#endif // R3BNEULANDNEUTRONRECONSTRUCTIONMON_H
