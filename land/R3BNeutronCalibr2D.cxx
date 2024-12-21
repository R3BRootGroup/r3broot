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

// -----------------------------------------------------------------------------
// -----                       R3BNeutronCalibr2D                          -----
// -----                 Created 27-02-2012 by D.Kresan                    -----
// -----------------------------------------------------------------------------

#include <iostream>

#include "TClonesArray.h"
#include "TH1F.h"
#include "TH2F.h"

#include "FairLogger.h"
#include "FairRootManager.h"

#include "R3BLandDigi.h"
#include "R3BNeutronCalibr2D.h"

using std::cout;
using std::endl;

ClassImp(R3BNeutronCalibr2D)

    // -----------------------------------------------------------------------------
    R3BNeutronCalibr2D::R3BNeutronCalibr2D()
    : FairTask("Calibration of NEULAND cuts")
    , fEventNo(0)
{
}
// -----------------------------------------------------------------------------

// -----------------------------------------------------------------------------
R3BNeutronCalibr2D::~R3BNeutronCalibr2D() {}
// -----------------------------------------------------------------------------

// -----------------------------------------------------------------------------
InitStatus R3BNeutronCalibr2D::Init()
{
    // Get input array
    FairRootManager* ioman = FairRootManager::Instance();
    if (!ioman)
        LOG(fatal) << "Init: No FairRootManager";
    fArrayDigi = dynamic_cast<TClonesArray*>(ioman->GetObject("LandDigi"));
    fArrayCluster = dynamic_cast<TClonesArray*>(ioman->GetObject("NeuLandCluster"));

    CreateHistograms();

    return kSUCCESS;
}
// -----------------------------------------------------------------------------

// -----------------------------------------------------------------------------
void R3BNeutronCalibr2D::Exec(Option_t* opt)
{
    //  if(0 == (fEventNo%10)) {
    //    cout << "-I- R3BNeutronCalibr2D : Event #: " << fEventNo << "." << endl;
    //  }
    fEventNo += 1;

    // Loop over digis
    Float_t eTot = 0;
    R3BLandDigi* digi;
    for (Int_t id = 0; id < fArrayDigi->GetEntries(); id++)
    {
        digi = dynamic_cast<R3BLandDigi*>(fArrayDigi->At(id));
        eTot += digi->GetQdc();
    }

    fh_etot->Fill(eTot);
    fh_ncl_etot->Fill(eTot, fArrayCluster->GetEntries());
}
// -----------------------------------------------------------------------------

// -----------------------------------------------------------------------------
void R3BNeutronCalibr2D::CreateHistograms()
{
    if (fBeamEnergy < 300)
    {
        fh_etot = new TH1F("h_etot", "Total light", 1500, 0., 1500.);
        fh_ncl_etot = new TH2F("h_ncl_etot", "Number of clusters vs. total light", 1500, 0., 1500., 150, -0.5, 149.5);
    }
    else
    {
        fh_etot = new TH1F("h_etot", "Total light", 150, 0., 1500.);
        fh_ncl_etot = new TH2F("h_ncl_etot", "Number of clusters vs. total light", 150, 0., 1500., 150, -0.5, 149.5);
    }
}
// -----------------------------------------------------------------------------

// -----------------------------------------------------------------------------
void R3BNeutronCalibr2D::Finish()
{
    fh_etot->Write();
    fh_ncl_etot->Write();
}
// -----------------------------------------------------------------------------
