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

// ROOT headers
#include "TClonesArray.h"
#include "TF1.h"
#include "TGraph.h"
#include "TH1D.h"
#include "TMath.h"
#include "TRandom.h"
#include "TVector3.h"

// FAIR headers
#include "FairLogger.h"
#include "FairRootManager.h"
#include "FairRunAna.h"
#include "FairRuntimeDb.h"

#include "R3BLogger.h"

// Musli headers
#include "R3BMusliCalPar.h"
#include "R3BMusliMapped2CalPar.h"
#include "R3BMusliMappedData.h"
#include "R3BMwpcHitData.h"
#include "R3BTGeoPar.h"

// R3BMusliMapped2CalPar: Default Constructor --------------------------
R3BMusliMapped2CalPar::R3BMusliMapped2CalPar()
    : R3BMusliMapped2CalPar("R3BMusliMapped2CalPar", 1, "Mwpc1", "Mwpc2")
{
}

// R3BMusliMapped2CalParPar: Standard Constructor --------------------------
R3BMusliMapped2CalPar::R3BMusliMapped2CalPar(const TString& name,
                                             Int_t iVerbose,
                                             const TString& namedeta,
                                             const TString& namedetb)
    : FairTask(name, iVerbose)
    , fNumSignals(18)
    , fNumGroupsAnodes(15)
    , fMaxMult(20)
    , fMinStatistics(1000)
    , fNumParamsEneFit(2)
    , fNumParamsPosFit(2)
    , fNameDetA(namedeta)
    , fNameDetB(namedetb)
    , fLim_DTrawMin(0)
    , fLim_DTrawMax(24000)
    , EneCalParams(NULL)
    , PosCalParams(NULL)
    , fMusliGeo_Par(NULL)
    , fMwAGeo_Par(NULL)
    , fMwBGeo_Par(NULL)
    , fMusliMappedDataCA(NULL)
    , fMwAHitDataCA(NULL)
    , fMwBHitDataCA(NULL)
{
}

// Virtual R3BMusliMapped2CalPar: Destructor
R3BMusliMapped2CalPar::~R3BMusliMapped2CalPar() { LOG(info) << "R3BMusliMapped2CalPar: Delete instance"; }

void R3BMusliMapped2CalPar::SetParContainers()
{
    // Parameter Container

    FairRuntimeDb* rtdb = FairRuntimeDb::instance();
    R3BLOG_IF(fatal, !rtdb, "R3BMusliMapped2CalPar::SetParContainers(), FairRuntimeDb not found");

    fMusliGeo_Par = dynamic_cast<R3BTGeoPar*>(rtdb->getContainer("MusliGeoPar"));
    R3BLOG_IF(
        error, !fMusliGeo_Par, "R3BMusliMapped2CalPar::SetParContainers() Couldn´t access to MusliGeoPar container.");

    fMwAGeo_Par = dynamic_cast<R3BTGeoPar*>(rtdb->getContainer(fNameDetA + "GeoPar"));
    R3BLOG_IF(error,
              !fMwAGeo_Par,
              "R3BMusliMapped2CalPar::SetParContainers() Couldn´t access to " + fNameDetA + "GeoPar container.");

    fMwBGeo_Par = dynamic_cast<R3BTGeoPar*>(rtdb->getContainer(fNameDetB + "GeoPar"));
    R3BLOG_IF(error,
              !fMwBGeo_Par,
              "R3BMusliMapped2CalPar::SetParContainers() Couldn´t access to " + fNameDetB + "GeoPar container.");
}

// -----   Public method Init   --------------------------------------------
InitStatus R3BMusliMapped2CalPar::Init()
{
    R3BLOG(info, "R3BMusliMapped2CalPar: Init");

    // INPUT DATA
    FairRootManager* rootManager = FairRootManager::Instance();
    if (!rootManager)
    {
        return kFATAL;
    }

    fMusliMappedDataCA = dynamic_cast<TClonesArray*>(rootManager->GetObject("MusliMappedData"));
    if (!fMusliMappedDataCA)
    {
        LOG(error) << "R3BMusliMapped2CalPar: MusliMappedData not found";
        return kFATAL;
    }

    // get access to hit data of the first Mwpc (A)
    fMwAHitDataCA = dynamic_cast<TClonesArray*>(rootManager->GetObject(fNameDetA + "HitData"));
    if (!fMwAHitDataCA)
    {
        LOG(error) << "R3BMusliMapped2CalPar: " + fNameDetA + "HitData not found";
        return kFATAL;
    }

    // get access to hit data of the second Mwpc (B)
    fMwBHitDataCA = dynamic_cast<TClonesArray*>(rootManager->GetObject(fNameDetB + "HitData"));
    if (!fMwBHitDataCA)
    {
        LOG(error) << "R3BMusliMapped2CalPar: " + fNameDetB + "HitData not found";
        return kFATAL;
    }

    // for the futur calibration parameters
    FairRuntimeDb* rtdb = FairRuntimeDb::instance();
    if (!rtdb)
    {
        return kFATAL;
    }

    fCal_Par = dynamic_cast<R3BMusliCalPar*>(rtdb->getContainer("musliCalPar"));
    if (!fCal_Par)
    {
        LOG(error) << "R3BMusliMapped2CalPar::Couldn't get handle on musliCalPar container";
        return kFATAL;
    }

    // define the 2D histograms for the position calibration of the MusLI
    char Name[255];
    fh2_XvsDT = new TH2D*[fNumGroupsAnodes];
    for (Int_t i = 0; i < fNumGroupsAnodes; i++)
    {
        sprintf(Name, "Musli_XvsDT_GroupAnode%02d", i);
        fh2_XvsDT[i] = new TH2D(Name, Name, 1200, 0, 24000, 1300, -20, 110);
    }

    return kSUCCESS;
}

// -----   Public method ReInit   ----------------------------------------------
InitStatus R3BMusliMapped2CalPar::ReInit()
{

    SetParContainers();
    return kSUCCESS;
}

// -----   Public method Execution   --------------------------------------------
void R3BMusliMapped2CalPar::Exec(Option_t* option)
{

    Int_t nHits = fMusliMappedDataCA->GetEntriesFast();
    if (nHits == 0)
        return;
    R3BMusliMappedData** mappedData = new R3BMusliMappedData*[nHits];

    for (Int_t i = 0; i < fNumSignals; i++)
    {
        multMap[i] = 0;
        for (Int_t j = 0; j < fMaxMult; j++)
        {
            fEneMap[j][i] = 0;
            fTimeMap[j][i] = 0;
        }
    }

    for (Int_t i = 0; i < nHits; i++)
    {
        mappedData[i] = dynamic_cast<R3BMusliMappedData*>(fMusliMappedDataCA->At(i));
        UInt_t signal_id = mappedData[i]->GetSignal() - 1; // signal_id is 0-based [0..17]
        if (multMap[signal_id] < fMaxMult)
        {
            fEneMap[multMap[signal_id]][signal_id] = mappedData[i]->GetEnergy();
            fTimeMap[multMap[signal_id]][signal_id] = mappedData[i]->GetTime();
            multMap[signal_id]++;
        }
    }

    // Continue only if there is only one TREF signal, otherwise return (->exit Exec())
    if (multMap[16] != 1)
        return;

    // Reading the Input -- MwpcA Hit Data -- only if one (X,Y) hit only
    nHits = fMwAHitDataCA->GetEntriesFast();
    if (nHits != 1)
        return;
    else
    {
        R3BMwpcHitData* hitDataMwA = dynamic_cast<R3BMwpcHitData*>(fMwAHitDataCA->At(0));
        fXA = hitDataMwA->GetX() + fMwAGeo_Par->GetPosX();
    }

    // Reading the Input -- MwpcB Hit Data --
    nHits = fMwBHitDataCA->GetEntriesFast();
    if (nHits != 1)
        return;
    else
    {
        R3BMwpcHitData* hitDataMwB = dynamic_cast<R3BMwpcHitData*>(fMwBHitDataCA->At(0));
        fXB = hitDataMwB->GetX() + fMwBGeo_Par->GetPosX();
    }

    // --- Calibration from Time at Mapped level to dt_cal in [mm] at CAL level
    // REMINDER :
    // for fSignal = [1..8]  8 groups of 2 anodes (width 2A = 50 mm)
    // for fSignal = [9..12] 4 groups of 4 anodes (width of 4A = 100 mm)
    // for fSignal = [13,14] 2 groups of 8 anodes (width of 8A = 200 mm)
    // for fSignal = [15]    1 group of 16 anodes (width of 16A = 400 mm = active volume)
    Double_t slope = (fXB - fXA) / (fMwBGeo_Par->GetPosZ() - fMwAGeo_Par->GetPosZ());
    Double_t musliZpos[fNumGroupsAnodes]; // position in Z of each group of anodes
    Double_t musliXpos[fNumGroupsAnodes]; // position in X extrapolated from MwA & MwB
    Double_t DTraw[fNumGroupsAnodes];     // raw drift time measured for each group of anodes

    for (Int_t i = 0; i < 8; i++)
        musliZpos[i] = fMusliGeo_Par->GetPosZ() * 10. - 175. + i * 50.;
    for (Int_t i = 8; i < 12; i++)
        musliZpos[i] = fMusliGeo_Par->GetPosZ() * 10. - 150. + (i - 8) * 100.;
    for (Int_t i = 12; i < 14; i++)
        musliZpos[i] = fMusliGeo_Par->GetPosZ() * 10. - 100. + (i - 12) * 200.;
    musliZpos[14] = fMusliGeo_Par->GetPosZ() * 10.;

    // Calculating the raw drift times and fill histograms and graphs
    for (Int_t i = 0; i < fNumGroupsAnodes; i++)
    {
        DTraw[i] = -10000.;
        if (multMap[i] == 1)
        {
            DTraw[i] = fTimeMap[0][i] - fTimeMap[0][16];
            musliXpos[i] = fXA + slope * (musliZpos[i] - fMwAGeo_Par->GetPosZ());
            fh2_XvsDT[i]->Fill(DTraw[i], musliXpos[i]);
        }
    }

    if (mappedData)
        delete[] mappedData;

    return;
}

void R3BMusliMapped2CalPar::CalculatePosCalPar()
{
    R3BLOG(info, "R3BMusliMapped2CalPar::CalculatePosCalPar()");
    TH1D* h1_pfx_XvsDT[fNumGroupsAnodes];
    TF1* fit_pol1_XvsDT[fNumGroupsAnodes];
    char Name[255];

    for (Int_t i = 0; i < fNumGroupsAnodes; i++)
    {
        if (fh2_XvsDT[i]->Integral() > fMinStatistics)
        {
            sprintf(Name, "pol1_XvsDT_GroupAnode%02d", i);
            fit_pol1_XvsDT[i] = new TF1(Name, "[0]+[1]*x", fLim_DTrawMin, fLim_DTrawMax);

            sprintf(Name, "pfx_XvsDT_GroupAnode%02d", i);
            fh2_XvsDT[i]->ProfileX(Name);
            h1_pfx_XvsDT[i] = dynamic_cast<TH1D*>(gDirectory->FindObject(Name));
            h1_pfx_XvsDT[i]->Fit(fit_pol1_XvsDT[i], "R+");

            Double_t offset_alignment = fit_pol1_XvsDT[i]->GetParameter(0);
            fCal_Par->SetPosParams(offset_alignment, i + 1, 0); // i+1 = signal_id 1-based

            Double_t drift_velocity = fit_pol1_XvsDT[i]->GetParameter(1);
            fCal_Par->SetPosParams(drift_velocity, i + 1, 1); // i+1 = signal_id 1-based
        }
    }
}

void R3BMusliMapped2CalPar::FinishTask()
{

    fCal_Par->SetNumSignals(fNumSignals);
    fCal_Par->SetNumSignals(fNumSignals);
    fCal_Par->SetMaxMult(fMaxMult);
    fCal_Par->SetNumParamsEneFit(fNumParamsEneFit);
    fCal_Par->SetNumParamsPosFit(fNumParamsPosFit);
    CalculatePosCalPar();

    fCal_Par->printParams();

    fCal_Par->setChanged();
}

ClassImp(R3BMusliMapped2CalPar);
