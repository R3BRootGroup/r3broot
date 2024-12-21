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

#include "R3BFrsSciMapped2TcalPar.h"

#include "FairLogger.h"
#include "FairRootManager.h"
#include "FairRuntimeDb.h"
#include "R3BFrsSciMappedData.h"
#include "R3BFrsSciTcalPar.h"
#include "TClonesArray.h"
#include "TF1.h"
#include "TH1F.h"
#include "TMath.h"
#include "TObjArray.h"
#include "TRandom.h"

#include <iostream>
#include <stdlib.h>

// R3BFrsSciMapped2TcalPar: Default Constructor --------------------------
R3BFrsSciMapped2TcalPar::R3BFrsSciMapped2TcalPar()
    : FairTask("R3BFrsSciMapped2TcalPar", 1)
    , fMapped(NULL)
    , fNumDets(3)
    , fNumPmts(3)
    , fNumPars(1000)
    , fMinStat(100000)
    , fTcalPar(NULL)
    , fOutputFile(NULL)
{
}

// R3BFrsSciMapped2TcalPar: Standard Constructor --------------------------
R3BFrsSciMapped2TcalPar::R3BFrsSciMapped2TcalPar(const char* name, Int_t iVerbose)
    : FairTask(name, iVerbose)
    , fMapped(NULL)
    , fNumDets(3)
    , fNumPmts(3)
    , fNumPars(1000)
    , fMinStat(100000)
    , fTcalPar(NULL)
    , fOutputFile(NULL)

{
}

// R3BFrsSciMapped2TcalPar: Destructor ----------------------------------------
R3BFrsSciMapped2TcalPar::~R3BFrsSciMapped2TcalPar()
{
    if (fMapped)
        delete fMapped;
    if (fTcalPar)
        delete fTcalPar;
}

// -----   Public method Init   --------------------------------------------
InitStatus R3BFrsSciMapped2TcalPar::Init()
{
    LOG(info) << "R3BFrsSciMapped2TcalPar::Init()";

    FairRootManager* rm = FairRootManager::Instance();
    if (!rm)
    {
        LOG(error) << "R3BFrsSciMapped2TcalPar::Init() FairRootManager not found";
        return kFATAL;
    }

    // input data
    fMapped = (TClonesArray*)rm->GetObject("FrsSciMappedData"); // see Instance->Register in R3BFrsSciReader.cxx
    if (!fMapped)
    {
        LOG(error) << "R3BFrsSciMapped2TcalPar::Init() FrsSciMappedData not found";
        return kFATAL;
    }

    // FrsSci Tcal parameters container
    FairRuntimeDb* rtdb = FairRuntimeDb::instance();
    if (!rtdb)
    {
        LOG(error) << "R3BFrsSciMapped2TcalPar::Init() FairRuntimeDb not found";
        return kFATAL;
    }
    fTcalPar = (R3BFrsSciTcalPar*)rtdb->getContainer("FrsSciTcalPar");
    if (!fTcalPar)
    {
        LOG(error) << "R3BFrsSciMapped2TcalPar::Init() Couldn't get handle on FrsSciTcalPar container";
        return kFATAL;
    }
    else
    {
        fTcalPar->SetNumDets((Int_t)fNumDets);
        fTcalPar->SetNumPmts((Int_t)fNumPmts);
        fTcalPar->SetNumPars((Int_t)fNumPars);
    }
    // histograms
    char name[100];
    fh_TimeFineBin = new TH1F*[fNumDets * fNumPmts];
    fh_TimeFineNs = new TH1D*[fNumDets * fNumPmts];
    for (UShort_t det = 0; det < fNumDets; det++)
    {
        for (UShort_t pmt = 0; pmt < fNumPmts; pmt++)
        {
            sprintf(name, "TimeFineBin_Sci%i_Pmt%i", det + 1, pmt + 1);
            fh_TimeFineBin[det * fNumPmts + pmt] = new TH1F(name, name, fNumPars, 0, fNumPars);
            sprintf(name, "TimeFineNs_Sci%i_Pmt%i", det + 1, pmt + 1);
            fh_TimeFineNs[det * fNumPmts + pmt] = new TH1D(name, name, fNumPars, 0, fNumPars);
        }
    }

    return kSUCCESS;
}

// -----   Public method ReInit   --------------------------------------------
InitStatus R3BFrsSciMapped2TcalPar::ReInit() { return kSUCCESS; }

// -----   Public method Exec   --------------------------------------------
void R3BFrsSciMapped2TcalPar::Exec(Option_t* opt)
{

    UInt_t nHitsSci = fMapped->GetEntries(); // can be very high especially for S2 detector
    UShort_t det;
    UShort_t pmt;
    for (UInt_t ihit = 0; ihit < nHitsSci; ihit++)
    {
        R3BFrsSciMappedData* hitSci = (R3BFrsSciMappedData*)fMapped->At(ihit);
        if (!hitSci)
        {
            LOG(warn) << "R3BFrsSciMapped2TcalPar::Exec() : could not get hitSci";
            continue; // should not happen
        }
        det = hitSci->GetDetector() - 1;
        pmt = hitSci->GetPmt() - 1;
        fh_TimeFineBin[det * fNumPmts + pmt]->Fill(hitSci->GetTimeFine());

    } // end of loop over the number of hits per event in MappedSci
}

// ---- Public method Reset   --------------------------------------------------
void R3BFrsSciMapped2TcalPar::Reset() {}

void R3BFrsSciMapped2TcalPar::FinishEvent() {}

// ---- Public method Finish   --------------------------------------------------
void R3BFrsSciMapped2TcalPar::FinishTask()
{
    CalculateVftxTcalParams();
    fTcalPar->printParams();
}

//------------------
void R3BFrsSciMapped2TcalPar::CalculateVftxTcalParams()
{
    LOG(info) << "R3BFrsSciMapped2TcalPar::CalculateVftxTcalParams()";

    fTcalPar->SetNumDets((Int_t)fNumDets);
    fTcalPar->SetNumPmts((Int_t)fNumPmts);
    fTcalPar->SetNumPars((Int_t)fNumPars);

    Double_t IntegralTot;
    Double_t IntegralPartial;
    Double_t Bin2Ns[fNumPars];

    LOG(info) << " min stat is set to " << fMinStat;

    for (Int_t sig = 0; sig < fNumDets * fNumPmts; sig++)
    {
        if (fh_TimeFineBin[sig]->Integral() > fMinStat)
        {
            LOG(info) << " Process calibration for signal " << sig;
            IntegralTot = (Double_t)fh_TimeFineBin[sig]->Integral();
            IntegralPartial = 0.;
            for (Int_t bin = 0; bin < fNumPars; bin++)
            {
                IntegralPartial += (Double_t)fh_TimeFineBin[sig]->GetBinContent(bin + 1);
                Bin2Ns[bin] = 5. * IntegralPartial / IntegralTot;
                fh_TimeFineNs[sig]->SetBinContent(bin + 1, Bin2Ns[bin]);
                fTcalPar->SetOneTcalParam(Bin2Ns[bin], sig * fNumPars + bin);
            }
        }
        else
        {
            for (Int_t bin = 0; bin < fNumPars; bin++)
            {
                fTcalPar->SetOneTcalParam(0, sig * fNumPars + bin);
            }
        }
        fh_TimeFineNs[sig]->Write(); // empty histo if stat <fMinStatistics
        fh_TimeFineBin[sig]->Write();
    }
    fTcalPar->setChanged();
    return;
}

ClassImp(R3BFrsSciMapped2TcalPar);
