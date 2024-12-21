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

// Fair headers
#include "FairLogger.h"
#include "FairRootManager.h"
#include "FairRunAna.h"
#include "FairRuntimeDb.h"

// FrsSci headers
#include "R3BFrsSciMapped2Tcal.h"
#include "R3BFrsSciMappedData.h"
#include "R3BFrsSciTcalData.h"

// --- Default Constructor
R3BFrsSciMapped2Tcal::R3BFrsSciMapped2Tcal()
    : FairTask("R3BFrsSciMapped2Tcal", 1)
    , fNevent(0)
    , fMapped(NULL)
    , fTcalPar(NULL)
    , fTcal(NULL)
    , fOnline(kFALSE)
{
}

// --- Standard Constructor
R3BFrsSciMapped2Tcal::R3BFrsSciMapped2Tcal(const char* name, Int_t iVerbose)
    : FairTask(name, iVerbose)
    , fNevent(0)
    , fMapped(NULL)
    , fTcalPar(NULL)
    , fTcal(NULL)
    , fOnline(kFALSE)
{
}

// --- Destructor
R3BFrsSciMapped2Tcal::~R3BFrsSciMapped2Tcal()
{
    LOG(info) << "R3BFrsSciMapped2Tcal: Delete instance";
    if (fMapped)
    {
        delete fMapped;
    }
    if (fTcal)
    {
        delete fTcal;
    }
}

// --- Parameter container : reading FrsSciTcalPar from FairRuntimeDb
void R3BFrsSciMapped2Tcal::SetParContainers()
{
    FairRuntimeDb* rtdb = FairRuntimeDb::instance();
    if (!rtdb)
    {
        LOG(error) << "FairRuntimeDb not opened!";
    }

    fTcalPar = (R3BFrsSciTcalPar*)rtdb->getContainer("FrsSciTcalPar");
    if (!fTcalPar)
    {
        LOG(error) << "R3BFrsSciMapped2Tcal::SetParContainers() : Could not get access to FrsSciTcalPar-Container.";
        return;
    }
    else
    {
        LOG(info) << "R3BFrsSciMapped2Tcal::SetParContainers() : FrsSciTcalPar-Container found";
    }
}

InitStatus R3BFrsSciMapped2Tcal::Init()
{
    LOG(info) << "R3BFrsSciMapped2Tcal::Init()";

    FairRootManager* rm = FairRootManager::Instance();
    if (!rm)
    {
        LOG(error) << "R3BFrsSciMapped2Tcal::Init() Couldn't instance the FairRootManager";
        return kFATAL;
    }

    // Input data
    fMapped = (TClonesArray*)rm->GetObject("FrsSciMappedData");
    if (!fMapped)
    {
        LOG(error) << "R3BFrsSciMapped2Tcal::Init() Couldn't get handle on FrsSciMappedData container";
        return kFATAL;
    }
    else
        LOG(info) << "R3BFrsSciMapped2Tcal::Init() FrsSciMappedData items found";

    // Register output array in tree
    fTcal = new TClonesArray("R3BFrsSciTcalData", 25);
    if (!fOnline)
    {
        rm->Register("FrsSciTcalData", "FrsSci", fTcal, kTRUE);
    }
    else
    {
        rm->Register("FrsSciTcalData", "FrsSci", fTcal, kFALSE);
    }

    return kSUCCESS;
}

InitStatus R3BFrsSciMapped2Tcal::ReInit()
{
    SetParContainers();
    return kSUCCESS;
}

void R3BFrsSciMapped2Tcal::Exec(Option_t* option)
{
    UShort_t det;
    UShort_t pmt;
    UInt_t tf;
    UInt_t tc;
    Double_t tns = 0.0;

    // Reset entries in output arrays, local arrays
    Reset();

    // Loop over the entries of the Mapped TClonesArray
    UInt_t nHitsPerEvent_FrsSci = fMapped->GetEntries();
    for (UInt_t ihit = 0; ihit < nHitsPerEvent_FrsSci; ihit++)
    {
        R3BFrsSciMappedData* hit = (R3BFrsSciMappedData*)fMapped->At(ihit);
        if (!hit)
            continue;
        det = hit->GetDetector();
        pmt = hit->GetPmt();
        tf = hit->GetTimeFine();
        tc = hit->GetTimeCoarse();

        if ((det < 1) || (det > (UShort_t)fTcalPar->GetNumDets()))
        {
            LOG(info) << "R3BFrsSciMapped2Tcal::Exec() : In FrsSciMappedData, process fDetetector = " << det
                      << "is out of range [1:" << fTcalPar->GetNumDets() << ", item skipped ";
            continue;
        }
        if ((pmt < 1) || (pmt > (UShort_t)fTcalPar->GetNumPmts()))
        {
            LOG(info) << "R3BFrsSciMapped2Tcal::Exec() : In FrsSciMappedData, process fPmt = " << pmt
                      << "is out of range [1:" << fTcalPar->GetNumPmts() << "], item skipped ";
            continue;
        }
        tns = CalculateTimeNs(det, pmt, tf, tc);
        AddTcalData(det, pmt, tns);
    }

    if (nHitsPerEvent_FrsSci != fTcal->GetEntries())
        LOG(warn) << "R3BFrsSciMapped2Tcal::Exec() mismatch between TClonesArray entries ";

    ++fNevent;
    return;
}

// -----   Public method Reset   ------------------------------------------------
void R3BFrsSciMapped2Tcal::Reset()
{
    LOG(debug) << "Clearing TcalData Structure";
    if (fTcal)
        fTcal->Clear();
}

// -----   Private method to calculate the time in ns  --------------------------
Double_t R3BFrsSciMapped2Tcal::CalculateTimeNs(UShort_t iDet, UShort_t iPmt, UInt_t iTf, UInt_t iTc)
{
    UInt_t rank = iTf + (UInt_t)fTcalPar->GetNumPars() * ((iDet - 1) * (UInt_t)fTcalPar->GetNumPmts() + (iPmt - 1));
    Double_t iPar = fTcalPar->GetOneTcalParam(rank);
    Double_t r = (Double_t)rand.Rndm() - 0.5;
    Double_t iTf_ns;
    Double_t iTc_ns = (Double_t)iTc * 5.;

    if (r < 0)
    {
        Double_t iParPrev = fTcalPar->GetOneTcalParam(rank - 1);
        iTf_ns = iPar + r * (iPar - iParPrev);
    }
    else
    {
        Double_t iParNext = fTcalPar->GetOneTcalParam(rank + 1);
        iTf_ns = iPar + r * (iParNext - iPar);
    }

    return (iTc_ns - iTf_ns);
}

// -----   Private method AddCalData  --------------------------------------------
R3BFrsSciTcalData* R3BFrsSciMapped2Tcal::AddTcalData(UShort_t det, UShort_t pmt, Double_t tns)
{
    // It fills the R3BFrsSciTcalData
    TClonesArray& clref = *fTcal;
    Int_t size = clref.GetEntriesFast();
    return new (clref[size]) R3BFrsSciTcalData(det, pmt, tns);
}

ClassImp(R3BFrsSciMapped2Tcal)
