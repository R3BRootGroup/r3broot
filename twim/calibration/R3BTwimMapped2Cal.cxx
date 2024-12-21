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

// ----------------------------------------------------------------
// -----          R3BTwimMapped2Cal source file               -----
// -----    Created 24/11/19 by J.L. Rodriguez-Sanchez        -----
// ----------------------------------------------------------------

// ROOT headers
#include "TClonesArray.h"
#include "TMath.h"
#include "TRandom.h"

// FAIR headers
#include "FairLogger.h"
#include "FairRootManager.h"
#include "FairRunAna.h"
#include "FairRuntimeDb.h"

// TWIM headers
#include "R3BEventHeader.h"
#include "R3BLogger.h"
#include "R3BTwimCalData.h"
#include "R3BTwimCalPar.h"
#include "R3BTwimMapped2Cal.h"
#include "R3BTwimMappedData.h"

constexpr int S444 = 444;
constexpr int S467 = 467;

// R3BTwimMapped2Cal::Default Constructor --------------------------
R3BTwimMapped2Cal::R3BTwimMapped2Cal()
    : R3BTwimMapped2Cal("R3BTwimMapped2Cal", 1)
{
}

// R3BTwimMapped2CalPar::Standard Constructor --------------------------
R3BTwimMapped2Cal::R3BTwimMapped2Cal(const TString& name, Int_t iVerbose)
    : FairTask(name, iVerbose)
    , fNumSec(4)
    , fNumAnodes(16)    // 16 anodes
    , fNumAnodesRef(2)  // 2 anode for TREF
    , fNumAnodesTrig(2) // 2 anode for TRIG
    , fMaxMult(20)
    , fMinDT(-1000)
    , fMaxDT(1000)
    , fNumEParams(2)
    , fNumPosParams(3)
    , fCal_Par(NULL)
    , fTwimMappedDataCA(NULL)
    , fTwimCalDataCA(NULL)
    , fExpId(0)
    , fOnline(kFALSE)
{
    CalEParams.resize(fNumSec, nullptr);
    PosParams.resize(fNumSec, nullptr);
    fTwimCal.resize(fNumSec);
}

void R3BTwimMapped2Cal::SetParContainers()
{
    // Parameter Container
    FairRuntimeDb* rtdb = FairRuntimeDb::instance();
    R3BLOG_IF(fatal, !rtdb, "FairRuntimeDb not found");

    fCal_Par = dynamic_cast<R3BTwimCalPar*>(rtdb->getContainer("twimCalPar"));
    if (!fCal_Par)
    {
        R3BLOG(error, "Couldn't get handle on twimCalPar container");
    }
    else
    {
        R3BLOG(info, "twimCalPar container opened");
    }
    return;
}

void R3BTwimMapped2Cal::SetParameter()
{
    //--- Parameter Container ---
    fNumSec = fCal_Par->GetNumSec();               // Number of sections
    fNumAnodesRef = fCal_Par->GetNumAnodesTRef();  // Anodes for TREF
    fNumAnodesTrig = fCal_Par->GetNumAnodesTrig(); // Anodes for Trig
    fMaxMult = fCal_Par->GetMaxMult();
    fMinDT = fCal_Par->GetMinDT();
    fMaxDT = fCal_Par->GetMaxDT();
    fNumAnodes = fCal_Par->GetNumAnodes();          // Number of anodes per section
    fNumEParams = fCal_Par->GetNumParamsEFit();     // Number of parameters for energy calibration
    fNumPosParams = fCal_Par->GetNumParamsPosFit(); // Number of parameters for position calibration

    R3BLOG(info, "Nb sections: " << fNumSec);
    R3BLOG(info, "Nb anodes: " << fNumAnodes);
    R3BLOG(info, "Nb of twim Ref anodes: " << fNumAnodesRef);
    R3BLOG(info, "Nb of twim Trigger anodes: " << fNumAnodesTrig);
    R3BLOG(info, "Nb of twim Max. multiplicity per anode: " << fMaxMult);
    R3BLOG(info, "Twim Min Drift Time accepted: " << fMinDT);
    R3BLOG(info, "Twim Max Drift Time accepted: " << fMaxDT);
    R3BLOG(info, "Nb parameters from energy cal fit: " << fNumEParams);
    R3BLOG(info, "Nb parameters from position cal fit: " << fNumPosParams);

    CalEParams.resize(fNumSec);
    PosParams.resize(fNumSec);
    for (Int_t s = 0; s < fNumSec; s++)
    {
        // Energy cal parameters
        CalEParams[s] = new TArrayF();
        Int_t array_e = fNumAnodes * fNumEParams;
        CalEParams[s]->Set(array_e);
        CalEParams[s] = fCal_Par->GetAnodeCalEParams(s + 1);
        // Position cal parameters
        PosParams[s] = new TArrayF();
        Int_t array_pos = fNumAnodes * fNumPosParams;
        PosParams[s]->Set(array_pos);
        PosParams[s] = fCal_Par->GetAnodeCalPosParams(s + 1);
    }

    // Count the number of dead anodes or not used
    for (Int_t secId = 0; secId < fNumSec; secId++)
    {
        Int_t numdeadanodes = 0;
        for (Int_t anodeId = 0; anodeId < fNumAnodes; anodeId++)
        {
            if (CalEParams[secId]->GetAt(fNumEParams * anodeId + 1) == -1 ||
                fCal_Par->GetInUse(secId + 1, anodeId + 1) == 0)
            {
                numdeadanodes++;
            }
        }
        R3BLOG(info, "Dead (or not used) anodes in section " << secId + 1 << ": " << numdeadanodes);
    }
    return;
}

// -----   Public method Init   --------------------------------------------
InitStatus R3BTwimMapped2Cal::Init()
{
    R3BLOG(info, "");

    auto* frm = FairRootManager::Instance();
    R3BLOG_IF(fatal, !frm, "FairRootManager not found.");

    header = dynamic_cast<R3BEventHeader*>(frm->GetObject("EventHeader."));
    R3BLOG_IF(error, !header, "EventHeader. not found.");

    if (fExpId == 0) // Obtain global ExpId if it's not set locally.
    {
        fExpId = header->GetExpId();
        R3BLOG(info, "fExpId :" << fExpId);
    }

    // INPUT DATA
    fTwimMappedDataCA = dynamic_cast<TClonesArray*>(frm->GetObject("TwimMappedData"));
    R3BLOG_IF(fatal, !fTwimMappedDataCA, "TwimMappedData not found.");

    // OUTPUT DATA
    fTwimCalDataCA = new TClonesArray("R3BTwimCalData");
    frm->Register("TwimCalData", "TWIM_Cal", fTwimCalDataCA, !fOnline);
    Reset();

    SetParameter();
    return kSUCCESS;
}

// -----   Public method ReInit   ----------------------------------------------
InitStatus R3BTwimMapped2Cal::ReInit()
{
    SetParContainers();
    SetParameter();
    return kSUCCESS;
}

void R3BTwimMapped2Cal::CalAnode::Init()
{
    fmult = 0;
    fE.clear();
    fDT.clear();
    fE.resize(MAXNUMANODE, NAN);
    fDT.resize(MAXNUMANODE, 0);
}

R3BTwimMapped2Cal::CalSection::CalSection()
{
    fAnode.clear();
    fAnode.resize(MAXNUMANODE);

    fTref.clear();
    fTref.resize(MAXNUMTREF);

    fTrig.clear();
    fTrig.resize(MAXNUMTRIG);
}
void R3BTwimMapped2Cal::CalSection::Init()
{
    for (auto& val : fAnode)
    {
        val.Init();
    }
    for (auto& val : fTref)
    {
        val.Init();
    }
    for (auto& val : fTrig)
    {
        val.Init();
    }
}

// -----   Public method Execution   --------------------------------------------
void R3BTwimMapped2Cal::Exec(Option_t* /*option*/) // NOLINT(readability-function-cognitive-complexity)
{
    // Reset entries in output arrays, local arrays
    Reset();

    // Reading the Input -- Mapped Data --
    Int_t nHits = fTwimMappedDataCA->GetEntriesFast();
    if (nHits == 0)
        return;

    for (Int_t sec = 0; sec < fNumSec; sec++)
    {
        fTwimCal.at(sec).Init();
    }
    for (Int_t i = 0; i < nHits; i++)
    {
        auto* mappedData = dynamic_cast<R3BTwimMappedData*>(fTwimMappedDataCA->At(i));
        auto secId = static_cast<Int_t>(mappedData->GetSecID()) - 1;
        auto anodeId = static_cast<Int_t>(mappedData->GetAnodeID()) - 1;
        if (anodeId < fNumAnodes && fCal_Par->GetInUse(secId + 1, anodeId + 1) == 1)
        {
            auto pedestal = CalEParams[secId]->GetAt(fNumEParams * anodeId);
            auto slope = CalEParams[secId]->GetAt(fNumEParams * anodeId + 1);
            fTwimCal.at(secId).GetAnode(anodeId).SetVal(pedestal + slope * mappedData->GetEnergy(),
                                                        mappedData->GetTime());
        }
        else if (anodeId >= fNumAnodes && anodeId < fNumAnodes + fNumAnodesRef)
        {
            auto TrefId = anodeId - fNumAnodes;
            if ((fExpId == S444 || fExpId == S467) && secId != 0)
            {
                secId = 0;
                TrefId++;
            }
            fTwimCal.at(secId).GetTref(TrefId).SetVal(NAN, mappedData->GetTime());
        }
    }

    auto calc_dtime = [&](Int_t sec, Int_t mul_anode, Int_t id_anode, Int_t mul_tref)
    {
        Int_t id_tref = 0;
        if ((fExpId == S444 || fExpId == S467) && id_anode >= fNumAnodes / 2)
        {
            // s444 and s467: 2020
            // anodes 1 to 16 : energy and time
            // anode 17 and 18 : reference time
            // will be cabled anode 19 and 20 : trigger time
            id_tref = 1;
        }
        const auto tdiff = static_cast<Double_t>(fTwimCal.at(sec).GetAnode(id_anode).GetDT(mul_anode) -
                                                 fTwimCal.at(sec).GetTref(id_tref).GetDT(mul_tref));
        //
        auto i_param = fNumPosParams * id_anode;
        auto par0 = PosParams[sec]->GetAt(i_param);
        auto par1 = PosParams[sec]->GetAt(i_param + 1);
        auto par2 = 0.0;
        if (fNumPosParams > 2)
        {
            par2 = PosParams[sec]->GetAt(i_param + 2);
        }
        return par0 + par1 * tdiff + par2 * pow(tdiff, 2);
    };

    // Fill data only if there is TREF signal
    for (Int_t sec = 0; sec < fNumSec; sec++)
    {
        if (fTwimCal.at(sec).GetTref().GetMult() != 1)
        {
            continue;
        }
        for (Int_t anodeId = 0; anodeId < fNumAnodes; anodeId++)
        {
            for (Int_t mult = 0; mult < fTwimCal.at(sec).GetAnode(anodeId).GetMult(); mult++)
            {
                if (fTwimCal.at(sec).GetAnode(anodeId).GetE(mult) <= 0.)
                {
                    continue;
                }
                const Int_t mult_Tref = 0; // mulanode != 1 hits are rejected already.
                auto dtime = calc_dtime(sec, mult, anodeId, mult_Tref);
                if (dtime > fMinDT && dtime < fMaxDT)
                {
                    AddCalData(sec + 1, anodeId + 1, dtime, fTwimCal.at(sec).GetAnode(anodeId).GetE(mult));
                }
            }
        }
    }
    return;
}

// -----   Public method Reset   ------------------------------------------------
void R3BTwimMapped2Cal::Reset()
{
    R3BLOG(debug1, "");
    if (fTwimCalDataCA)
    {
        fTwimCalDataCA->Clear();
    }
}

// -----   Private method AddCalData  --------------------------------------------
R3BTwimCalData* R3BTwimMapped2Cal::AddCalData(Int_t secID, Int_t anodeID, Double_t dtime, Double_t energy)
{
    // It fills the R3BTwimCalData
    TClonesArray& clref = *fTwimCalDataCA;
    Int_t size = clref.GetEntriesFast();
    return new (clref[size]) R3BTwimCalData(secID, anodeID, dtime, energy);
}

ClassImp(R3BTwimMapped2Cal)
