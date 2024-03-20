/******************************************************************************
 *   Copyright (C) 2019 GSI Helmholtzzentrum für Schwerionenforschung GmbH    *
 *   Copyright (C) 2019-2024 Members of R3B Collaboration                     *
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
#include "R3BFrsSciTcal2Cal.h"

#define SPEED_OF_LIGHT_MNS 0.299792458

// --- Default Constructor
R3BFrsSciTcal2Cal::R3BFrsSciTcal2Cal()
    : FairTask("R3BFrsSciTcal2Cal", 1)
    , fNevent(0)
    , fTcal(NULL)
    , fCalPar(NULL)
    , fPosCal(NULL)
    , fTofCal(NULL)
    , fOnline(kFALSE)
{
}

// --- Standard Constructor
R3BFrsSciTcal2Cal::R3BFrsSciTcal2Cal(const char* name, Int_t iVerbose)
    : FairTask(name, iVerbose)
    , fNevent(0)
    , fTcal(NULL)
    , fCalPar(NULL)
    , fPosCal(NULL)
    , fTofCal(NULL)
    , fOnline(kFALSE)
{
}

// --- Destructor
R3BFrsSciTcal2Cal::~R3BFrsSciTcal2Cal()
{
    LOG(info) << "R3BFrsSciTcal2Cal: Delete instance";
    if (fTcal)
    {
        delete fTcal;
    }
    if (fPosCal)
    {
        delete fPosCal;
    }
    if (fTofCal)
    {
        delete fTofCal;
    }
}

// --- Parameter container : reading FrsSciTcalPar from FairRuntimeDb
void R3BFrsSciTcal2Cal::SetParContainers()
{
    FairRuntimeDb* rtdb = FairRuntimeDb::instance();
    if (!rtdb)
    {
        LOG(error) << "FairRuntimeDb not opened!";
    }

    fCalPar = (R3BFrsSciCalPar*)rtdb->getContainer("FrsSciCalPar");
    if (!fCalPar)
    {
        LOG(error) << "R3BFrsSciTcal2Cal::SetParContainers() : Could not get access to FrsSciCalPar-Container.";
        return;
    }
    else
    {
        LOG(info) << "R3BFrsSciTcal2Cal::SetParContainers() : FrsSciCalPar-Container found ";
    }
}

InitStatus R3BFrsSciTcal2Cal::Init()
{
    LOG(info) << "R3BFrsSciTcal2Cal::Init()";

    FairRootManager* rm = FairRootManager::Instance();
    if (!rm)
    {
        LOG(error) << "R3BFrsSciTcal2Cal::Init() Couldn't instance the FairRootManager";
        return kFATAL;
    }

    // Input data
    fTcal = (TClonesArray*)rm->GetObject("FrsSciTcalData");
    if (!fTcal)
    {
        LOG(error) << "R3BFrsSciTcal2Cal::Init() Couldn't get handle on FrsSciTcalData container";
        return kFATAL;
    }
    else
        LOG(info) << "R3BFrsSciTcal2Cal::Init() FrsSciTcalData items found";

    // Register output arrays in tree
    fPosCal = new TClonesArray("R3BFrsSciPosCalData", 25);
    if (!fOnline)
    {
        rm->Register("FrsSciPosCalData", "FrsSci", fPosCal, kTRUE);
    }
    else
    {
        rm->Register("FrsSciPosCalData", "FrsSci", fPosCal, kFALSE);
    }

    fTofCal = new TClonesArray("R3BFrsSciTofCalData", 25);
    if (!fOnline)
    {
        rm->Register("FrsSciTofCalData", "FrsSci", fTofCal, kTRUE);
    }
    else
    {
        rm->Register("FrsSciTofCalData", "FrsSci", fTofCal, kFALSE);
    }

    return kSUCCESS;
}

InitStatus R3BFrsSciTcal2Cal::ReInit()
{
    SetParContainers();
    return kSUCCESS;
}

void R3BFrsSciTcal2Cal::Exec(Option_t* option)
{
    // Reset entries in output arrays, local arrays
    Reset();
    // Variables to read the Tcal data
    UShort_t nDets = (UShort_t)fCalPar->GetNumDets();
    UShort_t nPmts = (UShort_t)fCalPar->GetNumPmts();
    // UShort_t nDets = 2;
    // UShort_t nPmts = 3;
    UShort_t iDet = 0;
    UShort_t iPmt = 0;
    Double_t iTraw[nDets * nPmts][64]; // 64 hits max per VFTX channel
    UShort_t mult[nDets * nPmts];
    // ULong64_t maskR[nDets];
    // ULong64_t maskL[nDets];
    for (UShort_t i = 0; i < nDets * nPmts; i++)
    {
        mult[i] = 0;
        for (UShort_t j = 0; j < 64; j++)
        {
            iTraw[i][j] = 0;
        }
    }

    // Loop over the entries of the Tcal TClonesArray
    UInt_t nHitsPerEvent_FrsSci = fTcal->GetEntries();
    for (UInt_t ihit = 0; ihit < nHitsPerEvent_FrsSci; ihit++)
    {
        R3BFrsSciTcalData* hit = (R3BFrsSciTcalData*)fTcal->At(ihit);
        if (!hit)
            continue;
        iDet = hit->GetDetector() - 1;
        iPmt = hit->GetPmt() - 1;
        if (mult[iDet * nPmts + iPmt] < 64) // always true, max mult=64 per vftx channel
        {
            iTraw[iDet * nPmts + iPmt][mult[iDet * nPmts + iPmt]] = hit->GetRawTimeNs();
            mult[iDet * nPmts + iPmt]++;
        }
    }

    // Select the proper hit in the multihit tcal data
    // if nDets==1 -> Use the Position : not really selective, better use LOS (need to be coded)
    // if nDets>1  -> Use the Tof
    if (nDets == 1 && mult[0] > 0 && mult[1] > 0)
    {
        // Variables to fill PosCal and TofCal data
        Double_t iRawTime = NAN;
        Double_t iRawTime_wtref = NAN;
        Float_t iRawPos = NAN;
        Float_t iCalPos = NAN;
        // maskL[0] = 0x0;
        // maskR[0] = 0x0;
        for (UShort_t multR = 0; multR < mult[0]; multR++)
        {
            for (UShort_t multL = 0; multL < mult[1]; multL++)
            {
                // RawPos = TrawRIGHT - TrawLEFT corresponds to x increasing from RIGHT to LEFT
                iRawPos = iTraw[0][multR] - iTraw[1][multL];
                // if the raw position is outside the range: continue
                if (iRawPos < fCalPar->GetMinPosAtRank(0))
                    continue;
                if (iRawPos > fCalPar->GetMaxPosAtRank(0))
                    continue;
                // if the left or right hit has already been used, continue
                // if ((((maskR[0] >> multR) & (0x1)) == 1) || (((maskL[0] > multL) & (0x1)) == 1))
                //    continue;
                iRawTime = 0.5 * (iTraw[0][multR] + iTraw[1][multL]);

                // Since the searching window is set to wait for the time reference (LOS) to be arrived, the last hit
                // will be the correct one for the most cases. VFTX seems to be LIFO module?
                iRawTime_wtref = iRawTime - iTraw[2][0];
                // tag which hit is used
                // maskR[0] |= (0x1) << multR;
                // maskL[0] |= (0x1) << multL;

                iCalPos = fCalPar->GetPosCalGainAtRank(0) * iRawPos + fCalPar->GetPosCalOffsetAtRank(0);
                AddPosCalData(1, iRawTime, iRawTime_wtref, iRawPos, iCalPos);
            } // end of loop over the hits of the left PMTs
        }     // end of loop over the hits of the right PMTs
    }         // end of if (nDets==1)
    else if (nDets > 1)
    {
        // Variables to fill PosCal and TofCal data
        Double_t iRawTimeSta = NAN, iRawTimeSto = NAN;
        Float_t iRawPosSta = NAN;
        Float_t iRawPosSto = NAN;
        Double_t iRawTof = NAN, iCalTof = -1, iCalVelo = NAN;
        Double_t iBeta = NAN, iGamma = NAN, iBRho = NAN, iAoQ = NAN;
        Int_t selectLeftHit[nDets];
        Int_t selectRightHit[nDets];
        for (UShort_t det = 0; det < nDets; det++)
        {
            selectLeftHit[det] = -1;
            selectRightHit[det] = -1;
        }

        // stop detector is the last detector
        Int_t dSto = nDets - 1; // 0-based
        UShort_t selectRightHit_dSto[nDets - 1];
        UShort_t selectLeftHit_dSto[nDets - 1];
        for (UShort_t det = 0; det < nDets - 1; det++)
        {
            selectLeftHit_dSto[det] = -1;
            selectRightHit_dSto[det] = -1;
        }
        Bool_t selectHit = kFALSE;
        Bool_t selectSameHit = kTRUE;

        // if data in Sto
        if (mult[dSto * nPmts] > 0 && mult[dSto * nPmts + 1] > 0)
        {
            for (UShort_t multRsto = 0; multRsto < mult[dSto * nPmts]; multRsto++)
            {
                for (UShort_t multLsto = 0; multLsto < mult[dSto * nPmts + 1]; multLsto++)
                {
                    iRawPosSto = (Float_t)(iTraw[dSto * nPmts][multRsto] - iTraw[dSto * nPmts + 1][multLsto]);
                    // process only data with valid position
                    if ((iRawPosSto < fCalPar->GetMinPosAtRank(dSto)) || (iRawPosSto > fCalPar->GetMaxPosAtRank(dSto)))
                        continue;
                    // loop over all Sta detector
                    for (UShort_t dSta = 0; dSta < nDets - 1; dSta++)
                    {
                        // if data in Sta
                        if (mult[dSta * nPmts] > 0 && mult[dSta * nPmts + 1] > 0)
                        {
                            for (UShort_t multRsta = 0; multRsta < mult[dSta * nPmts]; multRsta++)
                            {
                                for (UShort_t multLsta = 0; multLsta < mult[dSta * nPmts + 1]; multLsta++)
                                {
                                    // process only data with vali position
                                    iRawPosSta =
                                        (Float_t)(iTraw[dSta * nPmts][multRsta] - iTraw[dSta * nPmts + 1][multLsta]);
                                    if ((iRawPosSta < fCalPar->GetMinPosAtRank(dSta)) ||
                                        (iRawPosSta > fCalPar->GetMaxPosAtRank(dSta)))
                                        continue;
                                    // if position in Sta and Sto are valid, selection with Tof
                                    iRawTimeSta =
                                        0.5 * (iTraw[dSta * nPmts][multRsta] + iTraw[dSta * nPmts + 1][multLsta]);
                                    iRawTimeSto =
                                        0.5 * (iTraw[dSto * nPmts][multRsto] + iTraw[dSto * nPmts + 1][multLsto]);
                                    iRawTof = iRawTimeSto - iRawTimeSta + iTraw[dSta * nPmts + 2][0] -
                                              iTraw[dSto * nPmts + 2][0];
                                    if ((fCalPar->GetMinTofAtRank(dSta) <= iRawTof) &&
                                        (iRawTof <= fCalPar->GetMaxTofAtRank(dSta)))
                                    {
                                        selectHit = kTRUE;

                                        selectLeftHit[dSta] = multLsta;
                                        selectLeftHit[dSto] = multLsto;

                                        selectRightHit[dSta] = multRsta;
                                        selectRightHit[dSto] = multRsto;

                                        selectLeftHit_dSto[dSta] = multLsto;
                                        selectRightHit_dSto[dSta] = multRsto;

                                        if (dSta > 0)
                                        {
                                            for (UShort_t d = dSta; d < nDets - 1; d++)
                                            {
                                                if (selectLeftHit_dSto[d] != selectLeftHit_dSto[d - 1])
                                                    selectSameHit = kFALSE;
                                                if (selectRightHit_dSto[d] != selectRightHit_dSto[d - 1])
                                                    selectSameHit = kFALSE;
                                            }
                                        }
                                    }

                                } // end of loop over the left start mult
                            }     // end of loop over the right start mult
                        }         // end of if data in Left and Right of Sta detector
                    }             // end of loop over the start detectors
                }                 // end of loop over the left stop mult
            }                     // end of loop over the right stop mult
        }                         // end of if data in Left and Right of Sto detector

        // Fill CalData levels
        Double_t iRawTime[nDets];
        Double_t iRawTime_wtref[nDets];
        Float_t iRawPos[nDets];
        Float_t iCalPos[nDets];
        for (UShort_t det = 0; det < nDets; det++)
        {
            iRawTime[det] = NAN;
            iRawTime_wtref[det] = NAN;
            iRawPos[det] = NAN;
            iCalPos[det] = NAN;
        }
        // if hit at stop detector is the same for all Sta detector
        if (selectSameHit == kTRUE)
        {
            for (UShort_t det = 0; det < nDets; det++)
            {
                iRawTime[det] =
                    0.5 * (iTraw[det * nPmts][selectRightHit[det]] + iTraw[det * nPmts + 1][selectLeftHit[det]]);
                if (nPmts >= 3 && mult[det * nPmts + 2] == 1)
                {
                    iRawTime_wtref[det] = static_cast<Float_t>(iRawTime[det] - iTraw[det * nPmts + 2][0]);
                }
                iRawPos[det] = static_cast<Float_t>(iTraw[det * nPmts][selectRightHit[det]] -
                                                    iTraw[det * nPmts + 1][selectLeftHit[det]]);
                iCalPos[det] = fCalPar->GetPosCalGainAtRank(det) * iRawPos[det] + fCalPar->GetPosCalOffsetAtRank(det);
                AddPosCalData(det + 1, iRawTime[det], iRawTime_wtref[det], iRawPos[det], iCalPos[det]);
            }
        }
        else // use the hit selection from the very first to the very last FrsSci)
        {
            for (UShort_t det = 0; det < nDets - 1; det++)
            {
                iRawTime[det] =
                    0.5 * (iTraw[det * nPmts][selectRightHit[det]] + iTraw[det * nPmts + 1][selectLeftHit[det]]);
                if (nPmts >= 3 && mult[det * nPmts + 2] == 1)
                {
                    iRawTime_wtref[det] = static_cast<Float_t>(iRawTime[det] - iTraw[det * nPmts + 2][0]);
                }
                iRawPos[det] = static_cast<Float_t>(iTraw[det * nPmts][selectRightHit[det]] -
                                                    iTraw[det * nPmts + 1][selectLeftHit[det]]);
                iCalPos[det] = fCalPar->GetPosCalGainAtRank(det) * iRawPos[det] + fCalPar->GetPosCalOffsetAtRank(det);
                AddPosCalData(det + 1, iRawTime[det], iRawTime_wtref[det], iRawPos[det], iCalPos[det]);
            }
            iRawTime[nDets - 1] = 0.5 * (iTraw[(nDets - 1) * nPmts][selectRightHit_dSto[0]] +
                                         iTraw[(nDets - 1) * nPmts + 1][selectLeftHit_dSto[0]]);
            if (nPmts >= 3 && mult[(nDets - 1) * nPmts + 2] == 1)
            {
                iRawTime_wtref[nDets - 1] =
                    static_cast<Float_t>(iRawTime[nDets - 1] - iTraw[(nDets - 1) * nPmts + 2][0]);
            }
            iRawPos[nDets - 1] = static_cast<Float_t>(iTraw[(nDets - 1) * nPmts][selectLeftHit_dSto[0]] -
                                                      iTraw[(nDets - 1) * nPmts + 1][selectLeftHit_dSto[0]]);
            iCalPos[nDets - 1] = fCalPar->GetPosCalGainAtRank(nDets - 1) * iRawPos[nDets - 1] +
                                 fCalPar->GetPosCalOffsetAtRank(nDets - 1);
            AddPosCalData(
                nDets, iRawTime[nDets - 1], iRawTime_wtref[nDets - 1], iRawPos[nDets - 1], iCalPos[nDets - 1]);
        }

        // no hit found
        if (selectHit == kFALSE)
            return;

        UShort_t rank = 0;
        for (UShort_t dSta = 0; dSta < nDets - 1; dSta++)
        {
            for (dSto = dSta + 1; dSto < nDets; dSto++)
            {
                iRawTof = iRawTime[dSto] - iRawTime[dSta] + iTraw[dSta * nPmts + 2][0] - iTraw[dSto * nPmts + 2][0];

                // Cal Tof
                iCalVelo = 1. / (fCalPar->GetTofCalOffsetAtRank(rank) + fCalPar->GetTofCalGainAtRank(rank) * iRawTof);
                iCalTof = 1. / (fCalPar->GetTofCalGainAtRank(rank) * iCalVelo);

                // Beta
                iBeta = iCalVelo / (Double_t)SPEED_OF_LIGHT_MNS;

                // Gamma
                iGamma = 1. / (sqrt(1. - pow(iBeta, 2)));

                // Brho : the signs take car of the X direction (increasing from right to left)
                iBRho = fCalPar->GetBRho0AtRank(rank) * (1 + iCalPos[dSta] / fCalPar->GetDispAtRank(dSta) -
                                                         iCalPos[dSto] / fCalPar->GetDispAtRank(dSto));
                // iBRho = fCalPar->GetBRho0AtRank(rank) * (1 + iCalPos[dSta] / fCalPar->GetDispAtRank(dSta));

                // AoQ
                iAoQ = iBRho / (3.10716 * iBeta * iGamma);
                AddTofCalData(
                    rank + 1, dSta + 1, dSto + 1, iCalPos[dSta], iCalPos[dSto], iRawTof, iCalTof, iBeta, iBRho, iAoQ);
                rank++;
            }
        }

    }    // end if else is (nDets>1)
    else // nDets==0, should not happen
    {
    }

    ++fNevent;
    return;
}

// -----   Public method Reset   ------------------------------------------------
void R3BFrsSciTcal2Cal::Reset()
{
    LOG(debug) << "Clearing TcalData Structure";
    if (fPosCal)
        fPosCal->Clear();
    if (fTofCal)
        fTofCal->Clear();
}

// -----   Private method AddPosCalData  --------------------------------------------
R3BFrsSciPosCalData* R3BFrsSciTcal2Cal::AddPosCalData(UShort_t det,
                                                      Double_t traw,
                                                      Double_t traw_wtref,
                                                      Float_t rawpos,
                                                      Float_t calpos)
{
    // It fills the R3BFrsSciTcalData
    TClonesArray& clref = *fPosCal;
    Int_t size = clref.GetEntriesFast();
    return new (clref[size]) R3BFrsSciPosCalData(det, traw, traw_wtref, rawpos, calpos);
}

// -----   Private method AddTofCalData  --------------------------------------------
R3BFrsSciTofCalData* R3BFrsSciTcal2Cal::AddTofCalData(UShort_t rank,
                                                      UShort_t detsta,
                                                      UShort_t detsto,
                                                      Float_t calpossta,
                                                      Float_t calpossto,
                                                      Double_t rawtof,
                                                      Double_t caltof,
                                                      Double_t beta,
                                                      Double_t brho,
                                                      Double_t aoq)
{
    // It fills the R3BFrsSciTcalData
    TClonesArray& clref = *fTofCal;
    Int_t size = clref.GetEntriesFast();
    return new (clref[size])
        R3BFrsSciTofCalData(rank, detsta, detsto, calpossta, calpossto, rawtof, caltof, beta, brho, aoq);
}

ClassImp(R3BFrsSciTcal2Cal);
