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

#include <FairTask.h>
#include <TClonesArray.h>
#include <list>

class TH1F;
class TH2F;
class R3BEventHeader;
class R3BCoarseTimeStitch;
class R3BFiberMAPMTCalData;
class R3BFiberMAPMTHitPar;
class R3BFiberMappingPar;
class R3BFiberMAPMTHitModulePar;

class R3BFiberMAPMTCal2Hit : public FairTask
{
  public:
    // This defines the fiber direction
    enum Direction
    {
        HORIZONTAL,
        VERTICAL
    };
    // This tells us where the first fiber is located:
    // VERTICAL direction ---------------------------------
    // STANDARD: fiber 1 wixhausen -> fiber 512 messel
    // INVERTED: fiber 1 messel -> fiber 512 wixhausen
    // HORIZONTAL direction -------------------------------
    // STANDARD: fiber 1 down -> fiber 512 up
    // INVERTED: fiber 1 up -> fiber 512 down
    enum Orientation
    {
        STANDARD,
        INVERTED
    };

    struct ToT
    {
        ToT(R3BFiberMAPMTCalData const*, R3BFiberMAPMTCalData const*, Double_t, Double_t, Double_t);
        R3BFiberMAPMTCalData const* lead;
        R3BFiberMAPMTCalData const* trail;
        Double_t lead_ns, tail_ns, tot_ns;
    };
    struct Channel
    {
        std::list<R3BFiberMAPMTCalData const*> lead_list;
        std::list<ToT> tot_list;
    };

    struct Fib_Hit
    {
        void Reset();
        Double_t time[2][2]; //[0-bottom, 1-top][0-lead, 1 -trail]
        int channel;
    };

    R3BFiberMAPMTCal2Hit(const char*,
                         Int_t v = 1, // verbosity
                         Direction dir = HORIZONTAL,
                         UInt_t fnb = 512,    // fiber number
                         Bool_t iscal = false // is calib
    );

    virtual ~R3BFiberMAPMTCal2Hit();

    virtual InitStatus Init();

    virtual InitStatus ReInit();

    virtual void SetParContainers();

    virtual void Exec(Option_t*);

    virtual void FinishEvent();

    virtual void FinishTask();

    inline void SetTofWindow(Double_t tmin, Double_t tmax)
    {
        ftofmin = tmin;
        ftofmax = tmax;
    }

    inline void SetDTimeWindow(double dtwin) { fDTime_window = dtwin; }

    inline void SetWriteHisto() { fWrite = kTRUE; }

    inline void SetTofMin(Double_t min) { ftofmin = min; }

    inline void SetTofMax(Double_t max) { ftofmax = max; }

    inline void SetGate(Double_t g) { fGate_ns = g; }

    inline void SetOrientation(Orientation opt) { fOrientation = opt; }

    // Accessor to select online mode
    inline void SetOnline(Bool_t option) { fOnline = option; }

  private:
    TString fName;
    Int_t fDetId;
    Int_t fnEvents;
    Int_t fNumFibers;
    Int_t maxevent;
    Int_t multi;
    Double_t fClockFreq;
    Double_t fGate_ns;

    Double_t tsync;
    Double_t gainUp;
    Double_t gainDown;
    // Double_t offsetUp;
    Double_t offsetDT;
    Bool_t fIsCalibrator;
    Double_t ftofmin, ftofmax;
    Bool_t fWrite;
    // Don't store data for online
    Bool_t fOnline;
    double fDTime_window = 20.; // ns

    R3BEventHeader* fHeader; /* Event header  */
    R3BCoarseTimeStitch* fTimeStitch;
    Direction fDirection;
    Orientation fOrientation;
    TClonesArray* fCalItems;
    TClonesArray* fCalTriggerItems;
    TClonesArray* fHitItems;
    R3BFiberMappingPar* fMapPar;
    R3BFiberMAPMTHitPar* fCalPar; /**< Parameter container. */
    R3BFiberMAPMTHitPar* fHitPar; /**< Hit parameter container. */
    Int_t fNofHitPars;            /**< Number of modules in parameter file. */
    // [0=bottom,1=top][Channel].
    std::vector<Channel> fChannelArray[2];

    // histograms for gain matching
    TH2F* fh_ToT_bottom_Fib_raw;
    TH2F* fh_ToT_top_Fib_raw;
    TH2F* fh_ToT_bottom_Fib;
    TH2F* fh_ToT_top_Fib;
    TH2F* fh_time_top_Fib;
    TH2F* fh_dt_Fib;
    TH2F* fh_Fib_ToF;
    TH2F* fh_dt_Fib_raw;
    TH2F* fh_Fib_ToF_raw;
    TH2F* fh_time_check_tsync;
    TH1F* fh_multi;
    TH2F* fh_time_bottom_Fib;

  public:
    ClassDef(R3BFiberMAPMTCal2Hit, 3)
};
