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

#ifndef R3BOnlineSyncCheck_H
#define R3BOnlineSyncCheck_H 1

#include "FairTask.h"

class TClonesArray;
class R3BEventHeader;
class TCanvas;
class TH1F;
class TH2F;

class R3BOnlineSyncCheck : public FairTask
{
  public:
    R3BOnlineSyncCheck();
    R3BOnlineSyncCheck(const TString& name, Int_t iVerbose = 1);
    InitStatus Init();
    void Exec(Option_t* option);
    void FinishEvent();
    void FinishTask();
    void Reset_Histo();
    void SetTpat(Int_t tpat) { fTpat = tpat; }
    void SetTrig(Int_t trig) { fTrig = trig; }
    static TString EnumName(int val);

  private:
    enum SyncValues
    {
        MUSIC_SC,
        RPC_SC,
        S2_SC,
        FOOT1_SC,
        FOOT2_SC,
        FOOT13_SC,
        FOOT4_SC,
        FOOT11_SC,
        FOOT6_SC,
        FOOT7_SC,
        FOOT12_SC,
        NUM_SC,
        MASTER_SC,
        MASTERREF_SC
    };

    TClonesArray* fCA;
    R3BEventHeader* header;
    Int_t fTpat;
    Int_t fTrig;
    TCanvas* canvas;
    std::vector<TH2F*> fh2_array;
    static constexpr auto XBINS = 16;
    static constexpr auto YBINS = 50;
    static constexpr auto XMIN = 0;
    static constexpr auto XMAX = 16;
    static constexpr auto YMIN = 0;
    static constexpr auto YMAX = 300;
    static constexpr auto NXPADS = 4;
    static constexpr auto NYPADS = 3;

  public:
    ClassDef(R3BOnlineSyncCheck, 1)
};

#endif /* R3BOnlineSyncCheck_H */
