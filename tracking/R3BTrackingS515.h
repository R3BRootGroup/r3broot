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

// Created on 8/11/2022 by V.Panin

#ifndef R3BTRACKINGS515
#define R3BTRACKINGS515

#include "FairTask.h"
#include "R3BEventHeader.h"
#include "R3BMDFWrapper.h"
#include "R3BTrack.h"
#include "TCanvas.h"
#include "TClonesArray.h"
#include "TMath.h"
#include "TString.h"
#include "TVector3.h"
#include <array>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <sstream>

class R3BTrackingS515 : public FairTask
{
  public:
    R3BTrackingS515();
    R3BTrackingS515(const char* name, Int_t iVerbose = 1);
    virtual ~R3BTrackingS515();
    virtual InitStatus Init();
    virtual void Exec(Option_t* option);
    virtual void FinishEvent();
    virtual void FinishTask();

    // GLAD current which was used for training MDF function
    // Should  be set form the steering macro
    void SetGladReferenceCurrent(Double_t cur) { GladReferenceCurrent = cur; }

    // GLAD current in the run being anlaysed, set from steering macro
    void SetGladCurrent(Double_t cur) { GladCurrent = cur; }

    // Set MDF functions from the steering macro
    void Set_MDF_PoQ(TString name) { MDF_PoQ_filename = name; }
    void Set_MDF_FlightPath(TString name) { MDF_FlightPath_filename = name; }
    void Set_MDF_TX0(TString name) { MDF_TX0_filename = name; }
    void Set_MDF_TY0(TString name) { MDF_TY0_filename = name; }
    void Set_MDF_TX1(TString name) { MDF_TX1_filename = name; }
    void Set_MDF_TY1(TString name) { MDF_TY1_filename = name; }
    void SetTrigger(Int_t trigger) { fTrigger = trigger; }
    void SetTpat(Int_t tpat) { fTpat = tpat; }
    void SetMaxEvent(Int_t nev) { maxevent = nev; }

    // Set lab positions and angles of the detectors from the steering macro
    void SetPositionLOS(double x, double y, double z) { los_position.SetXYZ(x, y, z); }   // cm
    void SetPositionMWPC(double x, double y, double z) { mwpc_position.SetXYZ(x, y, z); } // cm
    void SetPositionF10(double x, double y, double z) { f10_position.SetXYZ(x, y, z); }   // cm
    void SetPositionF11(double x, double y, double z) { f11_position.SetXYZ(x, y, z); }   // cm
    void SetPositionF12(double x, double y, double z) { f12_position.SetXYZ(x, y, z); }   // cm
    void SetPositionTofd(double x, double y, double z) { tofd_position.SetXYZ(x, y, z); } // cm

    void SetEulerAnglesLOS(double x, double y, double z) { los_angles.SetXYZ(x, y, z); }   // rad
    void SetEulerAnglesMWPC(double x, double y, double z) { mwpc_angles.SetXYZ(x, y, z); } // rad
    void SetEulerAnglesF10(double x, double y, double z) { f10_angles.SetXYZ(x, y, z); }   // rad
    void SetEulerAnglesF11(double x, double y, double z) { f11_angles.SetXYZ(x, y, z); }   // rad
    void SetEulerAnglesF12(double x, double y, double z) { f12_angles.SetXYZ(x, y, z); }   // rad
    void SetEulerAnglesTofd(double x, double y, double z) { tofd_angles.SetXYZ(x, y, z); } // rad

    // Get lab positions and angles (needed by alignment function)
    inline TVector3 GetPositionLOS() { return los_position; }   // cm
    inline TVector3 GetPositionMWPC() { return mwpc_position; } // cm
    inline TVector3 GetPositionF10() { return f10_position; }   // cm
    inline TVector3 GetPositionF11() { return f11_position; }   // cm
    inline TVector3 GetPositionF12() { return f12_position; }   // cm
    inline TVector3 GetPositionTofd() { return tofd_position; } // cm

    inline TVector3 GetEulerAnglesLOS() { return los_angles; }   // rad
    inline TVector3 GetEulerAnglesMWPC() { return mwpc_angles; } // rad
    inline TVector3 GetEulerAnglesF10() { return f10_angles; }   // rad
    inline TVector3 GetEulerAnglesF11() { return f11_angles; }   // rad
    inline TVector3 GetEulerAnglesF12() { return f12_angles; }   // rad
    inline TVector3 GetEulerAnglesTofd() { return tofd_angles; } // rad

    R3BMDFWrapper* Get_MDF_PoQ() { return MDF_PoQ; }

    void SetTofOffset(double offset) { tof_offset = offset; } // ns

    // Transofrming input detector hit (point) into laboratory system
    void TransformPoint(TVector3& point, TVector3 rotation, TVector3 translation);

    // Setup incoming cuts on FRS data and Music data
    void SetFrsEllipticPIDCut(double meanZ, double rZ, double meanAoZ, double rAoZ);
    bool IsInsideFrsEllipticPIDCut(double z, double aoz);
    void SetMusicZMinMax(double min, double max)
    {
        MusicZMin = min;
        MusicZMax = max;
    }
    void SetFRSBrhoMinMax(double min, double max)
    {
        FrsBrhoMin = min;
        FrsBrhoMax = max;
    }
    void SetFiberEnergyMinMax(double min, double max)
    {
        FiberEnergyMin = min;
        FiberEnergyMax = max;
    }

    // Setters for the alignment procedure
    void SetDoAlignment(bool flag) { DoAlignment = flag; }
    void SetReferencePoQ(Double_t val) { reference_PoQ = val; }
    Double_t GetReferencePoQ() { return reference_PoQ; }

    void Alignment();
    static double AlignmentErrorS515(const double* par);

    // Data containesr needed only for the Alignment() function
    struct det_points
    {
        TVector3 mwpc;
        TVector3 f10;
        TVector3 f11;
        TVector3 f12;
        TVector3 tofd;
    };
    std::vector<det_points> det_points_vec;

    TVector3 mwpc_point, f10_point, f11_point, f12_point, tofd_point;
    TVector3 mwpc_ang_offset, f10_ang_offset, f11_ang_offset, f12_ang_offset;
    TVector3 mwpc_pos_offset, f10_pos_offset, f11_pos_offset, f12_pos_offset;

  private:
    //-- Input hit data from the TClonesArray
    // do not change the order, add new det in the end
    enum DetectorInstances
    {
        DET_FI_FIRST,
        DET_FI10 = DET_FI_FIRST,
        DET_FI11,
        DET_FI12,
        DET_FI_LAST = DET_FI12,
        DET_TOFD,
        MWPC0_HITDATA,
        FRS_DATA,
        MUSIC_HITDATA,
        LOS_HITDATA,
        DET_MAX
    };

#define NOF_FIB_DET (DET_FI_LAST - DET_FI_FIRST + 1)

    // Names of essential branches in the input tree
    // do not change the order! add new data in the end
    const char* fDetectorNames[DET_MAX + 1] = { "Fi10Hit", "Fi11Hit",      "Fi12Hit", "TofdHit", "Mwpc0HitData",
                                                "FrsData", "MusicHitData", "LosHit",  NULL };

    R3BEventHeader* fHeader;
    std::vector<TClonesArray*> fDataItems; // input data
    TClonesArray* fTrackItems;             // output data

    TVector3 los_position;
    TVector3 mwpc_position;
    TVector3 f10_position;
    TVector3 f11_position;
    TVector3 f12_position;
    TVector3 tofd_position;

    TVector3 los_angles;
    TVector3 mwpc_angles;
    TVector3 f10_angles;
    TVector3 f11_angles;
    TVector3 f12_angles;
    TVector3 tofd_angles;

    R3BMDFWrapper* MDF_FlightPath;
    R3BMDFWrapper* MDF_PoQ;
    R3BMDFWrapper* MDF_TX0;
    R3BMDFWrapper* MDF_TX1;
    R3BMDFWrapper* MDF_TY0;
    R3BMDFWrapper* MDF_TY1;

    TString MDF_FlightPath_filename;
    TString MDF_PoQ_filename;
    TString MDF_TX0_filename;
    TString MDF_TY0_filename;
    TString MDF_TX1_filename;
    TString MDF_TY1_filename;

    Double_t mdf_data[8];   // data container for the MDF function
    unsigned long fNEvents; // Event counter
    Int_t fTrigger;
    Int_t fTpat;
    Int_t maxevent;
    Double_t GladCurrent;
    Double_t GladReferenceCurrent;
    Double_t reference_PoQ;
    Bool_t DoAlignment;
    Double_t tof_offset; // ns

    // frs eliptic cuts:
    Double_t FrsCutMeanZ;
    Double_t FrsCutMeanAoZ;
    Double_t FrsCutRadiusZ;
    Double_t FrsCutRadiusAoZ;

    // frs eliptic cuts:
    Double_t FrsBrhoMin;
    Double_t FrsBrhoMax;

    // Z cut ranges in MUSIC set by SetMusicZMinMax()
    Double_t MusicZMin;
    Double_t MusicZMax;

    // Cut on fiber hit energy set by SetFiberEnergyMinMax():
    Double_t FiberEnergyMin;
    Double_t FiberEnergyMax;

    Double_t PoQ;
    Double_t FlightPath;
    Double_t ToF;
    Double_t TX0;
    Double_t TX1;
    Double_t TY0;
    Double_t TY1;
    Double_t Beta;
    Double_t Gamma;
    Double_t mdf_AoZ;

    // Essential constants
    const Double_t SPEED_OF_LIGHT = 29.9792458; // cm/ns
    const Double_t AMU = 0.9314940038;          // GeV/c2

    // Private method to fill output track data
    R3BTrack* AddTrackData(TVector3 mw, TVector3 poq, Double_t charge, Double_t aoz);

  public:
    ClassDef(R3BTrackingS515, 1)
};

#endif
