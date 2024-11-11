/******************************************************************************
 *   Copyright (C) 2019 GSI Helmholtzzentrum für Schwerionenforschung GmbH    *
 *   Copyright (C) 2019 Members of R3B Collaboration                          *
 *                                                                            *
 *             This software is distributed under the terms of the            *
 *                 GNU General Public Licence (GPL) version 3,                *
 *                    copied verbatim in the file "LICENSE".                  *
 *                                                                            *
 * In applying this license GSI does not waive the privileges and immunities  *
 * granted to it by virtue of its status as an Intergovernmental Organization *
 * or submit itself to any jurisdiction.                                      *
 ******************************************************************************/
// ------------------------------------------------------------
// -----                 R3BPdcCal2Hit                   -----
// -----         Created June 2nd 2020 by M.Heil          -----
// ------------------------------------------------------------

#ifndef R3BPDCCAL2HIT
#define R3BPDCCAL2HIT

#include <list>
#include <map>

#include "FairTask.h"
#include "THnSparse.h"
#define N_WIRE_MAX 144
#define N_PLANE_MAX_PDC 4
#define wire_distance 6.92820323
#define xtc_points 300

class TClonesArray;
class R3BPdcHitModulePar;
class R3BPdcHitPar;
class R3BEventHeader;
class TH1F;
class TH2F;
class R3BPdcCalData;

/**
 * An analysis task to apply HIT calibration for Pdc.
 * it also applies walk corrections and Z-calibration
 */
class R3BPdcCal2Hit : public FairTask
{
  public:
    struct ToT
    {
        ToT(R3BPdcCalData const*, R3BPdcCalData const*, Double_t, Double_t, Double_t);
        R3BPdcCalData const* lead;
        R3BPdcCalData const* trail;
        Double_t lead_ns, tail_ns, tot_ns;
    };
    struct Channel
    {
        std::list<R3BPdcCalData const*> lead_list;
        std::list<ToT> tot_list;
    };

    /**
     * Default constructor.
     * Creates an instance of the task with default parameters.
     */
    R3BPdcCal2Hit();

    /**
     * Standard constructor.
     * Creates an instance of the task.
     * @param name a name of the task.
     * @param iVerbose a verbosity level.
     */
    R3BPdcCal2Hit(const char* name, Bool_t, Bool_t, Int_t iVerbose = 1);

    /**
     * Destructor.
     * Frees the memory used by the object.
     */
    virtual ~R3BPdcCal2Hit();

    /**
     * Method for task initialization.
     * This function is called by the framework before
     * the event loop.
     * @return Initialization status. kSUCCESS, kERROR or kFATAL.
     */
    virtual InitStatus Init();

    /**
     * Method for re-initialization of parameter containers
     * in case the Run ID has changed.
     */
    virtual InitStatus ReInit();

    /**
     * Method for event loop implementation.
     * Is called by the framework every time a new event is read.
     * @param option an execution option.
     */
    virtual void Exec(Option_t* option);

    /**
     * A method for finish of processing of an event.
     * Is called by the framework for each event after executing
     * the tasks.
     */
    virtual void FinishEvent();

    /**
     * Method for finish of the task execution.
     * Is called by the framework after processing the event loop.
     */
    virtual void FinishTask();

    virtual void SetParContainers();

    R3BPdcHitModulePar* GetModuleParAt(Int_t wire);

  private:
    TClonesArray* fCalItems;        /**< Array with Cal items - input data. */
    TClonesArray* fCalTriggerItems; /**< Array with trigger Cal items - input data. */
    TClonesArray* fHitItems;        /**< Array with Hit items - output data. */
    UInt_t fNofHitItems;            /**< Number of hit items for cur event. */
    R3BPdcHitPar* fHitPar;     /**< Hit parameter container. */
	R3BPdcHitPar* fCalPar; /**< Parameter container. */
    UInt_t fNofHitPars;             /**< Number of modules in parameter file. */
    R3BEventHeader* header;         /**< Event header - input data. */
    Double_t fClockFreq;            /**< Clock cycle in [ns]. */
    UInt_t maxevent;
    UInt_t fnEvents;
	Bool_t fIsCalibrator;
	Bool_t fIsSync;
	const Double_t * xtc_x;
	const Double_t * xtc_t;
	Double_t xt_xarray[xtc_points];
    Double_t xt_tarray[xtc_points];
    Double_t t_ns[4000];
    Double_t x_mm[4000];
    Double_t rebinned_x_mm[xtc_points];
    Double_t rebinned_t_ns[xtc_points];
    Int_t ipar_channels = N_PLANE_MAX_PDC*N_WIRE_MAX ;
    Double_t tsync_mem[N_PLANE_MAX_PDC*N_WIRE_MAX];
    Bool_t method2 = true;

    typedef std::vector<Channel> ChannelArray;
    std::vector<ChannelArray> fPlaneArray;

    // histograms for gain matching
    TH2F* fh_time_raw[N_PLANE_MAX_PDC];
    TH2F* fh_time_drift[N_PLANE_MAX_PDC];
    TH1F* fh_running_sum;
    TH1F* fh_xtc;
    TH1F* fh_error[N_PLANE_MAX_PDC];
    TH2F* fh_radius[N_PLANE_MAX_PDC];
    TH2F* fh_radius_dt[N_PLANE_MAX_PDC];
    TH2F* fh_mult[N_PLANE_MAX_PDC];    
    TH2F* fh_fish[N_PLANE_MAX_PDC];
    TH2F* fh_check[N_PLANE_MAX_PDC];
    TH2F* fh_rvspos[N_PLANE_MAX_PDC];
    TH2F* fh_tsync[N_PLANE_MAX_PDC];
    TH1F* fh_angle[N_PLANE_MAX_PDC];
    
	TH2F* fh_wire_cor;
	
  public:
    ClassDef(R3BPdcCal2Hit, 3)
};

#endif
