// ------------------------------------------------------------
// -----              R3BOnlineSpectraPdc             -----
// -----            Created July-2020 by M.Heil          -----
// -----               Fill online histograms             -----
// ------------------------------------------------------------

#ifndef R3BONLINESPECTRAPDC
#define R3BONLINESPECTRAPDC

#include "FairTask.h"
#include <array>
#include <fstream>
#include <iostream>
#include <sstream>
#include <vector>

#include "TClonesArray.h"
#include "TMath.h"
#include <cstdlib>

class TClonesArray;
class TH1F;
class TH2F;
class R3BEventHeader;

#define N_PLANE_MAX_PDC 4 
/**
 * This taks reads all detector data items and plots histograms
 * for online checks.
 */
class R3BOnlineSpectraPdc : public FairTask
{

  public:
    /**
     * Default constructor.
     * Creates an instance of the task with default parameters.
     */
    R3BOnlineSpectraPdc();

    /**
     * Standard constructor.
     * Creates an instance of the task.
     * @param name a name of the task.
     * @param iVerbose a verbosity level.
     */
    R3BOnlineSpectraPdc(const char* name, Int_t iVerbose = 1);

    /**
     * Destructor.
     * Frees the memory used by the object.
     */
    virtual ~R3BOnlineSpectraPdc();

    /**
     * Method for task initialization.
     * This function is called by the framework before
     * the event loop.
     * @return Initialization status. kSUCCESS, kERROR or kFATAL.
     */
    virtual InitStatus Init();

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

    inline void SetTrigger(Int_t trigger) { fTrigger = trigger; }
    inline void SetTpat(Int_t tpat) { fTpat = tpat; }

    void Reset_PDC_Histo();
    void Reset_Fi0_Histo();
    void Reset_Fi60_Histo();

  private:
    TClonesArray* fMappedItems; /**< Array with Cal items - input data. */
    TClonesArray* fCalItems;    /**< Array with Cal items - input data. */
    TClonesArray* fHitItems;    /**< Array with Cal items - input data. */

    TClonesArray* fMappedItems_fi0; /**< Array with Cal items - input data. */
    TClonesArray* fCalItems_fi0;    /**< Array with Cal items - input data. */
    TClonesArray* fHitItems_fi0;    /**< Array with Cal items - input data. */

    TClonesArray* fMappedItems_fi60; /**< Array with Cal items - input data. */
    TClonesArray* fCalItems_fi60;    /**< Array with Cal items - input data. */
    TClonesArray* fHitItems_fi60;    /**< Array with Cal items - input data. */

    // check for trigger should be done globablly (somewhere else)
    R3BEventHeader* header; /**< Event header. */
    Int_t fTrigger;         /**< Trigger value. */
    Int_t fTpat;
    Int_t fSamp;
    Int_t fNEvents;

    TH1F* fhTpat;
    TH1F* fh_spill_length;
    TH1F* fhTrigger;

    TH1F* fh_Pdc_Wire[N_PLANE_MAX_PDC];
    TH2F* fh_Pdc_Tot[N_PLANE_MAX_PDC];
    TH2F* fh_Pdc_mult[N_PLANE_MAX_PDC];
    TH2F* fh_Pdc_ToF[N_PLANE_MAX_PDC];
    TH1F* fh_Pdc_x[N_PLANE_MAX_PDC];
    TH1F* fh_Pdc_y[N_PLANE_MAX_PDC];
    TH2F* fh_Pdc_xy[N_PLANE_MAX_PDC];
    TH2F* fh_Pdc_Time[N_PLANE_MAX_PDC];
    TH2F* fh_Pdc_Wire_vs_Events[N_PLANE_MAX_PDC];
    TH2F* fh_Pdc_ToT_vs_Time[144];
    TH2F* fh_Pdc_ToT_vs_Hit[144];
    TH2F* fh_Pdc_xvsx_fi0;
    TH2F* fh_Pdc_yvsy_fi0;
    TH2F* fh_Pdc_xvsy_fi0;
    TH2F* fh_Pdc_xvsx_fi60;
    TH2F* fh_Pdc_yvsy_fi60;
    TH2F* fh_Pdc_xvsy_fi60;
    TH2F* fh_Pdc_yvsy;
    TH2F* fh_Pdc_xvsx;
    TH2F* fh_Pdc_wvsw;
    TH2F* fh_pdc_ebene4;

    TH1F* fh_fi0_fiber;
    TH1F* fh_fi0_mult;
    TH2F* fh_fi0_Tot;
    TH2F* fh_fi0_Time;
    TH1F* fh_fi0_fiber_mc;
    TH1F* fh_fi0_mult_mc;
    TH2F* fh_fi0_Tot_mc;
    TH2F* fh_fi0_Time_mc;
    
    TH1F* fh_fi60_fiber;
    TH1F* fh_fi60_mult;
    TH2F* fh_fi60_Tot;
    TH2F* fh_fi60_Time;
    TH1F* fh_fi60_fiber_mc;
    TH1F* fh_fi60_mult_mc;
    TH2F* fh_fi60_Tot_mc;
    TH2F* fh_fi60_Time_mc;

    TH2F* fh_fi0_pdc;
    TH2F* fh_fi0_pdc_time;
    TH2F* fh_fi0_pdc_eloss;
    TH2F* fh_fi0_tot2_tot1;
    TH2F* fh_fi0_cor;
    TH2F* fh_fi0_eff;
    
    TH2F* fh_fi60_pdc;
    TH2F* fh_fi60_pdc_time;
    TH2F* fh_fi60_pdc_eloss;
    TH2F* fh_fi60_tot2_tot1;
    TH2F* fh_fi60_cor;
    TH2F* fh_fi60_eff;

  public:
    ClassDef(R3BOnlineSpectraPdc, 2)
};

#endif
