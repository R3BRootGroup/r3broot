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

#ifndef R3BSci2Tcal2Hit_H
#define R3BSci2Tcal2Hit_H

#include "FairTask.h"
#include "R3BSci2HitData.h"
class TClonesArray;
class R3BSci2HitPar;

class R3BSci2Tcal2Hit : public FairTask
{

  public:
    /*
     * Default constructor.
     * Creates an instance of the task with default parameters.
     */
    R3BSci2Tcal2Hit();

    /*
     * Standard constructor.
     * Creates an instance of the task.
     * @param name a name of the task.
     * @param iVerbose a verbosity level.
     */
    R3BSci2Tcal2Hit(const char* name, Int_t iVerbose = 1);

    /*
     * Destructor.
     * Frees the memory used by the object.
     */
    virtual ~R3BSci2Tcal2Hit();

    /*
     * Method for task initialization.
     * This function is called by the framework before
     * the event loop.
     * @return Initialization status. kSUCCESS, kERROR or kFATAL.
     */
    virtual InitStatus Init();

    /*
     * Method for re-initialization of parameter containers
     * in case the Run ID has changed.
     */
    virtual InitStatus ReInit();

    /*
     * Method for event loop implementation.
     * Is called by the framework every time a new event is read.
     * @param option an execution option.
     */
    virtual void Exec(Option_t* option);

    /*
     * A method for finish of processing of an event.
     * Is called by the framework for each event after executing
     * the tasks.
     */
    virtual void FinishEvent();

    /*
     * Methods for setting position offset and effective velocity of light
     */
    inline void SetSci2ParamMCFD(Double_t offsetX, Double_t veffX)
    {
        fsci2OffsetX = offsetX;
        fsci2VeffX = veffX;
    }

    inline void SetSci2ParamTAMEX(Double_t offsetXT, Double_t veffXT)
    {
        fsci2OffsetXT = offsetXT;
        fsci2VeffXT = veffXT;
    }

    // Method to select online mode
    void SetOnline(Bool_t option) { fOnline = option; }

    // Method to set Coincidence window
    void SetCoincWindow(Double_t t_Diff) { fCoincWindow = t_Diff; }

    virtual void SetParContainers();

  private:
    void SetParameter();
    R3BSci2HitPar* fSci2Hit_Par; // Parameter container
    Bool_t fOnline;              // Don't store data for online
    TClonesArray* fCalItems;     /* < Array with Cal items - input data. */
    TClonesArray* fHitItems;     /* < Array with Hit items - output data. */
    Double_t fClockFreq;         /* < Clock cycle in [ns]. */
    Double_t fsci2VeffX;
    Double_t fsci2OffsetX;
    Double_t fsci2OffsetXT;
    Double_t fsci2VeffXT;
    Double_t fCoincWindow;
    Float_t fPos_p0, fPos_p1;

    Int_t Icount = 0;

    /** Private method AddHitData **/
    R3BSci2HitData* AddHitData(Int_t sci, Double_t x, Double_t Tmean, Double_t Tmean_w_Tref);

  public:
    ClassDef(R3BSci2Tcal2Hit, 1)
};

#endif /* R3BSci2Tcal2Hit_H */
