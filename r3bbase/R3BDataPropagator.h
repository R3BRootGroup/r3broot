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

#ifndef R3BDataPropagator_H
#define R3BDataPropagator_H

#include "FairTask.h"
#include <Rtypes.h>

class TClonesArray;

class R3BDataPropagator : public FairTask
{
  public:
    /**
     * Default constructor.
     * Creates an instance of the task with default parameters.
     */
    R3BDataPropagator();

    /**
     * Standard constructor.
     * Creates an instance of the task.
     * @param name a name of the task.
     * @param iVerbose a verbosity level.
     */
    R3BDataPropagator(const TString& name, Int_t iVerbose = 1, const TString& namebranch = "");

    /**
     * Destructor.
     * Frees the memory used by the object.
     */
    virtual ~R3BDataPropagator();

    /**
     * Method for task initialization.
     * This function is called by the framework before
     * the event loop.
     * @return Initialization status. kSUCCESS, kERROR or kFATAL.
     */
    virtual InitStatus Init() override;

    /**
     * Method for event loop implementation.
     * Is called by the framework every time a new event is read.
     * @param option an execution option.
     */
    virtual void Exec(Option_t* option) override;

  private:
    TString fNameBranch;
    TClonesArray* fInputData; // Array with input data

  public:
    ClassDefOverride(R3BDataPropagator, 0)
};

#endif // R3BDataPropagator_H
