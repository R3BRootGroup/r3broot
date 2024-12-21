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

#ifndef R3BEVENTHEADERPROPAGATOR_H
#define R3BEVENTHEADERPROPAGATOR_H 1

#include "FairTask.h"
#include <Rtypes.h>

#include "R3BEventHeader.h"
class FairSource;

class R3BEventHeaderPropagator : public FairTask
{
  public:
    /**
     * Default constructor.
     * Creates an instance of the task with default parameters.
     */
    R3BEventHeaderPropagator();

    /**
     * Standard constructor.
     * Creates an instance of the task.
     * @param name a name of the task.
     * @param iVerbose a verbosity level.
     */
    explicit R3BEventHeaderPropagator(const TString& name,
                                      Int_t iVerbose = 1,
                                      std::string_view nameheader = "EventHeader.");

    /**
     * Method for task initialization.
     * This function is called by the framework before
     * the event loop.
     * @return Initialization status. kSUCCESS, kERROR or kFATAL.
     */
    InitStatus Init() override;

    /**
     * Method for event loop implementation.
     * Is called by the framework every time a new event is read.
     * @param option an execution option.
     */
    void Exec(Option_t* /*option*/) override;

  private:
    std::string fNameHeader;
    R3BEventHeader* fHeader = nullptr;
    FairSource* fSource = nullptr;

  public:
    ClassDefOverride(R3BEventHeaderPropagator, 1)
};

#endif // R3BEVENTHEADERPROPAGATOR_H
