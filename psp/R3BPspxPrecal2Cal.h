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
// -----                    R3BPspxPrecal2Cal                 -----
// -----            Created  20-03-2017 by I. Syndikus	-----
// -----              Modified  Dec 2019  by M. Holl	        -----
// ----------------------------------------------------------------

#ifndef R3BPSPXPRECAL2CAL_H
#define R3BPSPXPRECAL2CAL_H

#include "FairTask.h"

class TClonesArray;
class R3BEventHeader;
class R3BPspxCalPar;

/**
 * Class to convert Precal data to Cal data for PSPX detector data.
 * Gains for energy matching of all strips are applied
 * Total deposited energt and position in internal coordinates is calculated
 * @author Ina Syndikus
 * @since March 20, 2016
 * Modified Dec 2019 by M.Holl
 * Modified April 2021 by J.L.Rodriguez
 */

class R3BPspxPrecal2Cal : public FairTask
{
  public:
    /** Default Constructor **/
    R3BPspxPrecal2Cal();

    /** Standard Constructor **/
    R3BPspxPrecal2Cal(const char* name, Int_t iVerbose);

    /** Destructor **/
    virtual ~R3BPspxPrecal2Cal();

    // Fair specific
    /** Virtual method Init **/
    virtual InitStatus Init();

    /** Virtual method ReInit **/
    virtual InitStatus ReInit();

    /** Virtual method Exec **/
    virtual void Exec(Option_t* option);

    /** Virtual method FinishEvent **/
    virtual void FinishEvent();

    /** Virtual method FinishTask **/
    virtual void FinishTask();

    /** Method SetParContainers **/
    void SetParContainers();

    /** Accessor to select online mode **/
    void SetOnline(Bool_t option) { fOnline = option; }

  private:
    void SetParameters();

    R3BEventHeader* fHeader;                 // do we need that?
    std::vector<TClonesArray*> fPrecalItems; /**< Arrays holding input (Precal) data */
    std::vector<TClonesArray*> fCalItems;    /**< Arrays holding output (Cal) data */

    R3BPspxCalPar* fCalPar; /**< Parameter instance holding gain for energy correction */
    std::vector<std::vector<Float_t>> gain;

    Bool_t fOnline; // Don't store data for online

  public:
    ClassDef(R3BPspxPrecal2Cal, 4)
};

#endif
