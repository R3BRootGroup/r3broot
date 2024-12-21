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

// -----------------------------------------------------------------------------
// -----                                                                   -----
// -----                     R3BRpcMapped2PreCal                              -----
// -----                Created 17/01/2021 by H. Alvarez-Pol               -----
// -----------------------------------------------------------------------------

#ifndef R3BRPCMAPPED2PRECAL_H
#define R3BRPCMAPPED2PRECAL_H

#include "FairTask.h"

#include "R3BRpcMapped2PreCalPar.h"

#include "R3BRpcMappedData.h"
#include "R3BTCalEngine.h"
#include "R3BTCalPar.h"

#include "R3BRpcPreCalData.h"

#include <TRandom.h>

class TClonesArray;
class R3BTCalModulePar;

class R3BRpcMapped2PreCal : public FairTask
{

  public:
    /** Default constructor **/
    R3BRpcMapped2PreCal();

    /** Destructor **/
    virtual ~R3BRpcMapped2PreCal();

    /** Virtual method Exec **/
    virtual void Exec(Option_t* option);

    /** Virtual method Reset **/
    virtual void Reset();

    // Fair specific
    /** Virtual method Init **/
    virtual InitStatus Init();

    /** Virtual method ReInit **/
    virtual InitStatus ReInit();

    /** Virtual method Finish **/
    virtual void Finish();

    void SetFpgaCorrelationFile(TString file) { fFpgaCorrelationFile = file; };

    /** Accessor to select online mode **/
    void SetOnline(Bool_t option) { fOnline = option; }

  private:
    UInt_t fNofTCalPars; /**< Number of modules in parameter file. */

    UInt_t fNofDetectors; /**< Number of detectors. */
    UInt_t fNofChannels;  /**< Number of channels per detector. */
    Double_t fClockFreq;  /**< Clock cycle in [ns]. */
    UInt_t fNEvent;
    // Don't store data for online
    Bool_t fOnline;
    TString fFpgaCorrelationFile;
    R3BTCalPar* fTCalPar;           /**< Parameter container. >*/
    std::vector<int> lut[46][2];    /**<look up table.>*/
    TClonesArray* fMappedDataCA;    /**< Array with RPC Mapped-input data. >*/
    TClonesArray* fRpcPreCalDataCA; /**< Array with RPC Cal- output data. >*/

    struct Entry
    {
        double time;
        R3BRpcMappedData const* Mapped;
    };

  public:
    // Class definition
    ClassDef(R3BRpcMapped2PreCal, 1)
};

#endif
