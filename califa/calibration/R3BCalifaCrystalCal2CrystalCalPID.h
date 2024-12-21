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

#ifndef R3BCalifaCrystalCal2CrystalCalPID_H
#define R3BCalifaCrystalCal2CrystalCalPID_H

#include "FairTask.h"
#include "R3BCalifaCrystalCalData.h"

class TClonesArray;

/**
 * Calibrate/Normalize PID information to energy (such that E = Nf + Ns)
 *
 * @remark Needs a file called qpidscale.txt with structure
 *  channel_nr slope delta
 * with E[channel_nr] = slope*(Nf[channel_nr] + Ns[channel_nr]) + delta
 *
 * @warning Currently only 128 channels (s438b demonstrator) supported.
 *
 * @todo Support full CALIFA
 */
class R3BCalifaCrystalCal2CrystalCalPID : public FairTask
{

  public:
    /** Default constructor **/
    R3BCalifaCrystalCal2CrystalCalPID();

    /** Destructor **/
    virtual ~R3BCalifaCrystalCal2CrystalCalPID();

    /** Virtual method Exec **/
    virtual void Exec(Option_t* option);

    /** Virtual method Reset **/
    virtual void Reset();

  protected:
    // Fair specific
    /** Virtual method Init **/
    virtual InitStatus Init();

    /** Virtual method ReInit **/
    virtual InitStatus ReInit();

    /** Virtual method SetParContainers **/
    virtual void SetParContainers();

    /** Virtual method Finish **/
    virtual void Finish();

    virtual void Register();

    Long64_t nEvents;

    double scale[128];
    double delta[128];

  private:
    TClonesArray* fCrystalCalDataCA;

  public:
    // Class definition
    ClassDef(R3BCalifaCrystalCal2CrystalCalPID, 0)
};

#endif
