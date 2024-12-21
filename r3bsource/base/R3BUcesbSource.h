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

/* R3BUcesbSource.h
 * R3BROOT
 *
 * Author: Bastian Löher <b.loeher@gsi.de>
 * Date: 29/05/15
 * */

#pragma once

#include "R3BReader.h"

#include <FairSource.h>
#include <Rtypes.h>
#include <TObjArray.h>
#include <TString.h>

#include <fstream>
#include <list>

/* External data client interface (ucesb) */
#include "ext_data_clnt.hh"
#include "ext_data_struct_info.hh"

/* External data structure
 * This structure is produced using ucesb in the following way:
 * ./ucesb --ntuple=RAW,STRUCT_HH,ext_h101.h
 * */
struct EXT_STR_h101_t;
typedef struct EXT_STR_h101_t EXT_STR_h101;

class R3BEventHeader;

class R3BUcesbSource : public FairSource
{
  public:
    R3BUcesbSource(const TString&, const TString&, const TString&, EXT_STR_h101*, size_t);
    virtual ~R3BUcesbSource();

    Source_Type GetSourceType() { return kONLINE; }

    /* Init() will fork a ucesb instance to deliver the unpacked data for
     * R3BRoot. It makes use of the ext_data_ interface of ucesb.
     * */
    virtual Bool_t Init();

    virtual Bool_t SpecifyRunId();

    virtual Bool_t InitUnpackers();
    virtual void SetParUnpackers();
    virtual Bool_t ReInitUnpackers();
    /* A wrapper for ext_data_fetch_event() */
    Int_t ReadEvent(UInt_t);
    /* Close the ext_data_ interface */
    virtual void Close();
    /* Reset the readers */
    virtual void Reset();
    /* The reader interface */
    void AddReader(R3BReader* a_reader) { fReaders->Add(a_reader); }
    /* Limit the number of events */
    void SetMaxEvents(int a_max) { fLastEventNo = a_max; }
    /* Get readers */
    const TObjArray* GetReaders() const { return fReaders; }

    virtual void FillEventHeader(FairEventHeader* feh);

    inline void SetInputFileName(TString tstr) { fInputFileName = std::move(tstr); }

    inline void SetFullDataStructureValidation() { fFullDataStructureValidation = true; }

  private:
    /* File descriptor returned from popen() */
    FILE* fFd;
    /* The ucesb interface class */
    ext_data_clnt fClient;
    /* The ucesb structure info class */
    ext_data_struct_info fStructInfo;
    /* The filename to read */
    const TString fFileName;
    /* Ucesb output options */
    const TString fNtupleOptions;
    /* The location of the unpacker */
    const TString fUcesbPath;
    /* The event counter */
    unsigned int fNEvent;
    /* The full event structure */
    EXT_STR_h101* fEvent;
    size_t fEventSize;
    /* Last event requested */
    int fLastEventNo;
    /* The array of readers */
    TObjArray* fReaders;
    /* R3B header */
    R3BEventHeader* fEventHeader;
    Int_t ReadIntFromString(const std::string& wholestr, const std::string& pattern);
    TString fInputFileName;
    std::ifstream fInputFile;
    Int_t fEntryMax;
    bool fFullDataStructureValidation = false;

  public:
    // Create dictionary
    ClassDef(R3BUcesbSource, 1)
};
