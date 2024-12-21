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

#include <iostream>
#include <sstream>
#include <string>

#include <FairLogger.h>
#include <FairRootManager.h>
#include <FairRunOnline.h>

#include "R3BEventHeader.h"
#include "R3BLogger.h"
#include "R3BUcesbSource.h"

#include "ext_data_client.h"

R3BUcesbSource::R3BUcesbSource(const TString& FileName,
                               const TString& NtupleOptions,
                               const TString& UcesbPath,
                               EXT_STR_h101* event,
                               size_t event_size)
    : FairSource()
    , fFd(nullptr)
    , fClient()
    , fStructInfo()
    , fFileName(FileName)
    , fNtupleOptions(NtupleOptions)
    , fUcesbPath(UcesbPath)
    , fNEvent(0)
    , fEvent(event)
    , fEventSize(event_size)
    , fLastEventNo(-1)
    , fEventHeader(nullptr)
    , fInputFile()
    , fEntryMax(0)
    , fReaders(new TObjArray())
{
}

R3BUcesbSource::~R3BUcesbSource()
{
    R3BLOG(debug1, "R3BUcesbSource destructor.");
    if (fReaders)
    {
        fReaders->Delete();
        delete fReaders;
    }
    R3BUcesbSource::Close();
}

Bool_t R3BUcesbSource::Init()
{
    // Register of R3BEventHeader in the output root file

    R3BLOG(debug, "checking whether R3BEventHeader has been defined in FairRun");
    auto run = FairRun::Instance();
    auto eventHeader = dynamic_cast<R3BEventHeader*>(run->GetEventHeader());
    if (eventHeader)
    {
        R3BLOG(info, "EventHeader. was defined properly");
    }
    else
    {
        eventHeader = new R3BEventHeader();
        run->SetEventHeader(eventHeader); // Implicit conversion and transfer ownership to FairRun
        R3BLOG(warn, "EventHeader. has been created from R3BEventHeader");
    }

    Bool_t status;
    std::ostringstream command;

    /* Call ucesb with this command */
    command << fUcesbPath << " " << fFileName << " "
            << "--ntuple=" << fNtupleOptions << ",STRUCT,-";

    if (fLastEventNo != -1)
    {
        command << " --max-events=" << fLastEventNo;
    }
    LOG(info) << "Calling ucesb with command: " << command.str();

    /* Fork off ucesb (calls fork() and pipe()) */
    fFd = popen(command.str().c_str(), "r");
    if (nullptr == fFd)
    {
        R3BLOG(fatal, "popen() failed");
        return kFALSE;
    }

    /* Connect to forked instance */
    status = fClient.connect(fileno(fFd));
    if (kFALSE == status)
    {
        R3BLOG(error, "ext_data_clnt::connect() failed");
        R3BLOG(fatal, "ucesb error: " << fClient.last_error());
        return kFALSE;
    }

    // Open configuration file with runid values if needed in this step
    fInputFile.open(fInputFileName.Data(), std::fstream::in);
    if (!fInputFile.is_open())
    {
        R3BLOG(warn, "Input file for RunIds could not be open, it is Ok!");
    }
    else
    {
        R3BLOG(info, "Input file for RunIds " << fInputFileName.Data() << " is open!");
        fInputFile.clear();
        fInputFile.seekg(0, std::ios::beg);
    }

    return kTRUE;
}

Bool_t R3BUcesbSource::InitUnpackers()
{
    // Register of R3BEventHeader in the output root file
    FairRootManager* frm = FairRootManager::Instance();
    R3BLOG_IF(fatal, !frm, "FairRootManager no found");

    R3BLOG(info, "Checking the register of R3BEventHeader");
    fEventHeader = dynamic_cast<R3BEventHeader*>(frm->GetObject("EventHeader."));
    if (fEventHeader)
    {
        R3BLOG(info, "EventHeader. was defined properly");
    }
    else
    {
        R3BLOG(error, "EventHeader. was not defined properly!");
    }

    /* Initialize all readers */
    for (int i = 0; i < fReaders->GetEntriesFast(); ++i)
    {
        if (!(dynamic_cast<R3BReader*>(fReaders->At(i)))->Init(&fStructInfo))
        {
            R3BLOG(fatal, "UCESB error: " << fClient.last_error());
            return kFALSE;
        }
    }

    /* Setup client */
    /* this is the version for ucesb setup with extended mapping info */
    uint32_t struct_map_success = 0;
    Bool_t status = fClient.setup(NULL, 0, &fStructInfo, &struct_map_success, fEventSize);
    if (status != 0)
    {
        // perror("ext_data_clnt::setup()");
        R3BLOG(error, "ext_data_clnt::setup() failed");
        R3BLOG(fatal, "UCESB error: " << fClient.last_error());
        return kFALSE;
    }
    /*
     * It is not needed, that *all* items are matched.
     * However, mapping should fail, if items are requested that don't exist
     * on the server, or if items are requested with wrong parameters.
     * See ucesb/hbook/ext_data_client.h for more information.
     */
    uint32_t map_ok = EXT_DATA_ITEM_MAP_OK | EXT_DATA_ITEM_MAP_NO_DEST;
    if (struct_map_success & ~(map_ok))
    {
        R3BLOG(warn, "ext_data_clnt::setup() failed to map all items:");
        ext_data_struct_info_print_map_success(fStructInfo, stderr, map_ok);
        /* FairRunOnline::Init() ignores the return value from
         * GetSource()->InitUnpackers(); so do a FATAL error.
         */
        if (fFullDataStructureValidation)
        {
            R3BLOG(fatal,
                   "ext_data_clnt::setup() mapping failure may "
                   "cause unexpected analysis results due to missing "
                   "data members. Unpacker needs fixing.");
        }
        else
        {
            R3BLOG(error,
                   "ext_data_clnt::setup() mapping failure may "
                   "cause unexpected analysis results due to missing "
                   "data members. Unpacker needs fixing.");
        }

        return kFALSE;
    }

    return kTRUE;
}

void R3BUcesbSource::SetParUnpackers()
{
    for (int i = 0; i < fReaders->GetEntriesFast(); ++i)
    {
        (dynamic_cast<R3BReader*>(fReaders->At(i)))->SetParContainers();
    }
}

Bool_t R3BUcesbSource::ReInitUnpackers()
{
    /* Initialize all readers */
    for (int i = 0; i < fReaders->GetEntriesFast(); ++i)
    {
        if (!(dynamic_cast<R3BReader*>(fReaders->At(i)))->ReInit())
        {
            R3BLOG(fatal, "ReInit of a reader failed.");
            return kFALSE;
        }
    }

    return kTRUE;
}

Int_t R3BUcesbSource::ReadIntFromString(const std::string& wholestr, const std::string& pattern)
{
    std::string tempstr = wholestr;
    tempstr.replace(0, tempstr.find(pattern) + pattern.length(), "");
    tempstr.replace(0, tempstr.find('=') + 1, "");
    return atoi(tempstr.c_str());
}

Int_t R3BUcesbSource::ReadEvent(UInt_t i)
{
    const void* raw;
    ssize_t raw_words;
    int ret;
    (void)i; /* Why is i not used? Outer loop seems not to use it. */

    R3BLOG(debug1, "ReadEvent " << fNEvent);

    fNEvent++;

    if (fNEvent > fEntryMax && fEntryMax != -1 && fInputFile.is_open())
    {

        R3BLOG(info, "ReadEvent()");

        std::string buffer;
        do
        {
            getline(fInputFile, buffer);
            LOG(info) << "read from file: \"" << buffer << "\"";
            if (buffer.find("EVENT BEGIN") == 0)
            {
                fRunId = ReadIntFromString(buffer, "RUNID");
                fEventHeader->SetRunId(fRunId);
            }
            if (buffer.find("EVENT") == 0)
                continue;
            Int_t fInit = atoi(buffer.c_str());
            buffer.erase(0, buffer.find(' ') + 1);
            fEntryMax = atoi(buffer.c_str());

        } while (fInputFile && buffer.compare("EVENT END"));
    }

    /* Need to initialize first */
    if (nullptr == fFd)
    {
        Init();
    }

    /* Fetch data */
    ret = fClient.fetch_event(fEvent, fEventSize);
    if (0 == ret)
    {
        LOG(info) << "End of input";
        return 1;
    }
    if (-1 == ret)
    {
        perror("ext_data_clnt::fetch_event()");
        R3BLOG(error, "ext_data_clnt::fetch_event() failed");
        R3BLOG(fatal, "UCESB error: " << fClient.last_error());
        return 0;
    }

    /* Get raw data, if any */
    ret = fClient.get_raw_data(&raw, &raw_words);
    if (0 != ret)
    {
        perror("ext_data_clnt::get_raw_data()");
        R3BLOG(fatal, "Failed to get raw data.");
        return 0;
    }

    /* Run detector specific readers */
    for (int r = 0; r < fReaders->GetEntriesFast(); ++r)
    {
        R3BReader* reader = dynamic_cast<R3BReader*>(fReaders->At(r));

        LOG(debug1) << "  Reading reader " << r << " (" << reader->GetName() << ")";
        reader->R3BRead();
    }

    /* Display raw data */
    if (raw)
    {
        int w, j;
        const uint32_t* u = reinterpret_cast<const uint32_t*>(raw);

        LOG(info) << "  Raw data:";
        for (w = 0; w < raw_words; w += 8)
        {
            printf("    RAW%4x:", w);
            for (j = 0; j < 8 && w + j < raw_words; j++)
                printf(" %08x", u[w + j]);
            printf("\n");
        }
    }

    return 0;
}

void R3BUcesbSource::Close()
{
    int ret;

    if (!fFd) // not open
        return;
    /* Close client connection */
    ret = fClient.close();
    if (0 != ret)
    {
        // perror("ext_data_clnt::close()");
        R3BLOG(fatal, "ext_data_clnt::close() failed");
    }

    /* Close pipe */
    int status;
    status = pclose(fFd);
    if (-1 == status)
    {
        // perror("pclose()");
        R3BLOG(fatal, "pclose() failed");
        abort();
    }
    fFd = nullptr;

    if (fInputFile.is_open())
        fInputFile.close();
}

void R3BUcesbSource::Reset()
{
    for (int i = 0; i < fReaders->GetEntriesFast(); ++i)
    {
        (dynamic_cast<R3BReader*>(fReaders->At(i)))->Reset();
    }
}

Bool_t R3BUcesbSource::SpecifyRunId()
{
    if (ReadEvent(0) == 0)
        return true;
    return false;
}

//_____________________________________________________________________________
void R3BUcesbSource::FillEventHeader(FairEventHeader* feh)
{
    // note: as of 2023-02-09, FairRootManager will always pass us a pointer to
    // FairEventHeader, never to R3BEventHeader. Investigation is under way.
    //    --Philipp
    feh->SetRunId(fRunId); // we can still set the run ID, though.
}

ClassImp(R3BUcesbSource);
