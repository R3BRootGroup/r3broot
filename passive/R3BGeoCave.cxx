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

#include "R3BGeoCave.h"

#include "FairGeoBasicShape.h"
#include "FairGeoMedia.h"
#include "FairGeoMedium.h"
#include "FairGeoNode.h"
#include "FairGeoShapes.h"
#include "FairLogger.h"

#include <TList.h>
#include <iostream>
#include <string.h>

using namespace std;

R3BGeoCave::R3BGeoCave()
    : FairGeoSet()
    , name("cave")
{
    // Constructor
    fName = "cave";
    maxModules = 1;
}

Bool_t R3BGeoCave::read(fstream& fin, FairGeoMedia* media)
{
    // Reads the geometry from file
    if (!media)
    {
        return kFALSE;
    }
    const Int_t maxbuf = 256;
    char buf[maxbuf];
    FairGeoNode* volu = 0;
    FairGeoMedium* medium;
    Bool_t rc = kTRUE;
    do
    {
        fin.getline(buf, maxbuf);
        if (buf[0] != '\0' && buf[0] != '/' && !fin.eof())
        {
            if (strcmp(buf, name) == 0)
            {
                volu = new FairGeoNode;
                volu->SetName(buf);
                volu->setVolumeType(kFairGeoTopNode);
                volu->setActive();
                fin.getline(buf, maxbuf);
                TString shape(buf);
                FairGeoBasicShape* sh = pShapes->selectShape(shape);
                if (sh)
                {
                    volu->setShape(sh);
                }
                else
                {
                    rc = kFALSE;
                }
                fin.getline(buf, maxbuf);
                medium = media->getMedium(buf);
                if (!medium)
                {
                    medium = new FairGeoMedium();
                    media->addMedium(medium);
                }
                volu->setMedium(medium);
                Int_t n = 0;
                if (sh)
                {
                    n = sh->readPoints(&fin, volu);
                }
                if (n <= 0)
                {
                    rc = kFALSE;
                }
            }
            else
            {
                rc = kFALSE;
            }
        }
    } while (rc && !volu && !fin.eof());
    if (volu && rc)
    {
        volumes->Add(volu);
        masterNodes->Add(new FairGeoNode(*volu));
    }
    else
    {
        delete volu;
        volu = 0;
        rc = kFALSE;
    }
    return rc;
}

void R3BGeoCave::addRefNodes()
{
    // Adds the reference node
    FairGeoNode* volu = getVolume(name);
    if (volu)
    {
        masterNodes->Add(new FairGeoNode(*volu));
    }
}

void R3BGeoCave::write(fstream& fout)
{
    // Writes the geometry to file
    fout.setf(ios::fixed, ios::floatfield);
    FairGeoNode* volu = getVolume(name);
    if (volu)
    {
        FairGeoBasicShape* sh = volu->getShapePointer();
        FairGeoMedium* med = volu->getMedium();
        if (sh && med)
        {
            fout << volu->GetName() << '\n' << sh->GetName() << '\n' << med->GetName() << '\n';
            sh->writePoints(&fout, volu);
        }
    }
}

void R3BGeoCave::print()
{
    // Prints the geometry
    FairGeoNode* volu = getVolume(name);
    if (volu)
    {
        FairGeoBasicShape* sh = volu->getShapePointer();
        FairGeoMedium* med = volu->getMedium();
        if (sh && med)
        {
            LOG(info) << volu->GetName() << '\n' << sh->GetName() << '\n' << med->GetName() << '\n';
            sh->printPoints(volu);
        }
    }
}

ClassImp(R3BGeoCave);
