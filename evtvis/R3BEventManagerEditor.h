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

#include <TGedFrame.h>

class R3BEventManager;
class TGNumberEntry;
class TGCheckButton;
class TEveGValuator;

class R3BEventManagerEditor : public TGedFrame
{
    R3BEventManagerEditor(const R3BEventManagerEditor&);            // Not implemented
    R3BEventManagerEditor& operator=(const R3BEventManagerEditor&); // Not implemented

  protected:
    TObject* fObject;
    R3BEventManager* fManager;
    TGNumberEntry *fCurrentEvent, *fCurrentPDG;
    TGCheckButton* fVizPri;
    TGCheckButton* fScaleE;
    TEveGValuator *fMinEnergy, *fMaxEnergy;

  public:
    R3BEventManagerEditor(const TGWindow* p = 0,
                          Int_t width = 170,
                          Int_t height = 30,
                          UInt_t options = kChildFrame,
                          Pixel_t back = GetDefaultFrameBackground());
    virtual ~R3BEventManagerEditor() {}
    void SetModel(TObject* obj);
    virtual void SelectEvent();
    virtual void SelectPDG();
    void DoVizPri();
    virtual void MaxEnergy();
    virtual void MinEnergy();
    virtual void Init();
    virtual void DoScaleE();
    ClassDef(R3BEventManagerEditor, 0); // Specialization of TGedEditor for proper update propagation to TEveManager.
};
