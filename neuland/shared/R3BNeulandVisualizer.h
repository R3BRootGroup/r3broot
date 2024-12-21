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

#ifndef R3BNEULANDVISUALIZER_H
#define R3BNEULANDVISUALIZER_H 1

#include "TCanvas.h"
#include "TFile.h"
#include "TH3D.h"
#include "TTree.h"
#include <memory>

class R3BNeulandVisualizer
{
  public:
    R3BNeulandVisualizer(const TString& input_file, const TString& what);
    virtual ~R3BNeulandVisualizer() = default; // to avoid a ROOT-based (?) warning

    inline void ShowEvent(const UInt_t i)
    {
        fIndex = i;
        Visualize();
    }

    inline void ShowNext()
    {
        fIndex++;
        Visualize();
    }

  protected:
    void Visualize();

  private:
    std::shared_ptr<TFile> fFile;
    TTree* fTree; // non-owning
    TH3D* fh3;    // non-owning
    std::shared_ptr<TCanvas> fCanvas;
    UInt_t fIndex;

    ClassDef(R3BNeulandVisualizer, 0);
};

#endif // R3BNEULANDVISUALIZER_H
