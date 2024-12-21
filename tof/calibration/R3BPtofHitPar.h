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

#ifndef R3BPTOFHITPAR_H
#define R3BPTOFHITPAR_H

#include "FairParGenericSet.h"
#include "FairParamList.h"
#include "TArrayD.h"

class R3BPtofHitPar : public FairParGenericSet
{
  public:
    R3BPtofHitPar(const char* name = "PtofHitPar",
                  const char* title = "Ptof Hit Calibration",
                  const char* context = "TestDefaultContext");
    ~R3BPtofHitPar(void){};
    void putParams(FairParamList* list);
    Bool_t getParams(FairParamList* list);

    virtual void Print(Option_t* option = "") const;
    /** Accessor functions **/

    /**
     * Method to set the parameter of a module.
     * @param paddle the paddle index [1,12]
     * @param tOffset1 the first time offset
     * @param tOffset2 the second time offset
     * @param vEff the effective speed of light
     * @param zScale the charge scale
     */
    void SetParAt(Int_t paddle, Double_t tOffset1, Double_t tOffset2, Double_t zScale, Double_t vEff)
    {
        fTOffset1[paddle - 1] = tOffset1;
        fTOffset2[paddle - 1] = tOffset2;
        fZScale[paddle - 1] = zScale;
        fVEff[paddle - 1] = vEff;
    }

    /**
     * Method to get the time offset of a module.
     * @param paddle the paddle index [1,12]
     * @param side the side {1;2}
     * @return the time offset
     */
    Double_t GetTOffsetAt(Int_t paddle, Int_t side)
    {
        if (side == 1)
            return fTOffset1[paddle - 1];
        else
            return fTOffset2[paddle - 1];
    }

    /**
     * Method to get the z-scale of a module.
     * @param paddle the paddle index [1,12]
     * @return the z-scale offset
     */
    Double_t GetZScaletAt(Int_t paddle) { return fZScale[paddle - 1]; }

    /**
     * Method to get the effective speed of light of a module.
     * @param paddle the paddle index [1,12]
     * @return the effective speed of light
     */
    Double_t GetVEfftAt(Int_t paddle) { return fVEff[paddle - 1]; }

  private:
    TArrayD fTOffset1; /**< an array with timeoffset parameters of modules */
    TArrayD fTOffset2; /**< an array with timeoffset parameters of modules */
    TArrayD fZScale;   /**< an array with zscale parameters of all modules */
    TArrayD fVEff;     /**< an array with veff parameters of all modules */

    ClassDef(R3BPtofHitPar, 1);
};

#endif /* !R3BPTOFHITPAR_H*/
