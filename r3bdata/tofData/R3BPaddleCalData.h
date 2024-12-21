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

#ifndef R3BPADDLECALITEM_H
#define R3BPADDLECALITEM_H

#include "TObject.h"
#include <fstream>
#include <iostream>
#include <sstream>
using namespace std;
// for the data analysis of the paddle-type detectors.
// Introduced by Ralf, Apr 2016

class R3BPaddleCalData : public TObject
{
  public:
    // Default Constructor
    R3BPaddleCalData();

    // Standard Constructor
    R3BPaddleCalData(UInt_t plane, UInt_t bar); //, Double_t r,Double_t t,Double_t l,Double_t b,Double_t tref);

    // Destructor
    virtual ~R3BPaddleCalData() {}

    // Getters
    // inline const UChar_t&  GetDetector()  const { return fDetector; }
    UInt_t GetPlane() const { return fPlane; }
    UInt_t GetBar() const { return fBar; }

    // tube: 0,1 edge: 0,1 time: in ns
    void SetTime(int tube, int edge, Double_t t)
    {
        if (tube == 0)
        {
            if (edge == 0)
                fTime1L_ns = t;
            else
                fTime1T_ns = t;
        }
        else
        {
            if (edge == 0)
                fTime2L_ns = t;
            else
                fTime2T_ns = t;
        }
    }

    void SetTimeL(int tube, int edge, Double_t t)
    {
        if (tube == 0)
            fTime1L_ns = t;
        else
            fTime2L_ns = t;
        // cout<< "PaddleCalData lead "<<tube<<", "<<edge<<", "<<fTime1L_ns<<"; "<<fTime2L_ns<<endl;
    }

    void SetTimeT(int tube, int edge, Double_t t)
    {
        if (tube == 0)
            fTime1T_ns = t;
        else
            fTime2T_ns = t;
        //	cout<< "PaddleCalData trail "<<tube<<", "<<edge<<", "<<fTime1T_ns<<"; "<<fTime2T_ns<<endl;
    }
    /*
    void SetTime(int tube, int edge,Double_t t)
    {
        if (edge==1) fTime1L_ns=t;
        if (edge==2) fTime1T_ns=t;
        if (edge==3) fTime2L_ns=t;
        if (edge==4) fTime2T_ns=t;
    }*/

  private:
    UInt_t fPlane;
    UInt_t fBar;

  public:
    // ToDo: Would it be better to have arrays here?
    Double_t fTime1L_ns; // both PMs of a paddle
    Double_t fTime1T_ns; // both PMs of a paddle
    Double_t fTime2L_ns;
    Double_t fTime2T_ns;

    ClassDef(R3BPaddleCalData, 2)
};

#endif
