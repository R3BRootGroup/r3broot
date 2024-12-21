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

// ----------------------------------------------------------------------------
// MUSLI stands for MUSIC LIGHT ION
// ----------------------------------------------------------------------------
//   for experiments with light ions, the VME MDPP16 modules should be operated
//   with a different firmware. This firmware is working with the sum of the
//   signals sent to the input of the MDPP16 module. Therefore, the PA signals
//   of the anodes are sent to the input of the MDPP16 module.
// ----------------------------------------------------------------------------
//   The outputs are the following:
//     - output 0 = E and T from the average of input 0 and input 1 signals
//     - output 1 = E and T from the average of input 2 and input 3 signals
//     - output 2 = E and T from the average of input 4 and input 5 signals
//     - output 3 = E and T from the average of input 6 and input 7 signals
//     - output 4 = E and T from the average of input 8 and input 9 signals
//     - output 5 = E and T from the average of input 10 and input 11 signals
//     - output 6 = E and T from the average of input 12 and input 13 signals
//     - output 7 = E and T from the average of input 14 and input 15 signals
//     - output 8 = E and T from the average of input 0 to 3 signals
//     - output 9 = E and T from the average of input 4 to 7 signals
//     - output 10 = E and T from the average of input 8 to 11 signals
//     - output 11 = E and T from the average of input 12 to 15 signals
//     - output 12 = E and T from the average of input 0 to 7 signals
//     - output 13 = E and T from the average of input 8 to 15 signals
//     - outputs 14 and 15 = E and T from the average of input 0 to 15 (SAME SIGNAL)
//     - Tref (used later on for the Tstart of the Drift Time calculation)
//     - Ttrig (time of the master start sent to trigger the module: no use)
// ----------------------------------------------------------------------------
//   NB: T from output 0 to 15 are used as Tstop of the Drift Time calculation
//       for the theta calculation, you should use the T from output 0 to 7
// ----------------------------------------------------------------------------
//   The Mapped Data has the following structure
//     UInt_t fType;
//            fType=t, with t = 2, 4, 8, 16 mean signal on x anodes (output 0 to 13)
//                          t = 17           Tref or Trig
//     UInt_t fIndex;
//            fIndex = i, with i = 1, 2, 3, 4, 5, 6, 7, 8 for fType = 2
//                        with i = 1, 2, 3, 4             for fType = 4
//                        with i = 1, 2                   for fType = 8
//                        with i = 1                      for fType = 16
//                        with i = 1                      for fType = 17, Tref
//                        with i = 2                      for fType = 17, Ttrig
//
//     Uint_t fSignal = [1..8], [9..12], [13..4], [15], [17..18]
//     UInt_t fTime;
//     UInt_t fEnergy;
//     Bool_t fPileup;
//     Bool_t fOverflow;
// ----------------------------------------------------------------------------
//  NB: we keep 17 and 18 for Tref and Ttrig, so we are coherent
//  with R3BMusic, Twim, Trim
// ----------------------------------------------------------------------------

#ifndef R3BMusliMappedData_H
#define R3BMusliMappedData_H 1

#include "TObject.h"

class R3BMusliMappedData : public TObject
{
  public:
    /** Default constructor **/
    R3BMusliMappedData();

    R3BMusliMappedData(UInt_t type, UInt_t index, UInt_t signal, UInt_t time, UInt_t energy, Bool_t pu, Bool_t ov);

    /** Destructor **/
    virtual ~R3BMusliMappedData() {}

    /** Accessors **/
    inline const UInt_t& GetType() const { return fType; }
    inline const UInt_t& GetIndex() const { return fIndex; }
    inline const UInt_t& GetSignal() const { return fSignal; }
    inline const UInt_t& GetTime() const { return fTime; }
    inline const UInt_t& GetEnergy() const { return fEnergy; }
    inline const Bool_t& GetPileupStatus() const { return fPileup; }
    inline const Bool_t& GetOverflowStatus() const { return fOverflow; }

    /** Modifiers **/
    void SetType(UInt_t type) { fType = type; };
    void SetIndex(UInt_t index) { fIndex = index; };
    void SetSignal(UInt_t signal) { fSignal = signal; };
    void SetTime(UInt_t time) { fTime = time; };
    void SetEnergy(UInt_t energy) { fEnergy = energy; };
    void SetPileup(Bool_t pu) { fPileup = pu; }
    void SetOverflow(Bool_t ov) { fOverflow = ov; }

  protected:
    UInt_t fType;
    UInt_t fIndex;
    UInt_t fSignal;
    UInt_t fTime;
    UInt_t fEnergy;
    Bool_t fPileup;
    Bool_t fOverflow;

  public:
    ClassDef(R3BMusliMappedData, 1)
};

#endif
