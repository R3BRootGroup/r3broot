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

#include <FairParGenericSet.h>

#define NCHMAX 5000

class FairParamList;

/**
 * Parameter container for time calibration of a single module. This class is used for
 * storage of time calibration parameters for a detector module. It contains
 * parametrisation of a table, used for TDC -> time [ns] conversion. Currently
 * supported systems: TACQUILA and VFTX.
 * @author D. Kresan
 * @since September 2, 2015
 */

class R3BTCalModulePar : public FairParGenericSet
{
  public:
    /**
     * Standard constructor.
     * @param name a name of container.
     * @param title a title of container.
     * @param context context/purpose for parameters and conditions.
     * @param own class ownership, if flag is kTRUE FairDB has the par.
     */
    R3BTCalModulePar(const char* name = "TCalModulePar",
                     const char* title = "TCAL calibration of a module",
                     const char* context = "TestDefaultContext",
                     Bool_t own = kTRUE);

    /**
     * Destructor.
     * Frees the memory allocated by the object.
     */
    virtual ~R3BTCalModulePar(void);

    /**
     * A method to reset the parameter values. Sets all parameters to 0.
     */
    void clear(void);

    /**
     * A method to write parameters using RuntimeDB.
     * @param list a list of parameters.
     */
    void putParams(FairParamList* list);

    /**
     * A method to read parameters using RuntimeDB.
     * @param list a list of parameters.
     * @return kTRUE if successfull, else kFALSE.
     */
    Bool_t getParams(FairParamList* list);

    /**
     * A method to print value of parameters to the standard
     * output using FairLogger.
     */
    void printParams();

    /**
     * A method to draw vaues of parameters on the current Canvas.
     */
    void DrawParams();

    /**
     * Member function for converting TDC value into time [ns]
     * using calibration parameters for clock TDC electronics.
     * @param tdc a TDC value.
     * @return time value in nanoseconds.
     */
    Double_t GetTimeClockTDC(Int_t tdc);

    /**
     * Member function for converting TDC value into time [ns]
     * using calibration parameters for TACQUILA electronics.
     * @param tdc a TDC value.
     * @return time value in nanoseconds.
     */
    Double_t GetTimeTacquila(Int_t tdc);

    /**
     * Member function for converting TDC value into time [ns]
     * using calibration parameters for VFTX electronics.
     * @param tdc a TDC value.
     * @return time value in nanoseconds.
     */
    Double_t GetTimeVFTX(Int_t tdc);

    /** Accessor functions **/
    Int_t GetPlane() const { return fPlane; }
    Int_t GetPaddle() const { return fPaddle; }
    Int_t GetSide() const { return fSide; }
    Int_t GetNofChannels() const { return fNofChannels; }
    Double_t GetSlopeAt(Int_t i) const { return fSlope[i]; }
    Double_t GetOffsetAt(Int_t i) const { return fOffset[i]; }
    Int_t GetBinLowAt(Int_t i) const { return fBinLow[i]; }
    Int_t GetBinUpAt(Int_t i) const { return fBinUp[i]; }
    void SetPlane(Int_t i) { fPlane = i; }
    void SetPaddle(Int_t i) { fPaddle = i; }
    void SetSide(Int_t i) { fSide = i; }
    void IncrementNofChannels() { fNofChannels += 1; }
    void SetBinLowAt(Int_t ch, Int_t i) { fBinLow[i] = ch; }
    void SetBinUpAt(Int_t ch, Int_t i) { fBinUp[i] = ch; }
    void SetSlopeAt(Double_t slope, Int_t i) { fSlope[i] = slope; }
    void SetOffsetAt(Double_t offset, Int_t i) { fOffset[i] = offset; }

  private:
    Int_t fPlane;             /**< Index of a plane. */
    Int_t fPaddle;            /**< Index of a paddle. */
    Int_t fSide;              /**< Side of a module: for NeuLAND - L/R PMT. */
    Int_t fNofChannels;       /**< Number of calibration parameters. */
    Int_t fBinLow[NCHMAX];    /**< Lower TDC range of a linear segment. */
    Int_t fBinUp[NCHMAX];     /**< Upper TDC range of a linear segment. */
    Double_t fSlope[NCHMAX];  /**< Slope of liear interpolation. */
    Double_t fOffset[NCHMAX]; /**< Offset of linear interpolation [ns]. */

  public:
    ClassDef(R3BTCalModulePar, 1);
};
