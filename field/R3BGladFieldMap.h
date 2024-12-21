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

#include "R3BFieldPar.h"

#include <FairField.h>
#include <TRotation.h>
#include <TString.h>
#include <TTree.h>
#include <TVector3.h>

class TArrayD;

class R3BGladFieldMap : public FairField
{
  public:
    /** Default constructor **/
    R3BGladFieldMap();

    /** Standard constructor
     ** @param name       Name of field map file without extension
     ** @param fileType   R = ROOT file, A = ASCII
     **/
    R3BGladFieldMap(const TString& mapName, const TString& fileType = "A");

    /** Constructor from R3BGladFieldPar **/
    R3BGladFieldMap(R3BFieldPar* fieldPar);

    /** Constructor from R3BGladFieldMapCreator **/
    // R3BGladFieldMap(R3BGladFieldMapCreator* creator);

    /** Destructor **/
    virtual ~R3BGladFieldMap();

    /** Initialisation (read map from file) **/
    virtual void Init();

    /** Get the field components at a certain point
     ** @param x,y,z     Point coordinates (global) [cm]
     ** @value Bx,By,Bz  Field components [kG]
     **/
    virtual Double_t GetBx(Double_t x, Double_t y, Double_t z);
    virtual Double_t GetBy(Double_t x, Double_t y, Double_t z);
    virtual Double_t GetBz(Double_t x, Double_t y, Double_t z);

    // special function for the field transformation
    TVector3 GetBtrans(Double_t x, Double_t y, Double_t z);

    /** Determine whether a point is inside the field map
     ** @param x,y,z              Point coordinates (global) [cm]
     ** @param ix,iy,iz (return)  Grid cell
     ** @param dx,dy,dz (return)  Distance from grid point [cm] if inside
     ** @value kTRUE if inside map, else kFALSE
     **/
    virtual Bool_t IsInside(Double_t x,
                            Double_t y,
                            Double_t z,
                            Int_t& ix,
                            Int_t& iy,
                            Int_t& iz,
                            Double_t& dx,
                            Double_t& dy,
                            Double_t& dz);

    /** Write the field map to an ASCII file **/
    void WriteAsciiFile(const TString& fileName);

    /** Set the position (in cm) of the field origin in lab **/
    virtual void SetPosition(Double_t x, Double_t y, Double_t z);

    /** Set a global field scaling factor **/
    virtual void SetScale(Double_t factor) { fScale = factor; }
    virtual void SetFieldfromCurrent(double current);

    /* Set Euler rotation angles of the field (in degrees)
     * default fYAngle = -14 deg, fXAngle=0, fZAngle=0 */
    virtual void SetXAngle(Double_t a) { fXAngle = a; };
    virtual void SetYAngle(Double_t a) { fYAngle = a; };
    virtual void SetZAngle(Double_t a) { fZAngle = a; };

    /** Accessors to field parameters in local coordinate system **/
    Double_t GetXmin() const { return fXmin; }
    Double_t GetYmin() const { return fYmin; }
    Double_t GetZmin() const { return fZmin; }
    Double_t GetXmax() const { return fXmax; }
    Double_t GetYmax() const { return fYmax; }
    Double_t GetZmax() const { return fZmax; }
    Double_t GetXstep() const { return fXstep; }
    Double_t GetYstep() const { return fYstep; }
    Double_t GetZstep() const { return fZstep; }
    Int_t GetNx() const { return fNx; }
    Int_t GetNy() const { return fNy; }
    Int_t GetNz() const { return fNz; }

    /** Accessor to field centre position in lab system **/
    Double_t GetPositionX() const { return fPosX; }
    Double_t GetPositionY() const { return fPosY; }
    Double_t GetPositionZ() const { return fPosZ; }

    /** Accessor to field rotation **/
    Double_t GetXAngle() const { return fXAngle; }
    Double_t GetYAngle() const { return fYAngle; }
    Double_t GetZAngle() const { return fZAngle; }

    /** Accessor to global scaling factor  **/
    Double_t GetScale() const { return fScale; }

    /** Accessors to the field value arrays **/
    TArrayD* GetBx() const { return fBx; }
    TArrayD* GetBy() const { return fBy; }
    TArrayD* GetBz() const { return fBz; }

    /** Accessor to field map file **/
    TString GetFileName() { return fFileName; }

    /** Screen output **/
    virtual void Print(Option_t* option = "") const;

  protected:
    /** Reset the field parameters and data **/
    void Reset();

    /** Read the field map from an ASCII file **/
    void ReadAsciiFile(const TString& fileName);

    /** Read field map from a ROOT file **/
    void ReadRootFile(const TString& fileName);

    /** Get field values by interpolation of the grid.
     ** @param dx,dy,dz  Relative distance from grid point [cell units]
     **/
    Double_t Interpolate(Double_t dx, Double_t dy, Double_t dz);

    /** Map file name **/
    TString fFileName;

    /** Global scaling factor (w.r.t. map on file) **/
    Double_t fScale;

    /** Field centre position in global coordinates  **/
    Double_t fPosX, fPosY, fPosZ;

    /** Euler rotations of the field around field local X, Y, Z axis **/
    Double_t fXAngle;
    Double_t fYAngle;
    Double_t fZAngle;

    /** Field limits in local coordinate system **/
    Double_t fXmin, fXmax, fXstep;
    Double_t fYmin, fYmax, fYstep;
    Double_t fZmin, fZmax, fZstep;

    /** Number of grid points  **/
    Int_t fNx, fNy, fNz; //

    /** Arrays of doubles with the field values  **/
    TArrayD* fBx; //!
    TArrayD* fBy; //!
    TArrayD* fBz; //!

    /** Variables for temporary storage
     ** Used in the very frequently called method GetFieldValue  **/
    Double_t fHa[2][2][2]; //! Field at corners of a grid cell
    Double_t fHb[2][2];    //! Interpolated field (2-dim)
    Double_t fHc[2];       //! Interpolated field (1-dim)

    // local transformation
    TRotation* gRot;  //!
    TVector3* gTrans; //!

    // TTree with the map data when reading a ROOT file
    TFile* fFile; // root file with the map data

  public:
    ClassDef(R3BGladFieldMap, 5)
};
