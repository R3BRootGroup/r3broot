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

#include "R3BFieldInterp.h"
#include "R3BFieldPar.h"
#include <FairField.h>
#include <TRotation.h>
#include <TVector3.h>
#include <map>

// ------------------  Land02 internal Structures -----------------------//
// <DB @ March 2010 >
// The syntax convention in Land02 is preserved, member variables
// begins with an _<xx>. It will help to separated between imported Land02
// code and native one.

struct fields_ALADiN
{
    R3BFieldInterp f[2][3];
};

typedef std::map<Double_t, fields_ALADiN*> map_fields_ALADiN;

struct coords_ALADiN
{
    // Note, these are in ALADiN measurement box axis system,
    // i.e. x is forward (usual z) and z is usual -x

    TVector3 fBox_pt[2];
    TVector3 fMag_pt[2];
    Double_t fCosa, fSina;
};

// -------------------- Standard R3BRoot Field Map API -----------------//

class TArrayF;

class R3BAladinFieldMap : public FairField
{

  public:
    /** Default constructor **/
    R3BAladinFieldMap();

    /** Standard constructor
     ** @param name       Name of field map
     ** @param fileType   R = ROOT file, A = ASCII
     **/
    R3BAladinFieldMap(const char* mapName, const char* fileType = "A");

    /** Constructor from R3BAladinFieldPar **/
    R3BAladinFieldMap(R3BFieldPar* fieldPar);

    /** Constructor from R3BAladinFieldMapCreator **/
    // R3BAladinFieldMap(R3BAladinFieldMapCreator* creator);

    /** Destructor **/
    virtual ~R3BAladinFieldMap();

    /** Initialisation (read map from file) **/
    virtual void Init();

    virtual void InitField();

    virtual void FillParContainer() {}

    /** Main GetField function */
    virtual void GetFieldValue(const Double_t point[3], Double_t* bField);

    /** Get the field components at a certain point
     ** @param x,y,z     Point coordinates (global) [cm]
     ** @value Bx,By,Bz  Field components [kG]
     **/
    virtual Double_t GetBx(Double_t x, Double_t y, Double_t z);
    virtual Double_t GetBy(Double_t x, Double_t y, Double_t z);
    virtual Double_t GetBz(Double_t x, Double_t y, Double_t z);

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
    void WriteAsciiFile(const char* fileName);

    /** Write field map data to a ROOT file **/
    // void WriteRootFile(const char* fileName, const char* mapName);

    /** Set the position of the field centre **/
    virtual void SetPosition(Double_t x, Double_t y, Double_t z);

    /** Set a global field scaling factor **/
    virtual void SetScale(Double_t factor) { fScale = factor; }

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

    /** Accessor to field centre position in global system **/
    Double_t GetPositionX() const { return fPosX; }
    Double_t GetPositionY() const { return fPosY; }
    Double_t GetPositionZ() const { return fPosZ; }

    /** Accessor to global scaling factor  **/
    Double_t GetScale() const { return fScale; }

    /** Accessors to the field value arrays **/
    TArrayF* GetBx() const { return fBx; }
    TArrayF* GetBy() const { return fBy; }
    TArrayF* GetBz() const { return fBz; }

    /** Accessor to field map file **/
    const char* GetFileName() { return fFileName.Data(); }

    void CalcFieldDiv(R3BFieldInterp f[3], Double_t d[3]);

    /** Screen output **/
    virtual void Print(Option_t* option = "") const;

    /** Set Current **/

    void SetCurrent(Double_t aCurrent) { fCurrent = aCurrent; }
    Double_t GetCurrent() { return fCurrent; }
    void SetFringeField(Bool_t set) { gFringeField = set; }

  protected:
    /** Reset the field parameters and data **/
    void Reset();

    /** Read the field map from an ASCII file **/
    void ReadAsciiFile(const char* fileName);

    /** Read field map from a ROOT file **/
    // void ReadRootFile(const char* fileName, const char* mapName);

    /** Set field parameters and data **/
    // void SetField(const R3BAladinFieldMapData* data);

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

    /** Field limits in local coordinate system **/
    Double_t fXmin, fXmax, fXstep;
    Double_t fYmin, fYmax, fYstep;
    Double_t fZmin, fZmax, fZstep;

    /** Number of grid points  **/
    Int_t fNx, fNy, fNz; //

    /** Arrays with the field values  **/
    TArrayF* fBx; //!
    TArrayF* fBy; //!
    TArrayF* fBz; //!

    /** Variables for temporary storage
     ** Used in the very frequently called method GetFieldValue  **/
    Double_t fHa[2][2][2]; //! Field at corners of a grid cell
    Double_t fHb[2][2];    //! Interpolated field (2-dim)
    Double_t fHc[2];       //! Interpolated field (1-dim)

    /** local transformation
     **/
    TRotation* gRot;  //!
    TVector3* gTrans; //!

    /** land02 imported variables
     **/

    static map_fields_ALADiN gMapIFieldOrig; //!
    static map_fields_ALADiN gMapIField;     //!

    static coords_ALADiN gCoords[2]; //!
    static Bool_t gInitialized;      //!

    fields_ALADiN* fCurField; //!
    Double_t fCurrent;        //!
    Double_t fFieldSign;      //!

    TVector3 af_box[2][2]; //!
    TVector3 af_mag[2][2]; //!
    Bool_t gFringeField;   //!

    ClassDef(R3BAladinFieldMap, 1)
};
