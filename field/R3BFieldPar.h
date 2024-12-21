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

// -------------------------------------------------------------------------
// -----                      R3BFieldPar header file                  -----
// -------------------------------------------------------------------------

/** R3BFieldPar.h
 ** @author
 ** @since 20.02.2006
 ** @version 1.0
 **
 ** Parameter set for the CBM magnetic field. For the runtime database.
 **/

#pragma once

#include <FairField.h>
#include <FairParGenericSet.h>

class FairParamList;

const int kMaxFieldMapType = 5;

class R3BFieldPar : public FairParGenericSet
{
  public:
    /** Standard constructor  **/
    R3BFieldPar(const char* name, const char* title, const char* context);

    /** default constructor  **/
    R3BFieldPar();

    /** Destructor **/
    ~R3BFieldPar();

    /** Put parameters **/
    virtual void putParams(FairParamList* list);

    /** Get parameters **/
    virtual Bool_t getParams(FairParamList* list);

    /** Set parameters from R3BField  **/
    void SetParameters(FairField* field);

    /** Accessors **/
    Int_t GetType() const { return fType; }
    Double_t GetXmin() const { return fXmin; }
    Double_t GetXmax() const { return fXmax; }
    Double_t GetYmin() const { return fYmin; }
    Double_t GetYmax() const { return fYmax; }
    Double_t GetZmin() const { return fZmin; }
    Double_t GetZmax() const { return fZmax; }
    Double_t GetBx() const { return fBx; }
    Double_t GetBy() const { return fBy; }
    Double_t GetBz() const { return fBz; }
    void MapName(TString& name) { name = fMapName; }
    TString GetFileName() { return fMapFileName; }
    Double_t GetPositionX() const { return fPosX; }
    Double_t GetPositionY() const { return fPosY; }
    Double_t GetPositionZ() const { return fPosZ; }
    Double_t GetXAngle() const { return fXAngle; }
    Double_t GetYAngle() const { return fYAngle; }
    Double_t GetZAngle() const { return fZAngle; }
    Double_t GetScale() const { return fScale; }
    Double_t GetCurrent() const { return fCurrent; }

  private:
    /** Field type
     ** 0 = constant field
     ** 1 = field map
     ** 2 = Glad field map
     **/
    Int_t fType;

    /** Field limits in case of constant field **/
    Double_t fXmin, fXmax;
    Double_t fYmin, fYmax;
    Double_t fZmin, fZmax;

    /** Field values in case of constant field [kG] **/
    Double_t fBx, fBy, fBz;

    /** Field map name in case of field map **/
    TString fMapName;
    TString fMapFileName;

    /** Field centre position for field map **/
    Double_t fPosX, fPosY, fPosZ;

    /** Euler rotations of the field map **/
    Double_t fXAngle, fYAngle, fZAngle;

    /** Scaling factor for field map **/
    Double_t fScale;

    Double_t fCurrent;

  public:
    ClassDef(R3BFieldPar, 2);
};
