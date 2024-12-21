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

#ifndef R3BFIELDINTERP_H
#define R3BFIELDINTERP_H

#include <stdlib.h>

// Class to interpolate one variable in three dimensions Using linear
// interpolation.  When given a point outside the valid map: produce
// values as at the boundary at that point, i.e.  give a continous
// value outside.  But not where and in what direction it went wrong.

class R3BFieldInterp
{
  public:
    R3BFieldInterp()
        : _m1(0)
        , _m2(0)
        , _n(0)
        , _data(NULL)
    {
        for (int i = 0; i < 3; i++)
            _np[i] = 0;
    }

    ~R3BFieldInterp() { free(_data); }

  private:
    const R3BFieldInterp& operator=(const R3BFieldInterp&);
    R3BFieldInterp(const R3BFieldInterp&);

  public:
    void interpolate(const R3BFieldInterp& s1, double w1, const R3BFieldInterp& s2, double w2);

    void init();

    bool expand();

  public:
    double interp(int ic[3], double dc[3] /*,int &outside*/);

    double interp3(int ic[3], double dc[3] /*,int &outside*/);

  public:
    int _np[3];
    int _max_ic[3]; // _max_ic[i] = _np[i] - 1
    int _m1, _m2;   // _m1 = _np[1] * _np[2] ; _m2 = _np[2]
    int _n;         // _n = _np[0] * _np[1] * _np[2]
    float* _data;

    float get_data_pt(int i0, int i1, int i2) { return _data[i0 * _m1 + i1 * _m2 + i2]; }

    void set_data_pt(int i0, int i1, int i2, float d) { _data[i0 * _m1 + i1 * _m2 + i2] = d; }
};

#endif // R3BFIELD_INTERP_H_
