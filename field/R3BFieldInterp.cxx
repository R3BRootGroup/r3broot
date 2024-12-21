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

#include "R3BFieldInterp.h"
#include <iostream>
#include <math.h>

using namespace std;

void R3BFieldInterp::interpolate(const R3BFieldInterp& s1, double w1, const R3BFieldInterp& s2, double w2)
{
    for (int i = 0; i < 3; i++)
    {
        if (s1._np[i] != s2._np[i])
            cout << "-E- Field interpolation error, s1.dim != s2.dim" << endl;
        ;
        _np[i] = s1._np[i];
    }

    init();

    for (int i = 0; i < _n; i++)
    {
        _data[i] = (float)(w1 * s1._data[i] + w2 * s2._data[i]);
    }
}

void R3BFieldInterp::init()
{

    for (int i = 0; i < 3; i++)
    {
        _max_ic[i] = _np[i] - 1;
    }

    _m1 = _np[1] * _np[2];
    _m2 = _np[2];
    _n = _np[0] * _np[1] * _np[2];

    float* d = reinterpret_cast<float*>(realloc(_data, sizeof(float) * _n));

    if (!d)
        cout << "-E- Field interpolation, memory allocation failure." << endl;

    _data = d;

    for (int i = 0; i < _n; i++)
        _data[i] = NAN;
}

bool R3BFieldInterp::expand()
{
    bool expanded = false;

    for (int i0 = 0; i0 < _np[0]; i0++)
        for (int i1 = 0; i1 < _np[1]; i1++)
            for (int i2 = 0; i2 < _np[2]; i2++)
                if (isnan(get_data_pt(i0, i1, i2)))
                {
                    double sum = 0;
                    double sum_w = 0;

                    for (int di0 = -1; di0 <= 1; di0++)
                        for (int di1 = -1; di1 <= 1; di1++)
                            for (int di2 = -1; di2 <= 1; di2++)
                            {
                                int i0_1 = i0 + di0;
                                int i1_1 = i1 + di1;
                                int i2_1 = i2 + di2;

                                if (i0_1 >= 0 && i0_1 < _np[0] && i1_1 >= 0 && i1_1 < _np[1] && i2_1 >= 0 &&
                                    i2_1 < _np[2])
                                {
                                    float v1 = get_data_pt(i0_1, i1_1, i2_1);

                                    if (!isnan(v1))
                                    {
                                        double w = 1 / sqrt(di0 * di0 + di1 * di1 + di2 * di2);
                                        sum += w * v1;
                                        sum_w += w;
                                    }
                                }
                            }
                    if (sum_w > 0)
                    {
                        set_data_pt(i0, i1, i2, (float)(sum / sum_w));
                        expanded = true;
                    }
                }

    return expanded;
}

double R3BFieldInterp::interp(int ic[3], double dc[3] /*,int &outside*/)
{
    int ic0[3], ic1[3];
    /*
    printf ("[%d %d %d / %d %d %d : %7.5f %7.5f %7.5f]",
        _np[0],_np[1],_np[2],ic[0],ic[1],ic[2],dc[0],dc[1],dc[2]);
    */
    // outside = 0;

    // int outside_mark = 1;

    for (int i = 0; i < 3; i++)
    {
        ic0[i] = ic[i];
        ic1[i] = ic[i] + 1;

// #define RETNAN return NAN
#define RETNAN

        if (ic0[i] < 0)
        {
            ic0[i] = 0;
            RETNAN; /*outside |= outside_mark; */
        }
        if (ic1[i] < 0)
        {
            ic1[i] = 0;
            RETNAN;
        }
        // outside_mark <<= 1;
        if (ic0[i] > _max_ic[i])
        {
            ic0[i] = _max_ic[i];
            RETNAN; /*outside |= outside_mark; */
        }
        if (ic1[i] > _max_ic[i])
        {
            ic1[i] = _max_ic[i];
            RETNAN;
        }
        // outside_mark <<= 1;
    }
    /*
    printf ("[%d %d %d - %d %d %d]",
        ic0[0],ic0[1],ic0[2],ic1[0],ic1[1],ic1[2]);
    */
    float v000 = get_data_pt(ic0[0], ic0[1], ic0[2]);
    float v001 = get_data_pt(ic0[0], ic0[1], ic1[2]);
    float v010 = get_data_pt(ic0[0], ic1[1], ic0[2]);
    float v011 = get_data_pt(ic0[0], ic1[1], ic1[2]);
    float v100 = get_data_pt(ic1[0], ic0[1], ic0[2]);
    float v101 = get_data_pt(ic1[0], ic0[1], ic1[2]);
    float v110 = get_data_pt(ic1[0], ic1[1], ic0[2]);
    float v111 = get_data_pt(ic1[0], ic1[1], ic1[2]);

    // Either we always operate with the raw data, which makes it rather
    // straight-forward to handle the out-of-bounds cases (by simply
    // giving the same values for the two interpolations to do)

    double vx00 = v000 + dc[0] * (v100 - v000);
    double vx01 = v001 + dc[0] * (v101 - v001);
    double vx10 = v010 + dc[0] * (v110 - v010);
    double vx11 = v011 + dc[0] * (v111 - v011);

    double vxy0 = vx00 + dc[1] * (vx10 - vx00);
    double vxy1 = vx01 + dc[1] * (vx11 - vx01);

    double vxyz = vxy0 + dc[2] * (vxy1 - vxy0);

    return vxyz;

#if 0
  // Or we precalculate the differences to the next points.  Then,
  // whenever coming outside range, we set the respective differences
  // to zero.

  // The precalculated differences does require 8 times the storage,
  // as for each cube, 8 values must be remembered.  The values would
  // be adjacent in memory.

  /*
  v = v0 + 
    dvdx * dx + dvdy * dy + dvdz * dz + 
    dvdxdy * dx * dy + dvdxdz * dx * dz + dvdydz * dy * dz + 
    dvdxdydz * dx * dy * dz;
  */
  v = v0 + 
    dx * (dvdx + dy * (dvdxdy + dz * dvdxdydz) + dz * dvdxdz) + 
    dy * (dvdy + dz * dvdydz) + 
    dz * dvdz;
#endif
}

inline void interp3_factors(double v0, double v1, double v2, double v3, double f[4])
{
    f[0] = v1;
    f[1] = -.5 * v0 + .5 * v2;
    f[2] = v0 - 2.5 * v1 + 2 * v2 - .5 * v3;
    f[3] = -.5 * v0 + 1.5 * v1 - 1.5 * v2 + .5 * v3;
}

inline double interp3_factors_mul(double v[4], double p1, double p2, double p3)
{
    double s = v[1];
    s += p1 * (-.5 * v[0] + .5 * v[2]);
    s += p2 * (v[0] - 2.5 * v[1] + 2 * v[2] - .5 * v[3]);
    s += p3 * (-.5 * v[0] + 1.5 * v[1] - 1.5 * v[2] + .5 * v[3]);
    return s;
}

struct interp3_cell
{
    double _f[4][4][4];
};

double R3BFieldInterp::interp3(int ic[3], double dc[3] /*,int &outside*/)
{
    // Make field interpolation that also takes neighbouring cells into
    // account, for a smoother field map.  Within each cell, use a
    // weighted average of the linear inter(extra)polations of the
    // points in this and neighbouring cells.  The weigths are
    // quadratic, being at the boundary half, and going to zero at the
    // other boundary.  One-dimensionally:
    //
    // We are in the mid cell, i.e. between t1 and t2.  At t1 values is
    // v1, the previous cell has index 0 of the other border and the
    // next has index 3.  I.e. (t0,v0) (t1,v1) (t2,v2) (t3,v3)
    //
    // Interpolation:  va = ( -x)*v0 + (x+1)*v1  (left cell)
    //                 vb = (1-x)*v1 + (x  )*v2  (this cell)
    //                 vc = (2-x)*v2 + (x-1)*v3  (right cell)
    //
    // Weights:        wa = 1/2 * (1-x)^2        (left cell)
    //                 wb = 3/4 - (1/2-x)^2      (this cell)
    //                 wc = 1/2 * x^2            (right cell)
    //
    // These are valid for 0<=x<=1, and should be summed:
    //
    // v(x) = va*wa+vb*wb+vc*wc =
    //      =       (          v1                  ) +
    //        x   * (-v0/2 +            v2/2       ) +
    //        x^2 * ( v0   - 5*v1/2 + 2*v2   - v3/2) +
    //        x^3 * (-v0/2 + 3*v1/2 - 3*v2/2 + v3/2)
    //

    int icj[3][4];

    // outside = 0;

    // int outside_mark = 1;

    for (int i = 0; i < 3; i++)
    {
        icj[i][0] = ic[0] - 1;
        icj[i][1] = ic[0];
        icj[i][2] = ic[0] + 1;
        icj[i][3] = ic[0] + 2;

        if (icj[i][0] < 0)
        {
            icj[i][0] = 0;
        }
        if (icj[i][1] < 0)
        {
            icj[i][1] = 0; /*outside |= outside_mark; */
        }
        if (icj[i][2] < 0)
        {
            icj[i][2] = 0;
        }
        if (icj[i][3] < 0)
        {
            icj[i][3] = 0;
        }
        /*outside_mark <<= 1;*/
        if (icj[i][0] > _max_ic[i])
        {
            icj[i][0] = _max_ic[i];
        }
        if (icj[i][1] > _max_ic[i])
        {
            icj[i][1] = _max_ic[i]; /*outside |= outside_mark; */
        }
        if (icj[i][2] > _max_ic[i])
        {
            icj[i][2] = _max_ic[i];
        }
        if (icj[i][3] > _max_ic[i])
        {
            icj[i][3] = _max_ic[i];
        }
        /*outside_mark <<= 1;*/
    }

    // First we get the values for the box.  That is to get 4x4x4 = 64 values.
    // We at the same time calculate the factors of them, i.e. rebase such that the next
    // step is to just apply the first multiplicative

    interp3_cell ip3c;

    for (int i0 = 0; i0 < 4; i0++)
        for (int i1 = 0; i1 < 4; i1++)
        {
            int offset = icj[0][i0] * _m1 + icj[1][i1] * _m2;

            double v[4];

            for (int i2 = 0; i2 < 4; i2++)
                v[i2] = _data[offset + icj[2][i2]];

            interp3_factors(v[0], v[1], v[2], v[3], ip3c._f[i0][i1]);
        }

    // First we cook down the 4x4x4 values to 4x4 by multiplying with the first set
    // of powers (of 3rd component)

    double ip2s[4][4];

    {
        double p1 = dc[2];
        double p2 = dc[2] * dc[2];
        double p3 = dc[2] * p2;

        for (int i0 = 0; i0 < 4; i0++)
            for (int i1 = 0; i1 < 4; i1++)
                ip2s[i0][i1] = interp3_factors_mul(ip3c._f[i0][i1], p1, p2, p3);
    }

    double ip1l[4];

    {
        double p1 = dc[1];
        double p2 = dc[1] * dc[1];
        double p3 = dc[1] * p2;

        for (int i0 = 0; i0 < 4; i0++)
            ip1l[i0] = interp3_factors_mul(ip2s[i0], p1, p2, p3);
    }

    double ip0;

    {
        double p1 = dc[0];
        double p2 = dc[0] * dc[0];
        double p3 = dc[0] * p2;

        ip0 = interp3_factors_mul(ip1l, p1, p2, p3);
    }

    return ip0;
}

/*

In[4]:= Wa=CC*(1-x) ^2

             2
      (1 - x)
Out[4]= --------
         2

In[5]:= Wc=CC*x^2

       2
      x
Out[5]= --
      2

In[6]:= Wb=3/4-(1/2-x)^2

      3    1     2
Out[6]= - - (- - x)
      4    2

In[7]:= Wa+Wb+Wc

                            2    2
      3    1     2   (1 - x)    x
Out[7]= - - (- - x)  + -------- + --
      4    2            2       2

In[8]:= Simplify[Wa+Wb+Wc]

Out[8]= 1

In[9]:= ya=(1+x)*y2-x*y1

Out[9]= -(x y1) + (1 + x) y2

In[10]:= yb=x*y3+(1-x)*y2

Out[10]= (1 - x) y2 + x y3

In[11]:= yc=(x-1)*y4+(2-x)*y3

Out[11]= (2 - x) y3 + (-1 + x) y4

In[12]:= ya/.x->0

Out[12]= y2

In[13]:= ya/.x->1

Out[13]= -y1 + 2 y2

In[14]:= ya/.x->-1

Out[14]= y1

In[15]:= yb/.x->0

Out[15]= y2

In[16]:= yb/.x->1

Out[16]= y3

In[17]:= yc/.x->1

Out[17]= y3

In[25]:= yc/.x->2

Out[25]= y4

In[24]:= Collect[Expand[ya*Wa+yb*Wb+yc*Wc],x]

               -y1   y3     2       5 y2          y4
Out[24]= y2 + x (--- + --) + x  (y1 - ---- + 2 y3 - --) +
                2    2               2            2

    3  -y1   3 y2   3 y3   y4
>    x  (--- + ---- - ---- + --)
        2     2      2     2

In[29]:= CForm[%24]

Out[29]//CForm=
 y2 + x*(-y1/2. + y3/2.) + Power(x,2)*(y1 - (5*y2)/2. + 2*y3 - y4/2.) +
  Power(x,3)*(-y1/2. + (3*y2)/2. - (3*y3)/2. + y4/2.)


*/
