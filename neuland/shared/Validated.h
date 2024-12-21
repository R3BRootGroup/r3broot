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

#ifndef R3BVALIDATED_H
#define R3BVALIDATED_H

//=============================================================================
//
// D A templatised class to represent data which can be validated, for instance
// D for data caching.
// D Original: https://bitbucket.org/davidcorne/mutable/
// D Modified by Jan Mayer

#include <assert.h>
#include <utility>

//=============================================================================
template <typename T>
class Validated
{
  public:
    Validated();
    // invalid constructor

    Validated(const T& data);
    Validated(T&& data);
    // valid constructors

    bool valid() const;
    // is this valid or not?

    void invalidate();
    // makes the class invalid again, so call in a function which will change the
    // cached value.

    void set(const T& data);
    void set(T&& data);
    // set the data, makes the instance valid.

    T get() const;
    // explicit conversion

    const T& getRef();
    // getter with const ref. Use this if you only want to check the value.

    operator T();
    operator T() const;
    // implicit conversions to the template type.

  private:
    bool m_valid;
    T m_data;
};

//=============================================================================
//
// D Include the source in the header file as it is a templatised class.
//

//=============================================================================
template <typename T>
Validated<T>::Validated()
    //
    // D Invalid constructor
    //
    : m_valid(false)
{
}

//=============================================================================
template <typename T>
Validated<T>::Validated(const T& data)
    //
    // D Valid constructor
    //
    : m_valid(true)
    , m_data(data)
{
}

//=============================================================================
template <typename T>
Validated<T>::Validated(T&& data)
    //
    // D Valid constructor
    //
    : m_valid(true)
    , m_data(std::move(data))
{
}

//=============================================================================
template <typename T>
bool Validated<T>::valid() const
//
// D is this valid or not?
//
{
    return m_valid;
}

//=============================================================================
template <typename T>
void Validated<T>::invalidate()
//
// D Makes the class invalid again, so call when the client calls a function
// D which will change the cached value.
//
{
    m_valid = false;
}

//=============================================================================
template <typename T>
void Validated<T>::set(const T& data)
//
// D set the data.
//
{
    m_data = data;
    m_valid = true;
}

//=============================================================================
template <typename T>
void Validated<T>::set(T&& data)
//
// D set the data.
//
{
    m_data = std::move(data);
    m_valid = true;
}

//=============================================================================
template <typename T>
T Validated<T>::get() const
{
    assert(valid() && "Cannot return invalid data.");
    return m_data;
}

//=============================================================================
template <typename T>
const T& Validated<T>::getRef()
{
    assert(valid() && "Cannot return invalid data.");
    return m_data;
}

//=============================================================================
template <typename T>
Validated<T>::operator T()
//
// D implicit conversion to the template type.
//
{
    assert(valid() && "Cannot return invalid data.");
    return m_data;
}

//=============================================================================
template <typename T>
Validated<T>::operator T() const
//
// D implicit conversion to the template type.
//
{
    // a nice way of returning a message with a standard assert.
    assert(valid() && "Cannot return invalid data.");
    return m_data;
}

#endif // R3BVALIDATED_H
