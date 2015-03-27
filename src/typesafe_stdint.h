/***************************************************************************
 *   Copyright (C) 2015 by Tim Molteno                                *
 *   tim@molteno.net                                                       *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 *   This program is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU General Public License for more details.                          *
 *                                                                         *
 *   You should have received a copy of the GNU General Public License     *
 *   along with this program; if not, write to the                         *
 *   Free Software Foundation, Inc.,                                       *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
 ***************************************************************************/

/*! \brief A Class to substitute for the stdint types in order to ensure
 * pedantic type safety. This will mean that adding two different integer
 * types will result in a compile-error!
 * 
 * \note Note that this class will not link. It does not define the operators that
 * it declares. So this should be used purely as a lint-like checking tool
 * to find problems with your code. For example, when changing an array
 * class to allow 64-bit lengths, all you have to do is declare the length
 * type as int64_t, and then define TYPESAFE_PEDANTIC and include this
 * class. All hidden type-promotions and coercions will be revealed.
 * 
 */
#if TYPESAFE_PEDANTIC 

template<typename T>
class typesafe_int {
public:
    typesafe_int(int x);
    
    bool operator==(const typesafe_int<T>);
    bool operator>=(const typesafe_int<T>);
    bool operator<(const typesafe_int<T>);
    bool operator>(const typesafe_int<T>);
    typesafe_int<T>& operator*(const typesafe_int<T>);
    typesafe_int<T>& operator/(const typesafe_int<T>);
    typesafe_int<T>& operator+(const typesafe_int<T>);
    typesafe_int<T>& operator-(const typesafe_int<T>);
    typesafe_int<T>& operator++(int);
};

typedef typesafe_int<long long> int64_t;
typedef typesafe_int<unsigned long long> uint64_t;

typedef typesafe_int<long> int32_t;
typedef typesafe_int<unsigned long> uint32_t;

#else
  #include <stdint.h>
#endif

