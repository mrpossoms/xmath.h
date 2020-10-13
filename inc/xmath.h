// Library extending the standard C math library to provide more advanced
// mathematical operations.
// Copyright (C) 2020 - Kirk Roerig

// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.

// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.

// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <https://www.gnu.org/licenses/>.
#ifndef _XMATH_H_
#define _XMATH_H_
#include <math.h>

#ifndef REAL
#define REAL float
#endif

static inline void vec_add(size_t n, REAL dst[n], REAL left[n], REAL right[n])
{
	for (size_t i = 0; i < n; i++)
	{
		dst[i] = left[i] + right[i];
	}
}


static inline void vec_sub(size_t n, REAL dst[n], REAL left[n], REAL right[n])
{
	for (size_t i = 0; i < n; i++)
	{
		dst[i] = left[i] - right[i];
	}
}


static inline void vec_scl(size_t n, REAL dst[n], REAL in[n], REAL scalar)
{
	for (size_t i = 0; i < n; i++)
	{
		dst[i] = in[i] * scalar;
	 }
}


static inline REAL vec_dot(size_t n, REAL left[n], REAL right[n])
{
	REAL dot = 0;
	for (size_t i = 0; i < n; i++)
	{
		dot += left[i] * right[i];
	}

	return dot;
}

static inline void vec_hadamard(size_t n, REAL dst[n], REAL left[n], REAL right[n])
{
	for (size_t i = 0; i < n; i++)
	{
		dst[i] = left[i] * right[i];
	}
}


static inline void vec_each_ele(size_t n, REAL dst[n], REAL src[n], REAL (*func)(REAL x))
{
	for (size_t i = 0; i < n; i++)
	{
		dst[i] = func(src[i]);
	}
}


static inline void mat_add()

#endif