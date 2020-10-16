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

#ifndef TYPE
#define TYPE float
#endif

/**
 * @brief      Compute vector addition between two vectors
 *
 * @param      n          Dimensionality of the vectors
 * @param      dst_arr    The destination array to contain the result of size n
 * @param      left_arr   The left array operand of the addition of size n
 * @param      right_arr  The right array operand of the addition of size n
 *
 */
#define VEC_ADD(n, dst_arr, left_arr, right_arr) \
{\
	for (size_t i = 0; i < (n); i++)\
	{\
		(dst_arr)[i] = (left_arr)[i] + (right_arr)[i];\
	}\
}\

/**
 * @brief      Compute vector subtraction between two vectors
 *
 * @param      n          Dimensionality of the vectors
 * @param      dst_arr    The destination array to contain the result of size n
 * @param      left_arr   The left array operand of the subtraction of size n
 * @param      right_arr  The right array operand of the subtraction of size n
 *
 */
#define VEC_SUB(n, dst_arr, left_arr, right_arr) \
{\
	for (size_t i = 0; i < (n); i++)\
	{\
		(dst_arr)[i] = (left_arr)[i] - (right_arr)[i];\
	}\
}\

/**
 * @brief      Scales each element of a vector by a scalar
 *
 * @param      n          Dimensionality of the vectors
 * @param      dst_arr    The destination array to contain the result of size n
 * @param      left_arr   The left array operand of the addition of size n
 * @param      right_scalar  The scalar by which each element is multiplied.
 *
 */
#define VEC_SCL(n, dst_arr, left_arr, right_scalar) \
{\
	for (size_t i = 0; i < (n); i++)\
	{\
		(dst_arr)[i] = (left_arr)[i] * (right_scalar);\
	}\
}\

/**
 * @brief      Computes the Hadamard product (element wise product) between two
 *             vectors.
 *
 * @param      n          Dimensionality of the vectors.
 * @param      dst_arr    The destination array to contain the result of size n.
 * @param      left_arr   The left array operand of size n of the operation.
 * @param      right_arr  The right array operand of size n of the operation.
 *
 */
#define VEC_HADAMARD(n, dst_arr, left_arr, right_arr)\
{\
	for (size_t i = 0; i < (n); i++)\
	{\
		(dst_arr)[i] = (left_arr)[i] * (right_arr)[i];\
	}\
}\

/**
 * @brief      Computes the dot product (inner product) between two vectors.
 *
 * @param      n          Dimensionality of the vectors.
 * @param      dst_arr    The destination array to contain the result of size n.
 * @param      left_arr   The left array operand of size n of the operation.
 * @param      right_arr  The right array operand of size n of the operation.
 *
 */
#define VEC_DOT(TYPE, n, left_arr, right_arr)\
{\
	TYPE dot = 0;\
	for (size_t i = 0; i < (n); i++)\
	{\
		dot += (left_arr)[i] * (right_arr)[i];\
	}\
	return dot;\
}\

/**
 * @brief      Computes the magnitude, or length of a vector.
 *
 * @param      TYPE  The type of the values in the vector
 * @param      n     Dimensionality of the vector.
 * @param      arr   The arr representing the vector
 *
 * @return     The magnitude of the vector.
 */
#define VEC_MAG(TYPE, n, arr)\
{\
	TYPE dot = 0;\
	for (size_t i = 0; i < (n); i++)\
	{\
		dot += (arr)[i] * (arr)[i];\
	}\
	return (TYPE)sqrt(dot);\
}\

static inline void vec_add(size_t n, TYPE dst[n], TYPE left[n], TYPE right[n])
{ VEC_ADD(n, dst, left, right) }

static inline void vec_sub(size_t n, TYPE dst[n], TYPE left[n], TYPE right[n])
{ VEC_SUB(n, dst, left, right) }

static inline void vec_scl(size_t n, TYPE dst[n], TYPE in[n], TYPE scalar)
{ VEC_SCL(n, dst, in, scalar) }

static inline void vec_hadamard(size_t n, TYPE dst[n], TYPE left[n], TYPE right[n])
{ VEC_HADAMARD(n, dst, left, right) }

static inline TYPE vec_dot(size_t n, TYPE left[n], TYPE right[n])
{ VEC_DOT(TYPE, n, left, right) }

static inline TYPE vec_mag(size_t n, TYPE left[n])
{ VEC_MAG(TYPE, n, left, left) }


#endif
