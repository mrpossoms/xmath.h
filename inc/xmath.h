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
#include <sys/types.h>

#ifndef XMTYPE
#define XMTYPE double
#endif

/**
 * @brief      Compute vector addition between two vectors
 *
 * @param      n          Dimensionality of the vectors
 * @param      dst_vec    The destination array to contain the result of size n
 * @param      left_vec   The left array operand of the addition of size n
 * @param      right_vec  The right array operand of the addition of size n
 *
 */
#define VEC_ADD(n, dst_vec, left_vec, right_vec) \
{\
	for (size_t __i = 0; __i < (n); __i++)\
	{\
		(dst_vec)[__i] = (left_vec)[__i] + (right_vec)[__i];\
	}\
}\

/**
 * @brief      Compute vector subtraction between two vectors
 *
 * @param      n          Dimensionality of the vectors
 * @param      dst_vec    The destination array to contain the result of size n
 * @param      left_vec   The left array operand of the subtraction of size n
 * @param      right_vec  The right array operand of the subtraction of size n

 */
#define VEC_SUB(n, dst_vec, left_vec, right_vec) \
{\
	for (size_t __i = 0; __i < (n); __i++)\
	{\
		(dst_vec)[__i] = (left_vec)[__i] - (right_vec)[__i];\
	}\
}\

/**
 * @brief      Scales each element of a vector by a scalar
 *
 * @param      n          Dimensionality of the vectors
 * @param      dst_vec    The destination array to contain the result of size n
 * @param      left_vec   The left array operand of the addition of size n
 * @param      right_scalar  The scalar by which each element is multiplied.
 *
 */
#define VEC_SCL(n, dst_vec, left_vec, right_scalar) \
{\
	for (size_t __i = 0; __i < (n); __i++)\
	{\
		(dst_vec)[__i] = (left_vec)[__i] * (right_scalar);\
	}\
}\

/**
 * @brief      Scales each element of a vector by a scalar
 *
 * @param      n          Dimensionality of the vectors
 * @param      dst_vec    The destination array to contain the result of size n
 * @param      left_vec   The left array operand of the addition of size n
 * @param      right_scalar  The scalar by which each element is multiplied.
 *
 */

/**
 * @brief Performs element-wise addition of left and right vectors after the
 * right vector has been scaled by a scalar operand.
 *
 * @param      n          Dimensionality of the vectors
 * @param      dst_vec    The destination array to contain the result of size n.
 * @param      left_vec   The left array operand of size n of the operation.
 * @param      right_vec  The right array operand of size n of the operation.
 * @param      scalar     The scalar by which each element of the right_vec is
 * multiplied.
 */
#define VEC_ADD_SCL(n, dst_vec, left_vec, right_vec, scalar) \
{\
	for (size_t __i = 0; __i < (n); __i++)\
	{\
		(dst_vec)[__i] = (left_vec)[__i] + ((right_vec)[__i] * (scalar));\
	}\
}\

/**
 * @brief      Computes the Hadamard product (element wise product) between two
 *             vectors.
 *
 * @param      n          Dimensionality of the vectors.
 * @param      dst_vec    The destination array to contain the result of size n.
 * @param      left_vec   The left array operand of size n of the operation.
 * @param      right_vec  The right array operand of size n of the operation.
 *
 */
#define VEC_HADAMARD(n, dst_vec, left_vec, right_vec)\
{\
	for (size_t __i = 0; __i < (n); __i++)\
	{\
		(dst_vec)[__i] = (left_vec)[__i] * (right_vec)[__i];\
	}\
}\

/**
 * @brief      Element wise division with the left vector elements as numerators
 *             and right vector elements as denominators
 *
 * @param      n          Dimensionality of the vectors.
 * @param      dst_vec    The destination array to contain the result of size n.
 * @param      left_vec   The left array operand of size n of the operation.
 * @param      right_vec  The right array operand of size n of the operation.
 *
 */
#define VEC_DIV(n, dst_vec, left_vec, right_vec)\
{\
	for (size_t __i = 0; __i < (n); __i++)\
	{\
		(dst_vec)[__i] = (left_vec)[__i] * (right_vec)[__i];\
	}\
}\

/**
 * @brief      Computes the dot product (inner product) between two vectors.
 *
 * @param      n          Dimensionality of the vectors.
 * @param      dst_vec    The destination array to contain the result of size n.
 * @param      left_vec   The left array operand of size n of the operation.
 * @param      right_vec  The right array operand of size n of the operation.
 *
 */
#define VEC_DOT(TYPE, n, left_vec, right_vec)\
{\
	TYPE dot = 0;\
	for (size_t __i = 0; __i < (n); __i++)\
	{\
		dot += (left_vec)[__i] * (right_vec)[__i];\
	}\
	return dot;\
}\

/**
 * @brief      Computes the magnitude, or length of a vector.
 *
 * @param      TYPE  The type of the values in the vector
 * @param      n     Dimensionality of the vector.
 * @param      vec   The vec representing the vector
 *
 * @return     The magnitude of the vector.
 */
#define VEC_MAG(TYPE, n, vec)\
{\
	TYPE dot = 0;\
	for (size_t __i = 0; __i < (n); __i++)\
	{\
		dot += (vec)[__i] * (vec)[__i];\
	}\
	return (TYPE)sqrt(dot);\
}\

/**
 * @brief      Performs a matrix multiplication between two matrices which are
 *             compatible.
 *
 * @param      TYPE  The type
 * @param      m_R   Rows in matrix 'm'
 * @param      m_C   Columns in matrix 'm'
 * @param      n_R   Rows in matrix 'n'
 * @param      n_C   Columns in matrix 'n'
 * @param      r     Resulting matrix of size (m_R, n_C)
 * @param      m     Left hand operand of size (m_R, m_C)
 * @param      n     Right hand operand of size (n_R, n_C)
 */
#define MAT_MUL(TYPE, m_R, m_C, n_R, n_C, r, m, n)\
{\
	for (int __row = (m_R); __row--;)\
	for (int __col = (n_C); __col--;)\
	{\
		for (int __i = (m_C); __i--;)\
		{\
			(r)[__row][__col] += (m)[__row][__i] * (n)[__i][__col];\
		}\
	}\
}\


/**
 * @brief      Performs a matrix addition. Each element of 'm' is added to its
 *             corresponding element in the matrix 'n' and stored in 'r'.
 *
 * @param      TYPE  The type
 * @param      mn_R   Rows in matrix 'm' and in 'n'
 * @param      mn_C   Columns in matrix 'm' and in 'n'
 * @param      r     Resulting matrix of size (mn_R, mn_C)
 * @param      m     Left hand operand of size (mn_R, mn_C)
 * @param      n     Right hand operand of size (mn_R, mn_C)
 */
#define MAT_ADD(TYPE, mn_R, mn_C, r, m, n)\
{\
	for (int __row = (mn_R); __row--;)\
	for (int __col = (mn_C); __col--;)\
	{\
		(r)[__row][__col] = (m)[__row][__col] + (n)[__row][__col];\
	}\
}\


/**
 * @brief      Performs a matrix subtraction. Each element of 'm' is subtracted
 *             from its corresponding element in the matrix 'n' and
 *             stored in 'r'.
 *
 * @param      TYPE  The type
 * @param      mn_R   Rows in matrix 'm' and in 'n'
 * @param      mn_C   Columns in matrix 'm' and in 'n'
 * @param      r     Resulting matrix of size (mn_R, mn_C)
 * @param      m     Left hand operand of size (mn_R, mn_C)
 * @param      n     Right hand operand of size (mn_R, mn_C)
 */
#define MAT_SUB(TYPE, mn_R, mn_C, r, m, n)\
{\
	for (int __row = (mn_R); __row--;)\
	for (int __col = (mn_C); __col--;)\
	{\
		(r)[__row][__col] = (m)[__row][__col] - (n)[__row][__col];\
	}\
}\


/**
 * @brief      Performs matrix-vector multiplication resulting in another
 *             vector.
 *
 * @param      TYPE  The type
 * @param      m_R   Rows in matrix 'm'
 * @param      m_C   Columns in matrix 'm'
 * @param      r     Resulting vector of size m_R.
 * @param      m     Left hand operand matrix of size (m_R, m_C).
 * @param      v     Right hand operand vector of size m_C.
 */
#define MAT_MUL_VEC(TYPE, m_R, m_C, r, m, v)\
{\
	for (TYPE __row = (m_R); __row--;)\
	{ \
		(r)[__row] = 0;\
		for (TYPE __col = (n_C); __col--;)\
		{\
			(r)[__row] += (m)[__row][__col] * (n)[__col];\
		}\
	{\
}\


/**
 * @brief      Scales each element of a matrix with a scalar and stores the
 *             result in another matrix of matching size.
 *
 * @param      TYPE  The type
 * @param      m_R   Rows in matrix 'm'
 * @param      m_C   Columns in matrix 'm'
 * @param      r     Resulting matrix with scale applied of size (m_R, m_C)
 * @param      m     Left hand operand matrix of size (m_R, m_C)
 * @param      s     Scalar which will be multiplied agains each element of 'm'
 */
#define MAT_MUL_E(TYPE, m_R, m_C, r, m, s)\
{\
	for (TYPE __row = (m_R); __row--;)\
	{ \
		for (TYPE __col = (n_C); __col--;)\
		{\
			(r)[__row][__col] = (m)[__row][__col] * (s);\
		}\
	{\
}\


/**
 * @brief      Performs an elementwise addition between two matrices
 *
 * @param      TYPE  The type
 * @param      m_R   Rows in matrix 'm'
 * @param      m_C   Columns in matrix 'm'
 * @param      r     Resulting matrix of size (m_R, m_C)
 * @param      m     Left hand operand in the addition of size (m_R, m_C)
 * @param      n     Right hand operand in the addition of size (m_R, m_C)
 */
#define MAT_ADD_E(TYPE, m_R, m_C, r, m, n)\
{\
	for (TYPE __row = (m_R); __row--;)\
	{ \
		for (TYPE __col = (n_C); __col--;)\
		{\
			(r)[__row][__col] = (m)[__row][__col] + (n)[__row][__col];\
		}\
	{\
}\


/**
 * @brief      Performs an elementwise subtraction between two matrices
 *
 * @param      TYPE  The type
 * @param      m_R   Rows in matrix 'm'
 * @param      m_C   Columns in matrix 'm'
 * @param      r     Resulting matrix of size (m_R, m_C)
 * @param      m     Left hand operand in the subtraction of size (m_R, m_C)
 * @param      n     Right hand operand in the subtraction of size (m_R, m_C)
 */
#define MAT_SUB_E(TYPE, m_R, m_C, r, m, n)\
{\
	for (TYPE __row = (m_R); __row--;)\
	{ \
		for (TYPE __col = (n_C); __col--;)\
		{\
			(r)[__row][__col] = (m)[__row][__col] + (n)[__row][__col];\
		}\
	{\
}\


/**
 * @brief      Swaps two rows of a matrix
 *
 * @param      TYPE  The storage type for each element.
 * @param      C     Number of columns in each row.
 * @param      M     Matrix (2d array) whose rows should be swapped.
 * @param      ri    Index of the first row to swap.
 * @param      rj    Index of the second row to swap.
 */
#define MAT_SWAP_ROWS(TYPE, C, M, ri, rj)\
{\
	for (size_t __ci = 0; __ci < (C); __ci++) \
	{ \
		TYPE t = M[(ri)][__ci];\
		M[(ri)][__ci] = M[(rj)][__ci];\
		M[(rj)][__ci] = t;\
	} \
}\


/**
 * @brief      Transforms augmented matrix to row-reduced echelon form.
 *
 * @param      TYPE  The storage type for each element.
 * @param      R     Number of rows in the matrix
 * @param      C     Number of columns in each row.
 * @param      M     Matrix (2d array) whose rows should be swapped.
 */
#define MAT_RREF(TYPE, R, C, M)\
{\
	size_t piv_c = 0;\
	for (size_t __r = 0; __r < (R); __r++)\
	{\
		if (M[__r][piv_c] == 0)\
		{\
			ssize_t swap_ri = -1;\
			for (size_t __ri = __r + 1; __ri < (R); __ri++)\
			{\
				if (M[__ri][piv_c] != 0)\
				{\
					swap_ri = __ri;\
					break;\
				}\
			}\
			if (swap_ri > -1)\
			{\
				MAT_SWAP_ROWS(TYPE, C, M, swap_ri, __r);\
			}\
		}\
		{\
			TYPE d = 1 / M[__r][piv_c];\
			for (size_t __c = piv_c; __c < (C); __c++) { M[__r][__c] *= d; }\
		}\
		for (size_t __ri = 0; __ri < (R); __ri++)\
		{\
			if (M[__ri][piv_c] == 0 || __ri == __r) { continue; }\
			TYPE d = M[__ri][piv_c];\
			for (size_t __c = piv_c; __c < (C); __c++)\
			{\
				M[__ri][__c] -= d * M[__r][__c];\
			}\
		}\
		++piv_c;\
	}\
}\


/**
 * @brief      Computes the inverse of an invertible matrix.
 *
 * @param      TYPE  The storage type for each element.
 * @param      R     Number of rows in the matrix
 * @param      C     Number of columns in each row.
 * @param      M     Matrix (2d array) whose rows should be swapped.
 */
#define MAT_INV_IMP(TYPE, R, C, M)\
{\
	TYPE aug[(R)][(C) << 1];\
	for (size_t __r = (R); __r--;)\
	{\
		for (size_t __c = (C); __c--;) { aug[__r][__c + (C)] = __c == __r ? 1 : 0; }\
		for (size_t __c = (C); __c--;) { aug[__r][__c] = (M)[__r][__c]; }\
	}\
	MAT_RREF(TYPE, (R), (C) << 1, aug);\
	for (size_t __r = (R); __r--;)\
	for (size_t __c = (C); __c--;)\
	{\
		(M)[__r][__c] = aug[__r][__c + (C)];\
	}\
}\

/**
 * @brief      Computes the transpose of a matrix.
 *
 * @param      TYPE   The storage type for each element.
 * @param      R      Number of rows in the matrix
 * @param      C      Number of columns in each row.
 * @param      IN     Original matrix (2d array).
 * @param      OUT    The trasnpose of the original matrix.
 */
#define MAT_TRANSPOSE(TYPE, R, C, IN, OUT)\
{\
	for (size_t __r = 0; __r < R; __r++)\
	{\
		for (size_t __c = 0; __c < C; __c++)\
		{\
			OUT[__c][__r] = IN[__r][__c];\
		}\
	}\
}\

/**
 * @brief      Sets an identity matrix.
 *
 * @param      TYPE   The storage type for each element.
 * @param      R      Number of rows in the matrix
 * @param      C      Number of columns in each row.
 * @param      MAT    The matrix in question.
 */
#define MAT_IDENTITY(R, C, MAT)\
{\
	for (size_t __r = 0; __r < (R); __r++)\
	{\
		for (size_t __c = 0; __c < (C); __c++)\
		{\
			if (__r == __c) (MAT)[__r][__c] = 1;\
			else (MAT)[__r][__c] = 0;\
		}\
	}\
}\

#ifndef __cplusplus
static inline void vec_add(size_t n, XMTYPE dst[n], const XMTYPE left[n], const XMTYPE right[n])
{ VEC_ADD(n, dst, left, right) }

static inline void vec_sub(size_t n, XMTYPE dst[n], const XMTYPE left[n], const XMTYPE right[n])
{ VEC_SUB(n, dst, left, right) }

static inline void vec_scl(size_t n, XMTYPE dst[n], const XMTYPE in[n], XMTYPE scalar)
{ VEC_SCL(n, dst, in, scalar) }

static inline void vec_add_scl(size_t n, XMTYPE dst[n], const XMTYPE left[n], const XMTYPE right[n], XMTYPE s)
{ VEC_ADD_SCL(n, dst, left, right, s) }

static inline void vec_hadamard(size_t n, XMTYPE dst[n], const XMTYPE left[n], const XMTYPE right[n])
{ VEC_HADAMARD(n, dst, left, right) }

static inline XMTYPE vec_dot(size_t n, const XMTYPE left[n], const XMTYPE right[n])
{ VEC_DOT(XMTYPE, n, left, right) }

static inline XMTYPE vec_mag(size_t n, const XMTYPE left[n])
{ VEC_MAG(XMTYPE, n, left) }

static inline void mat_transpose(size_t r, size_t c, const XMTYPE in[r][c], XMTYPE out[c][r])
{ MAT_TRANSPOSE(XMTYPE, r, c, in, out) }

static inline void mat_mul(size_t m_R, size_t m_C,  size_t n_R, size_t n_C, XMTYPE r[m_R][n_C], const XMTYPE m[m_R][m_C], const XMTYPE n[n_R][n_C])
{ MAT_MUL(XMTYPE, m_R, m_C, n_R, n_C, r, m, n) }

static inline void mat_add(size_t m_R, size_t m_C, XMTYPE r[m_R][m_C], const XMTYPE m[m_R][m_C], const XMTYPE n[m_R][m_C])
{ MAT_ADD(XMTYPE, m_R, m_C, r, m, n) }

static inline void mat_sub(size_t m_R, size_t m_C, XMTYPE r[m_R][m_C], const XMTYPE m[m_R][m_C], const XMTYPE n[m_R][m_C])
{ MAT_SUB(XMTYPE, m_R, m_C, r, m, n) }

static inline void mat_inv(size_t r, size_t c, const XMTYPE in[r][c], XMTYPE out[r][c])
{
	for (size_t i = 0; i < r; i++)
	{
		for (size_t j = 0; j < c; j++)
		{
			out[i][j] = in[i][j];
		}
	}
	MAT_INV_IMP(XMTYPE, r, c, out)
}

static inline void mat_identity(size_t r, size_t c, XMTYPE mat[r][c])
{ MAT_IDENTITY(r, c, mat) }
#endif

#ifdef __cplusplus
#include <initializer_list>
#include <functional>
#include <string>

namespace xmath
{

template <size_t N, typename S=XMTYPE>
struct vec
{
	vec()
	{
		for (auto i = N; i--;) { v[i] = {}; }
	}

	vec(const S* arr)
	{
		for (auto i = N; i--;) { v[i] = arr[i]; }
	}

	vec(std::initializer_list<S> init)
	{
		if (init.size() < N)
		{
			for (auto i = N; i--;)
			{
				v[i] = *init.begin();
			}
		}
		else
		{
			auto i = 0;
			for (auto s : init)
			{
				v[i++] = s;
			}
		}
	}

	inline S operator[](int i) const { return v[i]; }

	inline S& operator[](int i) { return v[i]; }

	inline S at(int i) const { return v[i]; }

	inline S* ptr() const { return v; }

	inline bool is_finite() const
	{
		for (auto i = N; i--;)
		{
			if (!isfinite(v[i])) { return false; }
		}

		return true;
	}

	inline vec<N,S> operator+(const vec<N,S>& v) const
	{
		vec<N,S> out;
		VEC_ADD(N, out.v, this->v, v.v)
		return out;
	}


	inline vec<N,S> operator-(const vec<N,S>& v) const
	{
		vec<N,S> out;
		VEC_SUB(N, out.v, this->v, v.v)
		return out;
	}


	inline vec<N,S> operator*(const vec<N,S>& v) const
	{
		vec<N,S> out;
		VEC_HADAMARD(N, out.v, this->v, v.v)
		return out;
	}


	inline vec<N,S> operator*(const S s) const
	{
		vec<N,S> out;
		VEC_SCL(N, out.v, this->v, s)
		return out;
	}


	inline vec<N,S>  operator/(const vec<N,S>& v) const
	{
		vec<N,S> out;
		VEC_DIV(N, out.v, this->v, v.v)
		return out;
	}


	inline vec<N,S> operator/(const S s) const
	{
		vec<N,S> out;
		VEC_SCL(N, out.v, this->v, 1.f/s)
		return out;
	}


	inline vec<N,S>& operator=(const vec<N,S>& v)
	{
		for (auto i = N; i--;)
		{
			this->v[i] = v.v[i];
		}
		return *this;
	}


	inline vec<N,S>& operator+=(const vec<N,S>& v) { return *this = *this + v; }

	inline vec<N,S>& operator-=(const vec<N,S>& v) { return *this = *this - v; }

	inline vec<N,S>& operator*=(const vec<N,S>& v) { return *this = *this * v; }

	inline vec<N,S>& operator*=(const S s) { return *this = *this * s; }

	inline vec<N,S>& operator/=(const vec<N,S>& v) { return *this = *this / v; }

	inline vec<N,S>& operator/=(const S s) { return *this = *this / s; }

	inline bool operator!=(const vec<N,S>& v) const { return !(*this == v); }

	inline bool operator==(const vec<N,S>& v) const
	{
		for (auto i = N; i--;)
		{
			if (v.at(i) != this->at(i)) { return false; }
		}

		return true;
	}

	template<typename T>
	inline vec<N,T> cast() const
	{
		vec<N,T> v;
		for (size_t i = 0; i < N; ++i)
		{
			v[i] = (T)this->v[i];
		}

		return v;
	}

	template<size_t NN>
	vec<NN,S> slice(size_t start = 0) const
	{
		vec<NN,S> r;
		for (size_t i = 0; i < NN; ++i) { r[i] = v[i + start]; }
		return r;
	}


	std::string to_string(bool parens=true) const
	{
		std::string str = parens ? "(" : "";
		for (size_t i = 0; i < N; ++i)
		{
			str += std::to_string(v[i]);
			if (i < N - 1) { str += ", "; }
		} str += parens ? ")" : "";

		return str;
	}


	inline S magnitude() const { VEC_MAG(S, N, this->v) }

	vec<N,S> unit() const { return *this / magnitude(); }

	S dot(vec<N,S> const& v) const { VEC_DOT(S, N, this->v, v.v) }

	vec<N,S> project_onto(vec<N,S> const& v) const
	{
		auto v_hat = v.norm();
		auto len = this->dot(v_hat);
		return *this - (v_hat * len);
	}

	vec<N,S> lerp(const vec<N,S>& to, S p)
	{
		return (*this * (static_cast<S>(1) - p)) + (to * p);
	}

	bool is_near(const vec<N,S>& v, S threshold=0.0001)
	{
		auto diff = *this - v;
		return diff.dot(diff) <= threshold;
	}

	vec<N,S>& take_min(const vec<N,S>& v)
	{
		for (size_t i = 0; i < N; ++i)
		{
			auto& cur = this->v[i];
			cur = v[i] < cur ? v[i] : cur;
		}

		return *this;
	}


	vec<N,S>& take_max(const vec<N,S>& v)
	{
		for (size_t i = 0; i < N; ++i)
		{
			auto& cur = this->v[i];
			cur = v[i] > cur ? v[i] : cur;
		}

		return *this;
	}

	static S cross(vec<2,S>& a, vec<2,S>& b)
	{
		return a[0]*b[1] - a[1]*b[0];
	}

	static vec<3,S> cross(vec<3,S> const& a, vec<3,S> const& b)
	{
		return {
			a[1]*b[2] - a[2]*b[1],
			a[2]*b[0] - a[0]*b[2],
			a[0]*b[1] - a[1]*b[0],
		};
	}

	S v[N]; // value store
};

template <size_t R, size_t C, typename S=XMTYPE>
struct mat
{
	mat()
	{
		for (auto row = R; row--;)
		for (auto col = C; col--;)
		{
			m[row][col] = {};
		}
	}

	mat(std::initializer_list<std::initializer_list<S>> init)
	{
		int ri = 0;
		for (auto row : init)
		{
			int ci = 0;
			for (auto c : row)
			{
				m[ri][ci] = c;
				ci += 1;
			}
			ri += 1;
		}

	}

	mat<R, C, S>& initialize(std::function<S (S r, S c)> init)
	{
		for (auto row = R; row--;)
		for (auto col = C; col--;)
		{
			m[row][col] = init(static_cast<S>(row), static_cast<S>(col));
		}

		return *this;
	}

	static inline mat<R, C, S> I()
	{
	mat<R, C, S> m;
	MAT_IDENTITY(R, C, m);
	return m;
	}

	void invert_inplace()
	{
		MAT_INV_IMP(S, R, C, m)
	}

	mat<R, C, S> invert()
	{
		mat<R, C, S> out = *this;
		MAT_INV_IMP(S, R, C, out);
		return out;
	}

	mat<C, R, S> transpose()
	{
		mat<C, R, S> out;
		MAT_TRANSPOSE(S, R, C, m, out);
	return out;
	}

	vec<C, S>& operator[](size_t r) { return m[r]; }

	mat<R, C, S> operator+ (const mat<R, C, S>& M)
	{
		mat<R, C, S> out;
		MAT_ADD(S, R, C, out, m, M.m);
		return out;
	}

	mat<R, C, S> operator- (const mat<R, C, S>& M)
	{
		mat<R, C, S> out;
		MAT_SUB(S, R, C, out, m, M.m);
		return out;
	}

	template<size_t O>
	mat<R, O, S> operator* (const mat<C, O, S>& N)
	{
		mat<R, O, S> out;
		MAT_MUL(S, R, C, C, O, out.m, m, N.m);
		return out;
	}

	const S* ptr() const { return m[0].v; }

	static mat<4, 4> rotation(vec<3> axis, float angle)
	{
		const auto a = axis;
		const auto c = cosf(angle);
		const auto s = sinf(angle);
		const auto omc = 1 - c;

		return {
			{c+a[0]*a[0]*omc,      a[1]*a[0]*omc+a[2]*s, a[2]*a[0]*omc-a[1]*s, 0},
			{a[0]*a[1]*omc-a[2]*s, c+a[1]*a[1]*omc,      a[2]*a[1]*omc+a[0]*s, 0},
			{a[0]*a[2]*omc+a[1]*s, a[1]*a[2]*omc-a[0]*s, c+a[2]*a[2]*omc,      0},
			{                   0,                    0,                    0, 1}
		};
	}

	static mat<4, 4> scale(vec<3> t)
	{
		return {
			{ t[0],    0,     0,    0    },
			{    0,  t[1],    0,    0    },
			{    0,    0,  t[2],    0    },
			{    0,    0,     0,    1.   }
		};
	}

	static mat<4, 4> translation(vec<3> t)
	{
		return {
			{    1,    0,    0,    0    },
			{    0,    1,    0,    0    },
			{    0,    0,    1,    0    },
			{  t[0], t[1], t[2],   1.   }
		};
	}

	static mat<4, 4> perspective(S near, S far, S fov, S aspect)
	{
		const auto a = tanf(M_PI * 0.5f - 0.5f * fov);
		const auto fsn = far - near;
		const auto fpn = far + near;
		const auto ftn = far * near;

		return {
			{  a/aspect,         0,          0,         0 },
			{         0,         a,          0,         0 },
			{         0,         0,   -fpn/fsn,        -1 },
			{         0,         0, -2*ftn/fsn,         0 }
		};
	}

	static mat<4, 4> orthographic(S near, S far, S left, S right, S top, S bottom)
	{
		const auto rml = right - left;
		const auto rpl = right + left;
		const auto tmb = top - bottom;
		const auto tpb = top + bottom;
		const auto fmn = far - near;
		const auto fpn = far + near;

		return {
			{2/rml,     0,      0,     0},
			{    0, 2/tmb,      0,     0},
			{    0,     0, -2/fmn,     0},
			{    0,     0,      0,     1},
		};
	}

	vec<C, S> m[R];
};

struct quat : public vec<4>
{
    quat() : vec({ 0, 0, 0, 1 })
    {
        // NOP
    }

    quat(const float* v) : vec({ v[0], v[1], v[2], v[3] })
    {
        // NOP
    }

    quat(float x, float y, float z, float w) : vec({ x, y, z, w })
    {
        // NOP
    }

    quat(vec<4> v) : vec(v)
    {
        // NOP
    }


    quat operator*(quat const& other) const
    {
        auto t3 = this->slice<3>(0);
        auto o3 = other.slice<3>(0);

        auto r = vec::cross(t3, o3);
        auto w = t3 * other[3];
        r += w;
        w = o3 * this->v[3];
        r += w;

        return {
            r[0], r[1], r[2],
            this->v[3] * other.v[3] - t3.dot(o3)
        };
    }


    quat operator*(float s) const
    {
        return { slice<4>(0) * s };
    }


    quat& operator*=(quat const& other)
    {
        *this = *this * std::move(other);
        return *this;
    }

    quat conjugate() const
    {
        auto& q = *this;
        return { -q[0], -q[1], -q[2], q[3] };
    }

    quat inverse() const
    {
        auto inv = conjugate();
        auto mag2 = dot(*this);
        static_cast<vec<4>>(inv) /= mag2;
        return inv;
    }

    inline float rotational_difference(quat const& q) const
    {
        auto q_d = q * this->inverse();
        return 2 * atan2(q_d.slice<3>(0).magnitude(), fabsf(q_d[3]));
    }

    quat slerp_to(quat const& p1, float t) const
    {
        const auto& p0 = *this;
        auto W = rotational_difference(p1);
        auto sin_W = sin(W);
        if (sin_W < 0.0001f) { sin_W = 0.0001f; }
        return p0 * (sin((1 - t) * W) / sin_W) + p1 * (sin(t * W) / sin_W);
    }

    vec<3> rotate(vec<3> const& v) const
    {
        auto q_xyz = this->slice<3>(0);

        auto t = vec::cross(q_xyz, v);
        t *= 2;

        auto u = vec::cross(q_xyz, t);
        t *= this->v[3];

        return v + t + u;
    }

    inline mat<4, 4> to_matrix() const
    {
        auto v = static_cast<vec<4>>(*this);
        auto a = v[3], b = v[0], c = v[1], d = v[2];
        auto a2 = a * a, b2 = b * b, c2 = c * c, d2 = d * d;

        return {
            { a2 + b2 - c2 - d2, 2 * (b*c - a*d)  , 2 * (b*d + a*c)  , 0},
            { 2 * (b*c + a*d)  , a2 - b2 + c2 - d2, 2 * (c*d - a*b)  , 0},
            { 2 * (b*d - a*c)  , 2 * (c*d + a*b)  , a2 - b2 - c2 + d2, 0},
            { 0                , 0                , 0                , 1},
        };
    }

    template <class S>
    vec<3, S> to_roll_pitch_yaw()
    {
        S roll, pitch, yaw;
        // roll (x-axis rotation)
        S sinr_cosp = +2.0 * (v[3] * v[0] + v[1] * v[2]);
        S cosr_cosp = +1.0 - 2.0 * (v[0] * v[0] + v[1] * v[1]);
        roll = atan2(sinr_cosp, cosr_cosp);

        // pitch (y-axis rotation)
        S sinp = +2.0 * (v[3] * v[1] - v[2] * v[0]);
        if (fabs(sinp) >= 1)
        {
            pitch = copysign(M_PI / 2, sinp); // use 90 degrees if out of range
        }
        else
        {
            pitch = asin(sinp);
        }

        // yaw (z-axis rotation)
        S siny_cosp = +2.0 * (v[3] * v[2] + v[0] * v[1]);
        S cosy_cosp = +1.0 - 2.0 * (v[1] * v[1] + v[2] * v[2]);
        yaw = atan2(siny_cosp, cosy_cosp);

        return { roll, pitch, yaw };
    }

    static quat from_axis_angle(vec<3> axis, float angle)
    {
        auto a_2 = angle / 2;
        auto a = sinf(a_2);

        axis *= a;

        return { axis[0], axis[1], axis[2], cosf(a_2) };
    }
};

} // namespace xmath end
#endif
#endif
