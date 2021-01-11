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
#define TYPE double
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
	for (size_t i = 0; i < (n); i++)\
	{\
		(dst_vec)[i] = (left_vec)[i] + (right_vec)[i];\
	}\
}\

/**
 * @brief      Compute vector subtraction between two vectors
 *
 * @param      n          Dimensionality of the vectors
 * @param      dst_vec    The destination array to contain the result of size n
 * @param      left_vec   The left array operand of the subtraction of size n
 * @param      right_vec  The right array operand of the subtraction of size n
 *
 */
#define VEC_SUB(n, dst_vec, left_vec, right_vec) \
{\
	for (size_t i = 0; i < (n); i++)\
	{\
		(dst_vec)[i] = (left_vec)[i] - (right_vec)[i];\
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
	for (size_t i = 0; i < (n); i++)\
	{\
		(dst_vec)[i] = (left_vec)[i] * (right_scalar);\
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
	for (size_t i = 0; i < (n); i++)\
	{\
		(dst_vec)[i] = (left_vec)[i] * (right_vec)[i];\
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
	for (size_t i = 0; i < (n); i++)\
	{\
		(dst_vec)[i] = (left_vec)[i] * (right_vec)[i];\
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
	for (size_t i = 0; i < (n); i++)\
	{\
		dot += (left_vec)[i] * (right_vec)[i];\
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
	for (size_t i = 0; i < (n); i++)\
	{\
		dot += (vec)[i] * (vec)[i];\
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
    for (int row = (m_R); row--;)\
    for (int col = (n_C); col--;)\
    {\
        for (int i = (m_C); i--;)\
        {\
            (r)[row][col] += (m)[row][i] * (n)[i][col];\
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
    for (int row = (mn_R); row--;)\
    for (int col = (mn_C); col--;)\
    {\
        (r)[row][col] = (m)[row][col] + (n)[row][col];\
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
    for (int row = (mn_R); row--;)\
    for (int col = (mn_C); col--;)\
    {\
        (r)[row][col] = (m)[row][col] - (n)[row][col];\
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
    for (TYPE row = (m_R); row--;)\
    { \
    	(r)[row] = 0;\
	    for (TYPE col = (n_C); col--;)\
	    {\
	        (r)[row] += (m)[row][col] * (n)[col];\
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
    for (TYPE row = (m_R); row--;)\
    { \
	    for (TYPE col = (n_C); col--;)\
	    {\
	        (r)[row][col] = (m)[row][col] * (s);\
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
    for (TYPE row = (m_R); row--;)\
    { \
	    for (TYPE col = (n_C); col--;)\
	    {\
	        (r)[row][col] = (m)[row][col] + (n)[row][col];\
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
    for (TYPE row = (m_R); row--;)\
    { \
	    for (TYPE col = (n_C); col--;)\
	    {\
	        (r)[row][col] = (m)[row][col] + (n)[row][col];\
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
	for (size_t ci = 0; ci < (C); ci++) \
	{ \
		TYPE t = M[(ri)][ci];\
		M[(ri)][ci] = M[(rj)][ci];\
		M[(rj)][ci] = t;\
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
	for (size_t r = 0; r < (R); r++)\
	{\
		if (M[r][piv_c] == 0)\
		{\
			ssize_t swap_ri = -1;\
			for (size_t ri = r + 1; ri < (R); ri++)\
			{\
				if (M[ri][piv_c] != 0)\
				{\
					swap_ri = ri;\
					break;\
				}\
			}\
			if (swap_ri > -1)\
			{\
				MAT_SWAP_ROWS(TYPE, C, M, swap_ri, r);\
			}\
		}\
		{\
			TYPE d = 1 / M[r][piv_c];\
			for (size_t c = piv_c; c < (C); c++) { M[r][c] *= d; }\
		}\
		for (size_t ri = 0; ri < (R); ri++)\
		{\
			if (M[ri][piv_c] == 0 || ri == r) { continue; }\
			TYPE d = M[ri][piv_c];\
			for (size_t c = piv_c; c < (C); c++)\
			{\
				M[ri][c] -= d * M[r][c];\
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
	for (size_t r = (R); r--;)\
	{\
		for (size_t c = (C); c--;) { aug[r][c + (C)] = c == r ? 1 : 0; }\
		for (size_t c = (C); c--;) { aug[r][c] = (M)[r][c]; }\
	}\
	MAT_RREF(TYPE, (R), (C) << 1, aug);\
	for (size_t r = (R); r--;)\
	for (size_t c = (C); c--;)\
	{\
		(M)[r][c] = aug[r][c + (C)];\
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
	for (size_t r = 0; r < R; r++)\
	{\
		for (size_t c = 0; c < C; c++)\
		{\
			OUT[c][r] = IN[r][c];\
		}\
	}\
}\


#ifndef __cplusplus
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
{ VEC_MAG(TYPE, n, left) }


#endif

#ifdef __cplusplus
#include <initializer_list>
#include <functional>
#include <string>

namespace xmath
{

template <size_t N, typename S=TYPE>
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

	inline S* ptr() { return v; }

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
	vec<NN,S> slice(size_t start = 0)
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

template <size_t R, size_t C, typename S=TYPE>
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

    void transpose_inplace()
    {
    	//MAT_INV_IMP(S, R, C, m)
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

	vec<C, S> m[R];
};

} // namespace xmath end
#endif
#endif
