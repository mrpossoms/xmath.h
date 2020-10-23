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
 * @brief      Element wise division with the left vector elements as numerators
 *             and right vector elements as denominators
 *
 * @param      n          Dimensionality of the vectors.
 * @param      dst_arr    The destination array to contain the result of size n.
 * @param      left_arr   The left array operand of size n of the operation.
 * @param      right_arr  The right array operand of size n of the operation.
 *
 */
#define VEC_DIV(n, dst_arr, left_arr, right_arr)\
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
#include <string>

namespace xmath
{

template<size_t D, typename S=TYPE>
struct vec
{
	vec()
	{
		for (auto i = D; i--;) { v[i] = {}; }
	}

	vec(const S* arr)
	{
		for (auto i = D; i--;) { v[i] = arr[i]; }
	}

	vec(std::initializer_list<S> init)
	{
		if (init.size() < D)
		{
			for (auto i = D; i--;)
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
		for (auto i = D; i--;)
		{
			if (!isfinite(v[i])) { return false; }
		}

		return true;
	}

	inline vec<D,S> operator+(const vec<D,S>& v) const
	{
		vec<D,S> out;
		VEC_ADD(D, out.v, this->v, v.v)
		return out;
	}


	inline vec<D,S> operator-(const vec<D,S>& v) const
	{
		vec<D,S> out;
		VEC_SUB(D, out.v, this->v, v.v)
		return out;
	}


	inline vec<D,S> operator*(const vec<D,S>& v) const
	{
		vec<D,S> out;
		VEC_HADAMARD(D, out.v, this->v, v.v)
		return out;
	}


	inline vec<D,S> operator*(const S s) const
	{
		vec<D,S> out;
		VEC_SCL(D, out.v, this->v, s)
		return out;
	}


	inline vec<D,S>  operator/(const vec<D,S>& v) const
	{
		vec<D,S> out;
		VEC_DIV(D, out.v, this->v, v.v)
		return out;
	}


	inline vec<D,S> operator/(const S s) const
	{
		vec<D,S> out;
		VEC_SCL(D, out.v, this->v, 1.f/s)
		return out;
	}


	inline vec<D,S>& operator=(const vec<D,S>& v)
	{
		for (auto i = D; i--;)
		{
			this->v[i] = v.v[i];
		}
		return *this;
	}


	inline vec<D,S>& operator+=(const vec<D,S>& v) { return *this = *this + v; }

	inline vec<D,S>& operator-=(const vec<D,S>& v) { return *this = *this - v; }

	inline vec<D,S>& operator*=(const vec<D,S>& v) { return *this = *this * v; }

	inline vec<D,S>& operator*=(const S s) { return *this = *this * s; }

	inline vec<D,S>& operator/=(const vec<D,S>& v) { return *this = *this / v; }

	inline vec<D,S>& operator/=(const S s) { return *this = *this / s; }

	inline bool operator!=(const vec<D,S>& v) const { return !(*this == v); }

	inline bool operator==(const vec<D,S>& v) const
	{
		for (auto i = D; i--;)
		{
			if (v.at(i) != this->at(i)) { return false; }
		}

		return true;
	}

	template<typename T>
	inline vec<D,T> cast() const
	{
		vec<D,T> v;
		for (size_t i = 0; i < D; ++i)
		{
			v[i] = (T)this->v[i];
		}

		return v;
	}

	template<size_t ND>
	vec<ND,S> slice(size_t start = 0)
	{
		vec<ND,S> r;
		for (size_t i = 0; i < ND; ++i) { r[i] = v[i + start]; }
		return r;
	}


	std::string to_string(bool parens=true) const
	{
		std::string str = parens ? "(" : "";
		for (size_t i = 0; i < D; ++i)
		{
			str += std::to_string(v[i]);
			if (i < D - 1) { str += ", "; }
		} str += parens ? ")" : "";

		return str;
	}


	inline S magnitude() const { VEC_MAG(S, D, this->v) }


	vec<D,S> unit() const { return *this / magnitude(); }

	S dot(vec<D,S> const& v) const { VEC_DOT(S, D, this->v, v.v) }

	vec<D,S> project_onto(vec<D,S> const& v) const
	{
		auto v_hat = v.norm();
		auto len = this->dot(v_hat);
		return *this - (v_hat * len);
	}

	vec<D,S> lerp(const vec<D,S>& to, S p)
	{
		return (*this * (static_cast<S>(1) - p)) + (to * p);
	}

	bool is_near(const vec<D,S>& v, S threshold=0.0001)
	{
		auto diff = *this - v;
		return diff.dot(diff) <= threshold;
	}

	vec<D,S>& take_min(const vec<D,S>& v)
	{
		for (size_t i = 0; i < D; ++i)
		{
			auto& cur = this->v[i];
			cur = v[i] < cur ? v[i] : cur;
		}

		return *this;
	}


	vec<D,S>& take_max(const vec<D,S>& v)
	{
		for (size_t i = 0; i < D; ++i)
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

	S v[D]; // value store
};

} // namespace xmath end
#endif

#endif
