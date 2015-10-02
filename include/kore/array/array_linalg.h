#pragma once
#include <complex>
#include <tuple>
#include "array.h"

namespace kore {

namespace array {

namespace linalg {

std::tuple<Array<kore::float64_t, 1>, Array<kore::float64_t, 2> >
eigensystem_hermitian(const Array<kore::float64_t, 2>& matrix);

std::tuple<Array<kore::float64_t, 1>, Array<kore::complex128_t, 2> >
eigensystem_hermitian(const Array<kore::complex128_t, 2>& matrix);


/// ----- Transpose
template <typename ValueType>
Array<ValueType, 2> transpose(const Array<ValueType, 2>& arr);

//template <typename ValueType>
//Array<ValueType, 2> transpose(Array<ValueType, 2>&& arr);

template <typename ValueType>
Array<ValueType, 2> transpose(const ConstArray<ValueType, 2>& arr);


/// ----- Inverse
Array<kore::float64_t, 2> inverse(const Array<kore::float64_t, 2>& mat);
Array<kore::complex128_t, 2> inverse(const Array<kore::complex128_t, 2>& mat);

/// ----- Determinant

//! Compute determinant.
//! Destroy the contents inside.
inline
kore::complex128_t determinant(Array<kore::complex128_t, 2>& mat);

inline
kore::complex128_t log_determinant(Array<kore::complex128_t, 2>& mat);

//! Compute determinant.
//! Keep the contents inside.
inline
kore::complex128_t determinant(const Array<kore::complex128_t, 2>& mat);

//! Compute determinant.
//! Keep the contents inside.
inline
kore::complex128_t determinant(const ConstArray<kore::complex128_t, 2>& mat);

//! Compute determinant.
//! Keep the contents inside.
inline
kore::complex128_t log_determinant(const Array<kore::complex128_t, 2>& mat);

//! Compute determinant.
//! Keep the contents inside.
inline
kore::complex128_t log_determinant(const ConstArray<kore::complex128_t, 2>& mat);


//! Compute determinant.
//! Destroy the contents inside.
inline
kore::float64_t determinant(Array<kore::float64_t, 2>& mat);

//! Compute determinant.
//! Keep the contents inside.
inline
kore::float64_t determinant(const Array<kore::float64_t, 2>& mat);

//! Compute determinant.
//! Keep the contents inside.
inline
kore::float64_t determinant(const ConstArray<kore::float64_t, 2>& mat);


} // namespace linalg

} // namespace array

} // namespace kore

#include "array_linalg_impl.h"
