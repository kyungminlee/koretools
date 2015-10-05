#pragma once
#include <complex>
#include <tuple>
#include "array.h"

namespace kore {

namespace array {

namespace linalg {

//! \name Binary
//@{
template <typename T1, typename T2>
Array<typename std::remove_const<T1>::type, 2>
dot(Array<T1, 2> const & a,
    Array<T2, 2> const & b);
//@}

//! \name Unary
//@{
/// ----- Transpose
template <typename ValueType>
Array<typename std::remove_const<ValueType>::type, 2> transpose(Array<ValueType, 2> const & arr);

//template <typename ValueType>
//Array<ValueType, 2> transpose(Array<ValueType, 2>&& arr);

/// ----- Inverse
template <typename ValueType, typename>
Array<typename std::remove_const<ValueType>::type, 2> inverse(const Array<ValueType, 2>& mat);
//Array<kore::float64_t, 2> inverse(const Array<kore::float64_t, 2>& mat);
//Array<kore::complex128_t, 2> inverse(const Array<kore::complex128_t, 2>& mat);

/// ----- Determinant

//! Compute determinant.
//! Destroy the contents inside.


// READ ONLY ARRAYS
complex128_t determinant(Array<complex128_t const, 2> const & mat);
float64_t determinant(Array<float64_t const, 2> const& mat);

complex128_t log_determinant(Array<complex128_t const, 2> const & mat);
complex128_t log_determinant(Array<float64_t const, 2> const & mat);

float64_t log_determinant_abs(Array<complex128_t const, 2> const & mat);
float64_t log_determinant_abs(Array<float64_t const, 2> const & mat);


// copied
std::tuple<Array<float64_t, 1>, Array<float64_t, 2> >
eigensystem_symmetric(Array<float64_t const, 2> const & matrix);

std::tuple<Array<float64_t, 1>, Array<complex128_t, 2> >
eigensystem_hermitian(Array<complex128_t const, 2> const & matrix);



//@}

} // namespace linalg
} // namespace array
} // namespace kore

#include "array_linalg_impl.h"
