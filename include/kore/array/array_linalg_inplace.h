#pragma once

#include "array.h"

namespace kore {
namespace array {
namespace linalg {
namespace inplace {

lapack_int inverse(Array<float64_t, 2> & mat);
lapack_int inverse(Array<complex128_t, 2> & mat);

// READ WRITE ARRAYS
complex128_t determinant(Array<complex128_t, 2> & mat);
float64_t    determinant(Array<float64_t, 2> & mat);

complex128_t log_determinant(Array<complex128_t, 2> & mat);
complex128_t log_determinant(Array<float64_t, 2> & mat);

float64_t    log_determinant_abs(Array<complex128_t, 2> & mat);
float64_t    log_determinant_abs(Array<float64_t, 2> & mat);

// compute eigenvalues and eigenvectors.
// eigenvectors will be stored in the original matrix

lapack_int eigensystem_symmetric(Array<float64_t, 2> & matrix_eigvecs,
                                 Array<float64_t, 1> & eigvals);
lapack_int eigensystem_hermitian(Array<complex128_t, 2> & matrix_eigvecs,
                                 Array<float64_t, 1> & eigvals);


} // namespace inplace
} // namespace linalg
} // namespace array
} // namespace kore

#include "array_linalg_inplace_impl.h"
