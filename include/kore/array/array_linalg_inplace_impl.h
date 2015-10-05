#pragma once

#include "array_linalg_inplace.h"

namespace kore {
namespace array {
namespace linalg {
namespace inplace {


inline lapack_int
inverse(Array<float64_t, 2> & mat)
{
  lapack_int m   = static_cast<lapack_int>(mat.shape(0));
  lapack_int n   = static_cast<lapack_int>(mat.shape(1));
  lapack_int lda = static_cast<lapack_int>(mat.stride(0));
  assert(static_cast<size_t>(m)   == mat.shape(0));
  assert(static_cast<size_t>(n)   == mat.shape(1));
  assert(static_cast<size_t>(lda) == mat.stride(0));
  assert(m == n);

  std::vector<lapack_int> ipiv(n);

  lapack_int info = 0;
  info = LAPACKE_dgetrf(LAPACK_ROW_MAJOR, m, n, mat.begin(), lda, ipiv.data());
  info = LAPACKE_dgetri(LAPACK_ROW_MAJOR, n, mat.begin(), lda, ipiv.data());
  return info;
}


inline lapack_int
inverse(Array<complex128_t, 2> & mat)
{
  lapack_int m   = static_cast<lapack_int>(mat.shape(0));
  lapack_int n   = static_cast<lapack_int>(mat.shape(1));
  lapack_int lda = static_cast<lapack_int>(mat.stride(0));
  assert(static_cast<size_t>(m)   == mat.shape(0));
  assert(static_cast<size_t>(n)   == mat.shape(1));
  assert(static_cast<size_t>(lda) == mat.stride(0));
  assert(m == n);
  std::vector<lapack_int> ipiv(n);

  lapack_int info = 0;
  info = LAPACKE_zgetrf(LAPACK_ROW_MAJOR, m, n,
    reinterpret_cast<lapack_complex_double*>(mat.begin()), lda, ipiv.data());
  info = LAPACKE_zgetri(LAPACK_ROW_MAJOR, n, reinterpret_cast<lapack_complex_double*>(mat.begin()), lda, ipiv.data());
  return info;
}

//! Compute determinant.
//! Destroy the contents inside.
inline float64_t
determinant(Array<float64_t, 2> & mat)
{
  lapack_int m = static_cast<lapack_int>(mat.shape(0));
  lapack_int n = static_cast<lapack_int>(mat.shape(1));
  lapack_int lda = static_cast<lapack_int>(mat.stride(0));
  assert(static_cast<size_t>(m) == mat.shape(0));
  assert(static_cast<size_t>(n) == mat.shape(1));
  assert(static_cast<size_t>(lda) == mat.stride(0));
  assert(m == n);

  std::vector<lapack_int> ipiv(n);
  LAPACKE_dgetrf(LAPACK_ROW_MAJOR, m, n, mat.begin(), lda, ipiv.data());
  auto sgn = false;
  for (lapack_int i = 0; i < n; ++i) {
    if (ipiv[i] != i + 1) { sgn = !sgn; }
  }
  float64_t det = (sgn ? -1.0 : 1.0);
  for (lapack_int i = 0; i < n; ++i) {
    det *= mat(i, i);
  }
  return det;
}


//! Compute determinant.
//! Destroy the contents inside.
inline complex128_t
determinant(Array<complex128_t, 2> & mat)
{
  lapack_int m = static_cast<lapack_int>(mat.shape(0));
  lapack_int n = static_cast<lapack_int>(mat.shape(1));
  lapack_int lda = static_cast<lapack_int>(mat.stride(0));
  assert(static_cast<size_t>(m) == mat.shape(0));
  assert(static_cast<size_t>(n) == mat.shape(1));
  assert(static_cast<size_t>(lda) == mat.stride(0));
  assert(m == n);

  std::vector<lapack_int> ipiv(n);
  LAPACKE_zgetrf(LAPACK_ROW_MAJOR, m, n,
    reinterpret_cast<lapack_complex_double*>(mat.begin()), lda, ipiv.data());
  int sgn = 0;
  for (lapack_int i = 0; i < n; ++i) {
    if (ipiv[i] != i + 1) { sgn = !sgn; }
  }
  complex128_t det = (sgn ? -1.0 : 1.0);
  for (lapack_int i = 0; i < n; ++i) {
    det *= mat(i, i);
  }
  return det;
}


inline complex128_t
log_determinant(Array<float64_t, 2> & mat)
{
  lapack_int m = static_cast<lapack_int>(mat.shape(0));
  lapack_int n = static_cast<lapack_int>(mat.shape(1));
  lapack_int lda = static_cast<lapack_int>(mat.stride(0));
  assert(static_cast<size_t>(m) == mat.shape(0));
  assert(static_cast<size_t>(n) == mat.shape(1));
  assert(static_cast<size_t>(lda) == mat.stride(0));
  assert(m == n);

  std::vector<lapack_int> ipiv(n);
  LAPACKE_dgetrf(LAPACK_ROW_MAJOR, m, n,
                 mat.begin(), lda, ipiv.data());
  int sgn = 0;
  for (lapack_int i = 0; i < n; ++i) {
    if (ipiv[i] != i + 1) { sgn = !sgn; }
  }
  kore::complex128_t log_det = (sgn ? kore::complex128_t(0.0, M_PI) : 0.0);
  for (lapack_int i = 0; i < n; ++i) {
    log_det += std::log(mat(i, i));
  }
  return log_det;
}


inline complex128_t
log_determinant(Array<complex128_t, 2> & mat)
{
  lapack_int m = static_cast<lapack_int>(mat.shape(0));
  lapack_int n = static_cast<lapack_int>(mat.shape(1));
  lapack_int lda = static_cast<lapack_int>(mat.stride(0));
  assert(static_cast<size_t>(m) == mat.shape(0));
  assert(static_cast<size_t>(n) == mat.shape(1));
  assert(static_cast<size_t>(lda) == mat.stride(0));
  assert(m == n);

  std::vector<lapack_int> ipiv(n);
  LAPACKE_zgetrf(LAPACK_ROW_MAJOR, m, n,
                 reinterpret_cast<lapack_complex_double*>(mat.begin()), lda, ipiv.data());
  int sgn = 0;
  for (lapack_int i = 0; i < n; ++i) {
    if (ipiv[i] != i + 1) { sgn = !sgn; }
  }
  kore::complex128_t log_det = (sgn ? kore::complex128_t(0.0, M_PI) : 0.0);
  for (lapack_int i = 0; i < n; ++i) {
    log_det += std::log(mat(i, i));
  }
  return log_det;
}


inline float64_t
log_determinant_abs(Array<float64_t, 2> & mat)
{
  lapack_int m = static_cast<lapack_int>(mat.shape(0));
  lapack_int n = static_cast<lapack_int>(mat.shape(1));
  lapack_int lda = static_cast<lapack_int>(mat.stride(0));
  assert(static_cast<size_t>(m) == mat.shape(0));
  assert(static_cast<size_t>(n) == mat.shape(1));
  assert(static_cast<size_t>(lda) == mat.stride(0));
  assert(m == n);

  std::vector<lapack_int> ipiv(n);
  LAPACKE_dgetrf(LAPACK_ROW_MAJOR, m, n,
                 mat.begin(), lda, ipiv.data());
  int sgn = 0;
  for (lapack_int i = 0; i < n; ++i) {
    if (ipiv[i] != i + 1) { sgn = !sgn; }
  }
  float64_t log_det = 0.0;
  for (lapack_int i = 0; i < n; ++i) {
    log_det += std::log(std::abs(mat(i, i)));
  }
  return log_det;
}


inline float64_t
log_determinant_abs(Array<complex128_t, 2> & mat)
{
  lapack_int m = static_cast<lapack_int>(mat.shape(0));
  lapack_int n = static_cast<lapack_int>(mat.shape(1));
  lapack_int lda = static_cast<lapack_int>(mat.stride(0));
  assert(static_cast<size_t>(m) == mat.shape(0));
  assert(static_cast<size_t>(n) == mat.shape(1));
  assert(static_cast<size_t>(lda) == mat.stride(0));
  assert(m == n);

  std::vector<lapack_int> ipiv(n);
  LAPACKE_zgetrf(LAPACK_ROW_MAJOR, m, n,
    reinterpret_cast<lapack_complex_double*>(mat.begin()), lda, ipiv.data());
  int sgn = 0;
  for (lapack_int i = 0; i < n; ++i) {
    if (ipiv[i] != i + 1) { sgn = !sgn; }
  }
  float64_t log_det = 0.0;// = (sgn ? kore::complex128_t(0.0, M_PI) : 0.0);
  for (lapack_int i = 0; i < n; ++i) {
    log_det += std::log(std::abs(mat(i, i)));
  }
  return log_det;
}


inline lapack_int
eigensystem_symmetric(Array<float64_t, 2> & matrix,
                      Array<float64_t, 1> & eigvals)
{
  lapack_int n   = static_cast<lapack_int>(eigvals.shape(0));
  lapack_int n1  = static_cast<lapack_int>(matrix.shape(0));
  lapack_int n2  = static_cast<lapack_int>(matrix.shape(1));
  lapack_int lda = static_cast<lapack_int>(matrix.stride(0));
  assert(static_cast<size_t>(n)   == eigvals.shape(0));
  assert(static_cast<size_t>(n1)  == matrix.shape(0));
  assert(static_cast<size_t>(n2)  == matrix.shape(1));
  assert(static_cast<size_t>(lda) == matrix.stride(0));
  assert(n == n1 && n == n2);

  lapack_int info = LAPACKE_dsyevd(LAPACK_ROW_MAJOR, 'V', 'L', n,
                                   matrix.begin(), lda, eigvals.begin());
  assert(info == 0);
  return info;
}


inline lapack_int
eigensystem_hermitian(Array<complex128_t, 2> & matrix,
                      Array<float64_t, 1> & eigvals)
{
  lapack_int n   = static_cast<lapack_int>(eigvals.shape(0));
  lapack_int n1  = static_cast<lapack_int>(matrix.shape(0));
  lapack_int n2  = static_cast<lapack_int>(matrix.shape(1));
  lapack_int lda = static_cast<lapack_int>(matrix.stride(0));
  assert(static_cast<size_t>(n)   == eigvals.shape(0));
  assert(static_cast<size_t>(n1)  == matrix.shape(0));
  assert(static_cast<size_t>(n2)  == matrix.shape(1));
  assert(static_cast<size_t>(lda) == matrix.stride(0));
  assert(n == n1 && n == n2);

  lapack_int info = LAPACKE_zheevd(LAPACK_ROW_MAJOR, 'V', 'L',
                                   n,
                                   reinterpret_cast<lapack_complex_double*>(matrix.begin()),
                                   lda,
                                   eigvals.begin());
  assert(info == 0);
  return info;
}
















} // namespace inplace
} // namespace linalg
} // namespace array
} // namespace kore
