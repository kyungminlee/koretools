#pragma once

#include "array_linalg.h"

#ifdef USE_MKL
#include <mkl.h>
#include <mkl_lapacke.h>
#else
#include <lapacke.h>
#endif

namespace kore {

namespace array {

namespace linalg {

/// ----- Matrix Multiplication
/// C = A * B
inline
Array<float64_t, 2> dot(const Array<float64_t, 2>& a, const Array<float64_t, 2>& b)
{
  auto & fortran_a = b;
  auto & fortran_b = a;

  lapack_int m = static_cast<lapack_int>(a.shape(0));
  lapack_int k1 = static_cast<lapack_int>(a.shape(1));
  lapack_int k2 = static_cast<lapack_int>(b.shape(0));
  lapack_int n = static_cast<lapack_int>(b.shape(1));
  lapack_int lda = static_cast<lapack_int>(a.stride(0));
  lapack_int ldb = static_cast<lapack_int>(b.stride(0));

  assert(m == a.shape(0));
  assert(k1 == a.shape(1));
  assert(k2 == b.shape(0));
  assert(n == b.shape(1));
  assert(lda == a.stride(0));
  assert(ldb == b.stride(0));

  assert(k1 == k2);
  lapack_int k = k1;

  Array<float64_t, 2> c(n, m);
  c.fill(0.0);

  lapack_int ldc = static_cast<lapack_int>(c.stride(0));
  assert(ldc == c.stride(0));

  float64_t alpha = 1.0;
  cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, 
    m, n, k, alpha, a.begin(), lda, b.begin(), ldb, 0.0, c.begin(), ldc);

  return c;
}


template <typename ValueType>
Array<ValueType, 2> transpose(const Array<ValueType, 2>& arr)
{
  size_t n1 = arr.shape(0), n2 = arr.shape(1);

  Array<ValueType, 2> out(n2, n1);
  for (size_t i1 = 0; i1 < n1; ++i1) {
    for (size_t i2 = 0; i2 < n2; ++i2) {
      out(i2, i1) = arr(i1, i2);
    } // for i2
  } // for i1
  return out;
}

template <typename ValueType>
Array<ValueType, 2> transpose(const ConstArray<ValueType, 2>& arr)
{
  size_t n1 = arr.shape(0), n2 = arr.shape(1);

  Array<ValueType, 2> out(n2, n1);
  for (size_t i1 = 0; i1 < n1; ++i1) {
    for (size_t i2 = 0; i2 < n2; ++i2) {
      out(i2, i1) = arr(i1, i2);
    } // for i2
  } // for i1
  return out;
}

inline
Array<kore::complex128_t, 2> inverse(const Array<kore::complex128_t, 2>& mat)
{
  lapack_int m = static_cast<lapack_int>(mat.shape(0));
  lapack_int n = static_cast<lapack_int>(mat.shape(1));
  lapack_int lda = static_cast<lapack_int>(mat.stride(0));
  assert(m == mat.shape(0));
  assert(n == mat.shape(1));
  assert(lda == mat.stride(0));
  assert(m == n);

  Array<kore::complex128_t, 2> out(mat.clone());
  std::vector<lapack_int> ipiv(n);

  LAPACKE_zgetrf(LAPACK_ROW_MAJOR, m, n,
    reinterpret_cast<lapack_complex_float64_t*>(out.begin()), lda, ipiv.data());
  LAPACKE_zgetri(LAPACK_ROW_MAJOR, n, reinterpret_cast<lapack_complex_float64_t*>(out.begin()), lda, ipiv.data());
  return out;
}

inline
Array<float64_t, 2> inverse(const Array<float64_t, 2>& mat)
{
  lapack_int m = static_cast<lapack_int>(mat.shape(0));
  lapack_int n = static_cast<lapack_int>(mat.shape(1));
  lapack_int lda = static_cast<lapack_int>(mat.stride(0));
  assert(m == mat.shape(0));
  assert(n == mat.shape(1));
  assert(lda == mat.stride(0));
  assert(m == n);

  Array<float64_t, 2> out(mat.clone());
  std::vector<lapack_int> ipiv(n);

  LAPACKE_dgetrf(LAPACK_ROW_MAJOR, m, n, out.begin(), lda, ipiv.data());
  LAPACKE_dgetri(LAPACK_ROW_MAJOR, n, out.begin(), lda, ipiv.data());
  return out;
}




std::tuple<Array<float64_t, 1>, Array<float64_t, 2> >
eigensystem_hermitian(const Array<float64_t, 2>& matrix)
{
#ifndef NDEBUG
  if (matrix.shape(0) != matrix.shape(1)) {
    throw std::length_error("kore::linalg::eigh() matrix not square");
  }
#endif
  lapack_int n = static_cast<lapack_int>(matrix.shape(0));
  lapack_int lda = static_cast<lapack_int>(matrix.stride(0));
  assert(n == matrix.shape(0));
  assert(lda == matrix.stride(0));

  Array<float64_t, 1> eigvals(n);
  Array<float64_t, 2> eigvecs = matrix.clone();
  {
    lapack_int info = LAPACKE_dsyevd(LAPACK_ROW_MAJOR, 'V', 'L', n, eigvecs.begin(), lda, eigvals.begin());
    assert(info == 0);
  }
  return std::make_tuple(std::move(eigvals), std::move(eigvecs));
}

std::tuple<Array<float64_t, 1>, Array<kore::complex128_t, 2> >
eigensystem_hermitian(const Array<kore::complex128_t, 2>& matrix)
{
#ifndef NDEBUG
  if (matrix.shape(0) != matrix.shape(1)) {
    throw std::length_error("kore::linalg::eigh() matrix not square");
  }
#endif
  lapack_int n = static_cast<lapack_int>(matrix.shape(0));
  lapack_int lda = static_cast<lapack_int>(matrix.stride(0));
  assert(n == matrix.shape(0));
  assert(lda == matrix.stride(0));

  Array<float64_t, 1> eigvals(n);
  Array<kore::complex128_t, 2> eigvecs = matrix.clone();
  {
    lapack_int info = LAPACKE_zheevd(LAPACK_ROW_MAJOR, 'V', 'L',
      n, reinterpret_cast<lapack_complex_float64_t*>(eigvecs.begin()),
      lda, eigvals.begin());
    assert(info == 0);
  }
  return std::make_tuple(std::move(eigvals), std::move(eigvecs));
}


std::tuple<Array<float64_t, 1>, Array<kore::complex128_t, 2> >
eigensystem_hermitian(const ConstArray<kore::complex128_t, 2>& matrix)
{
#ifndef NDEBUG
  if (matrix.shape(0) != matrix.shape(1)) {
    throw std::length_error("kore::linalg::eigh() matrix not square");
  }
#endif
  lapack_int n = static_cast<lapack_int>(matrix.shape(0));
  lapack_int lda = static_cast<lapack_int>(matrix.stride(0));
  assert(n == matrix.shape(0));
  assert(lda == matrix.stride(0));

  Array<float64_t, 1> eigvals(n);
  Array<kore::complex128_t, 2> eigvecs(matrix);
  {
    lapack_int info = LAPACKE_zheevd(LAPACK_ROW_MAJOR, 'V', 'L',
      n, reinterpret_cast<lapack_complex_float64_t*>(eigvecs.begin()),
      lda, eigvals.begin());
    assert(info == 0);
  }
  return std::make_tuple(std::move(eigvals), std::move(eigvecs));
}


//! Compute determinant.
//! Destroy the contents inside.
inline
kore::complex128_t determinant(Array<kore::complex128_t, 2>& mat)
{
  lapack_int m = static_cast<lapack_int>(mat.shape(0));
  lapack_int n = static_cast<lapack_int>(mat.shape(1));
  lapack_int lda = static_cast<lapack_int>(mat.stride(0));
  assert(m == mat.shape(0));
  assert(n == mat.shape(1));
  assert(lda == mat.stride(0));
  assert(m == n);

  std::vector<lapack_int> ipiv(n);
  LAPACKE_zgetrf(LAPACK_ROW_MAJOR, m, n,
    reinterpret_cast<lapack_complex_float64_t*>(mat.begin()), lda, ipiv.data());
  int sgn = 0;
  for (size_t i = 0; i < n; ++i) {
    if (ipiv[i] != i + 1) { sgn = !sgn; }
  }
  kore::complex128_t det = (sgn ? -1.0 : 1.0);
  for (size_t i = 0; i < n; ++i) {
    det *= mat(i, i);
  }
  return det;
}

inline
kore::complex128_t log_determinant(Array<kore::complex128_t, 2>& mat)
{
  lapack_int m = static_cast<lapack_int>(mat.shape(0));
  lapack_int n = static_cast<lapack_int>(mat.shape(1));
  lapack_int lda = static_cast<lapack_int>(mat.stride(0));
  assert(m == mat.shape(0));
  assert(n == mat.shape(1));
  assert(lda == mat.stride(0));
  assert(m == n);

  std::vector<lapack_int> ipiv(n);
  LAPACKE_zgetrf(LAPACK_ROW_MAJOR, m, n,
    reinterpret_cast<lapack_complex_float64_t*>(mat.begin()), lda, ipiv.data());
  int sgn = 0;
  for (size_t i = 0; i < n; ++i) {
    if (ipiv[i] != i + 1) { sgn = !sgn; }
  }
  kore::complex128_t log_det = (sgn ? kore::complex128_t(0.0, M_PI) : 0.0);
  for (size_t i = 0; i < n; ++i) {
    log_det += std::log(mat(i, i));
  }
  return log_det;
}


//! Compute determinant.
//! Destroy the contents inside.
inline
float64_t determinant(Array<float64_t, 2>& mat)
{
  lapack_int m = static_cast<lapack_int>(mat.shape(0));
  lapack_int n = static_cast<lapack_int>(mat.shape(1));
  lapack_int lda = static_cast<lapack_int>(mat.stride(0));
  assert(m == mat.shape(0));
  assert(n == mat.shape(1));
  assert(lda == mat.stride(0));
  assert(m == n);

  std::vector<lapack_int> ipiv(n);
  LAPACKE_dgetrf(LAPACK_ROW_MAJOR, m, n, mat.begin(), lda, ipiv.data());
  auto sgn = false;
  for (size_t i = 0; i < n; ++i) {
    if (ipiv[i] != i + 1) { sgn = !sgn; }
  }
  float64_t det = (sgn ? -1.0 : 1.0);
  for (size_t i = 0; i < n; ++i) {
    det *= mat(i, i);
  }
  return det;
}

//! Compute determinant.
//! Keep the contents inside.
inline
kore::complex128_t determinant(const Array<kore::complex128_t, 2>& mat)
{
  lapack_int m = static_cast<lapack_int>(mat.shape(0));
  lapack_int n = static_cast<lapack_int>(mat.shape(1));
  lapack_int lda = static_cast<lapack_int>(mat.stride(0));
  assert(m == mat.shape(0));
  assert(n == mat.shape(1));
  assert(lda == mat.stride(0));
  assert(m == n);

  std::vector<kore::complex128_t> cache_mat(n*n);
  std::copy(mat.begin(), mat.begin() + n*n, cache_mat.begin());
  std::vector<lapack_int> ipiv(n);

  LAPACKE_zgetrf(LAPACK_ROW_MAJOR, m, n,
    reinterpret_cast<lapack_complex_float64_t*>(cache_mat.data()), lda, ipiv.data());
  auto sgn = false;
  for (size_t i = 0; i < n; ++i) {
    if (ipiv[i] != i + 1) { sgn = !sgn; }
  }
  kore::complex128_t det = (sgn ? -1.0 : 1.0);
  for (size_t i = 0; i < n; ++i) {
    det *= cache_mat[i*(n + 1)];
  }
  return det;
}

//! Compute determinant.
//! Keep the contents inside.
inline
kore::complex128_t determinant(const ConstArray<kore::complex128_t, 2>& mat)
{
  lapack_int m = static_cast<lapack_int>(mat.shape(0));
  lapack_int n = static_cast<lapack_int>(mat.shape(1));
  lapack_int lda = static_cast<lapack_int>(mat.stride(0));
  assert(m == mat.shape(0));
  assert(n == mat.shape(1));
  assert(lda == mat.stride(0));
  assert(m == n);
  std::vector<kore::complex128_t> cache_mat(n*n);
  std::copy(mat.begin(), mat.begin() + n*n, cache_mat.begin());
  std::vector<lapack_int> ipiv(n);

  LAPACKE_zgetrf(LAPACK_ROW_MAJOR, m, n,
    reinterpret_cast<lapack_complex_float64_t*>(cache_mat.data()), lda, ipiv.data());
  int sgn = 0;
  for (size_t i = 0; i < n; ++i) {
    if (ipiv[i] != i + 1) { sgn = !sgn; }
  }
  kore::complex128_t det = (sgn ? -1.0 : 1.0);
  for (size_t i = 0; i < n; ++i) {
    det *= cache_mat[i*(n + 1)];
  }
  return det;
}




//! Compute determinant.
//! Keep the contents inside.
inline
float64_t determinant(const Array<float64_t, 2>& mat)
{
  lapack_int m = static_cast<lapack_int>(mat.shape(0));
  lapack_int n = static_cast<lapack_int>(mat.shape(1));
  assert(m == mat.shape(0));
  assert(n == mat.shape(1));
  assert(m == n);

  std::vector<float64_t> cache_mat(n*n);
  std::copy(mat.begin(), mat.begin() + n*n, cache_mat.begin());
  std::vector<lapack_int> ipiv(n);

  LAPACKE_dgetrf(LAPACK_ROW_MAJOR, n, n,
    cache_mat.data(), n, ipiv.data());
  int sgn = 0;
  for (size_t i = 0; i < n; ++i) {
    if (ipiv[i] != i + 1) { sgn = !sgn; }
  }
  float64_t det = (sgn ? -1.0 : 1.0);
  for (size_t i = 0; i < n; ++i) {
    det *= cache_mat[i*(n + 1)];
  }
  return det;
}


//! Compute determinant.
//! Keep the contents inside.
inline
float64_t determinant(const ConstArray<float64_t, 2>& mat)
{
  lapack_int m = static_cast<lapack_int>(mat.shape(0));
  lapack_int n = static_cast<lapack_int>(mat.shape(1));
  assert(m == mat.shape(0));
  assert(n == mat.shape(1));
  assert(m == n);

  std::vector<float64_t> cache_mat(n*n);
  std::copy(mat.begin(), mat.begin() + n*n, cache_mat.begin());
  std::vector<lapack_int> ipiv(n);

  LAPACKE_dgetrf(LAPACK_ROW_MAJOR, n, n,
    cache_mat.data(), n, ipiv.data());
  int sgn = 0;
  for (size_t i = 0; i < n; ++i) {
    if (ipiv[i] != i + 1) { sgn = !sgn; }
  }
  float64_t det = (sgn ? -1.0 : 1.0);
  for (size_t i = 0; i < n; ++i) {
    det *= cache_mat[i*(n + 1)];
  }
  return det;
}

} // namespace linalg

} // namespace array

} // namespace kore
