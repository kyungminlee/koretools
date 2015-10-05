#pragma once

#include "array_linalg.h"

#ifdef USE_MKL
#include <mkl.h>
#include <mkl_lapacke.h>
#else
#include <cblas.h>
#include <lapacke.h>
#endif

#include "array_linalg_inplace.h"

namespace kore {
namespace array {
namespace linalg {

/// ----- Matrix Multiplication
/// C = A * B
Array<float64_t, 2> dot(Array<float64_t const, 2> const & a, 
                        Array<float64_t const, 2> const & b)
{
  lapack_int m   = static_cast<lapack_int>(a.shape(0));
  lapack_int k1  = static_cast<lapack_int>(a.shape(1));
  lapack_int k2  = static_cast<lapack_int>(b.shape(0));
  lapack_int n   = static_cast<lapack_int>(b.shape(1));
  lapack_int lda = static_cast<lapack_int>(a.stride(0));
  lapack_int ldb = static_cast<lapack_int>(b.stride(0));

  assert(static_cast<size_t>(m)   == a.shape(0));
  assert(static_cast<size_t>(k1)  == a.shape(1));
  assert(static_cast<size_t>(k2)  == b.shape(0));
  assert(static_cast<size_t>(n)   == b.shape(1));
  assert(static_cast<size_t>(lda) == a.stride(0));
  assert(static_cast<size_t>(ldb) == b.stride(0));

  assert(k1 == k2);
  lapack_int k = k1;

  Array<float64_t, 2> c(n, m);
  c.fill(0.0);

  lapack_int ldc = static_cast<lapack_int>(c.stride(0));
  assert(static_cast<size_t>(ldc) == c.stride(0));

  float64_t alpha = 1.0, beta = 0.0;
  cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
              m, n, k,
              alpha,
              a.cbegin(), lda,
              b.cbegin(), ldb,
              beta,
              c.begin(), ldc);

  return c;
}

/// ----- Matrix Multiplication
/// C = A * B
Array<complex128_t, 2> dot(Array<complex128_t const, 2> const & a,
                           Array<complex128_t const, 2> const & b)
{
  lapack_int m   = static_cast<lapack_int>(a.shape(0));
  lapack_int k1  = static_cast<lapack_int>(a.shape(1));
  lapack_int k2  = static_cast<lapack_int>(b.shape(0));
  lapack_int n   = static_cast<lapack_int>(b.shape(1));
  lapack_int lda = static_cast<lapack_int>(a.stride(0));
  lapack_int ldb = static_cast<lapack_int>(b.stride(0));

  assert(static_cast<size_t>(m)   == a.shape(0));
  assert(static_cast<size_t>(k1)  == a.shape(1));
  assert(static_cast<size_t>(k2)  == b.shape(0));
  assert(static_cast<size_t>(n)   == b.shape(1));
  assert(static_cast<size_t>(lda) == a.stride(0));
  assert(static_cast<size_t>(ldb) == b.stride(0));

  assert(k1 == k2);
  lapack_int k = k1;

  Array<complex128_t, 2> c(n, m);
  c.fill(0.0);

  lapack_int ldc = static_cast<lapack_int>(c.stride(0));
  assert(static_cast<size_t>(ldc) == c.stride(0));

  complex128_t alpha = 1.0, beta = 0.0;
  cblas_zgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
              m, n, k,
              reinterpret_cast<const double*>(&alpha),
              reinterpret_cast<const double*>(a.cbegin()), lda,
              reinterpret_cast<const double*>(b.cbegin()), ldb,
              reinterpret_cast<const double*>(&beta),
              reinterpret_cast<      double*>(c.begin()), ldc);
  return c;
}


template <typename ValueType> inline
Array<typename std::remove_const<ValueType>::type, 2>
transpose(Array<ValueType const, 2> const & arr)
{
  size_t n1 = arr.shape(0), n2 = arr.shape(1);
  Array<typename std::remove_const<ValueType>::type, 2> out(n2, n1);
  for (size_t i1 = 0; i1 < n1; ++i1) {
    for (size_t i2 = 0; i2 < n2; ++i2) {
      out(i2, i1) = arr(i1, i2);
    } // for i2
  } // for i1
  return out;
}


inline Array<complex128_t, 2>
inverse(Array<complex128_t const, 2> const & mat)
{
  Array<complex128_t, 2> out = mat.clone();
  inplace::inverse(out);
  return out;
}


inline Array<float64_t, 2>
inverse(Array<float64_t const, 2> const & mat)
{
  Array<float64_t, 2> out = mat.clone();
  inplace::inverse(out);
  return out;
}


//! \name READ ONLY ARRAYS
//@{
inline complex128_t
determinant(Array<complex128_t const, 2> const & mat)
{
  Array<complex128_t, 2> m2 = mat.clone();
  return inplace::determinant(m2);
}

inline float64_t
determinant(Array<float64_t const, 2> const& mat)
{
  Array<float64_t, 2> m2 = mat.clone();
  return inplace::determinant(m2);
}

inline complex128_t
log_determinant(Array<complex128_t const, 2> const & mat)
{
  Array<complex128_t, 2> m2 = mat.clone();
  return inplace::log_determinant(m2);
}

inline complex128_t
log_determinant(Array<float64_t const, 2> const & mat)
{
  Array<float64_t, 2> m2 = mat.clone();
  return inplace::log_determinant(m2);
}

inline float64_t
log_determinant_abs(Array<float64_t const, 2> const & mat)
{
  Array<float64_t, 2> m2 = mat.clone();
  return inplace::log_determinant_abs(m2);
}

inline float64_t
log_determinant_abs(Array<complex128_t const, 2> const & mat)
{
  Array<complex128_t, 2> m2 = mat.clone();
  return inplace::log_determinant_abs(m2);
}

//@}


inline std::tuple<Array<float64_t, 1>, Array<float64_t, 2> >
eigensystem_symmetric(Array<float64_t const, 2> const & matrix)
{
  std::tuple<Array<float64_t, 1>, Array<float64_t, 2> > ret
    = std::make_tuple(Array<float64_t, 1>(matrix.shape(0)), matrix.clone());
  auto & eigvals = std::get<0>(ret);
  auto & eigvecs = std::get<1>(ret);
  //Array<float64_t, 2> eigvecs = matrix.clone();
  //Array<float64_t, 1> eigvals(matrix.shape(0));
  lapack_int info = inplace::eigensystem_symmetric(eigvecs, eigvals);
  //return std::move(std::make_tuple(std::move(eigvals), std::move(eigvecs)));
  return ret;
}



inline std::tuple<Array<float64_t, 1>, Array<complex128_t, 2> >
eigensystem_hermitian(Array<complex128_t const, 2> const & matrix)
{
  std::tuple<Array<float64_t, 1>, Array<complex128_t, 2> > ret
    = std::make_tuple(Array<float64_t, 1>(matrix.shape(0)), matrix.clone());
  auto & eigvals = std::get<0>(ret);
  auto & eigvecs = std::get<1>(ret);
  lapack_int info = inplace::eigensystem_hermitian(eigvecs, eigvals);
  //return std::move(std::make_tuple(std::move(eigvals), std::move(eigvecs)));
  return ret;
}






} // namespace linalg
} // namespace array
} // namespace kore
