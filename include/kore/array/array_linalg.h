#pragma once
#include <complex>
#include <tuple>
#include "array.h"

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

namespace internal {
  template <typename T1, typename TgtType>
  struct true_if_match_1
  {
    static const bool value =
      std::is_same<typename std::remove_const<T1>::type, TgtType>::value;
  };

  template <typename T1, typename T2, typename TgtType>
  struct true_if_match_2
  {
    static const bool value =
      std::is_same<typename std::remove_const<T1>::type, TgtType>::value &&
      std::is_same<typename std::remove_const<T2>::type, TgtType>::value;
  };

  template <typename T1, typename TgtType, typename RetType = void>
  struct type_if_match_1
  {
    using type = typename std::enable_if<
      std::is_same<typename std::remove_const<T1>::type, TgtType>::value,
      RetType>::type;
  };

  template <typename T1, typename T2, typename TgtType, typename RetType = void>
  struct type_if_match_2
  {
    using type = typename std::enable_if<
      std::is_same<typename std::remove_const<T1>::type, TgtType>::value &&
      std::is_same<typename std::remove_const<T2>::type, TgtType>::value,
      RetType>::type;
  };

  template <typename T1, typename TgtType, typename RetType>
  using type_if_match_1_t = typename type_if_match_1<T1, TgtType, RetType>::type;

  template <typename T1, typename T2, typename TgtType, typename RetType>
  using type_if_match_2_t = typename type_if_match_2<T1, T2, TgtType, RetType>::type;
} // namespace internal


//! Matrix Multiplication
//! \f$ C = A \cdot B \f$
template <typename T1, typename T2,
  typename = typename std::enable_if<
        internal::true_if_match_2<T1, T2, float64_t>::value>::type >
  inline Array<float64_t, 2>
  dot(Array<T1, 2> const & a,
      Array<T2, 2> const & b)
{
  lapack_int m = static_cast<lapack_int>(a.shape(0));
  lapack_int k1 = static_cast<lapack_int>(a.shape(1));
  lapack_int k2 = static_cast<lapack_int>(b.shape(0));
  lapack_int n = static_cast<lapack_int>(b.shape(1));
  lapack_int lda = static_cast<lapack_int>(a.stride(0));
  lapack_int ldb = static_cast<lapack_int>(b.stride(0));

  {
    KORE_ASSERT(static_cast<size_t>(m) == a.shape(0), "overflow");
    KORE_ASSERT(static_cast<size_t>(k1) == a.shape(1), "overflow");
    KORE_ASSERT(static_cast<size_t>(k2) == b.shape(0), "overflow");
    KORE_ASSERT(static_cast<size_t>(n) == b.shape(1), "overflow");
    KORE_ASSERT(static_cast<size_t>(lda) == a.stride(0), "overflow");
    KORE_ASSERT(static_cast<size_t>(ldb) == b.stride(0), "overflow");
    KORE_ASSERT(k1 == k2, "dimensions do not match");
  }
  lapack_int k = k1;

  Array<float64_t, 2> c(n, m);
  c.fill(0.0);

  lapack_int ldc = static_cast<lapack_int>(c.stride(0));
  KORE_ASSERT(static_cast<size_t>(ldc) == c.stride(0), "overflow");

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


//! Matrix - Matrix Multiplication
//! \f$ C = A \cdot B \f$
template <typename T1, typename T2,
  typename = typename std::enable_if<
  internal::true_if_match_2<T1, T2, complex128_t>::value >::type>
  inline Array<complex128_t, 2>
  dot(Array<T1, 2> const & a,
      Array<T2, 2> const & b)
{
  lapack_int m = static_cast<lapack_int>(a.shape(0));
  lapack_int k1 = static_cast<lapack_int>(a.shape(1));
  lapack_int k2 = static_cast<lapack_int>(b.shape(0));
  lapack_int n = static_cast<lapack_int>(b.shape(1));
  lapack_int lda = static_cast<lapack_int>(a.stride(0));
  lapack_int ldb = static_cast<lapack_int>(b.stride(0));

  { //ASSERTIONS
    KORE_ASSERT(static_cast<size_t>(m) == a.shape(0), "Overflow");
    KORE_ASSERT(static_cast<size_t>(k1) == a.shape(1), "Overflow");
    KORE_ASSERT(static_cast<size_t>(k2) == b.shape(0), "Overflow");
    KORE_ASSERT(static_cast<size_t>(n) == b.shape(1), "Overflow");
    KORE_ASSERT(static_cast<size_t>(lda) == a.stride(0), "Overflow");
    KORE_ASSERT(static_cast<size_t>(ldb) == b.stride(0), "Overflow");
    KORE_ASSERT(k1 == k2, "dimensions do not match");
  }
  lapack_int k = k1;

  Array<complex128_t, 2> c(n, m);
  c.fill(0.0);

  lapack_int ldc = static_cast<lapack_int>(c.stride(0));
  KORE_ASSERT(static_cast<size_t>(ldc) == c.stride(0), "Overflow");

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


template <typename T1, typename T2,
  typename = typename std::enable_if<
  !internal::true_if_match_2<T1, T2, float64_t>::value &&
  !internal::true_if_match_2<T1, T2, complex128_t>::value
>::type >
auto dot(Array<T1, 2> const & a,
         Array<T2, 2> const & b)
  ->Array<typename std::remove_const<typename std::remove_reference<
  decltype(*a.cbegin() * *b.cbegin() + *a.cbegin() * *b.cbegin())
  >::type>::type , 2>
{
  using ValueType =
    typename std::remove_const<typename std::remove_reference<
    decltype(*a.cbegin() * *b.cbegin() + *a.cbegin() * *b.cbegin())
    >::type>::type ;

  size_t m = a.shape(0);
  size_t k1 = a.shape(1);
  size_t k2 = b.shape(0);
  size_t n = b.shape(1);
  KORE_ASSERT(k1 == k2, "dimensions do not match");

  Array<ValueType, 2> foo(m, n);
  foo.fill(ValueType(0));
  for (size_t i = 0; i < m; ++i) {
    for (size_t j = 0; j < n; ++j) {
      for (size_t k = 0; k < k1; ++k) {
        foo(i, j) += a(i, k) * b(k, j);
      }
    }
  }
  return foo;
}

//@}


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

template <typename T>
inline typename std::enable_if<internal::true_if_match_1<T, complex128_t>::value, Array<complex128_t, 2>>::type
inverse(Array<T, 2> const & mat)
{
  Array<complex128_t, 2> out = mat.clone();
  inplace::inverse(out);
  return out;
}

template <typename T>
inline typename std::enable_if<internal::true_if_match_1<T, float64_t>::value, Array<float64_t, 2>>::type
inverse(Array<T, 2> const & mat)
{
  Array<float64_t, 2> out = mat.clone();
  inplace::inverse(out);
  return out;
}

//! \name determinant
//@{
template <typename T>
inline typename std::enable_if<internal::true_if_match_1<T, complex128_t>::value, complex128_t>::type
determinant(Array<T, 2> const & mat)
{
  Array<complex128_t, 2> m2 = mat.clone();
  return inplace::determinant(m2);
}

template <typename T>
inline typename std::enable_if<internal::true_if_match_1<T, float64_t>::value, float64_t>::type
determinant(Array<T, 2> const& mat)
{
  Array<float64_t, 2> m2 = mat.clone();
  return inplace::determinant(m2);
}

//@}

//! \name log_determinant
//@{

template <typename T>
inline typename std::enable_if<internal::true_if_match_1<T, float64_t>::value, complex128_t>::type
log_determinant(Array<T, 2> const & mat)
{
  Array<float64_t, 2> m2 = mat.clone();
  return inplace::log_determinant(m2);
}

template <typename T>
inline typename std::enable_if<internal::true_if_match_1<T, complex128_t>::value, complex128_t>::type
log_determinant(Array<T, 2> const & mat)
{
  Array<complex128_t, 2> m2 = mat.clone();
  return inplace::log_determinant(m2);
}

//@}

//! \name log_determinant_abs
//@{

template <typename T>
inline typename std::enable_if<internal::true_if_match_1<T, float64_t>::value, float64_t>::type
log_determinant_abs(Array<T, 2> const & mat)
{
  Array<float64_t, 2> m2 = mat.clone();
  return inplace::log_determinant_abs(m2);
}

template <typename T>
inline typename std::enable_if<internal::true_if_match_1<T, complex128_t>::value, float64_t>::type
log_determinant_abs(Array<T, 2> const & mat)
{
  Array<complex128_t, 2> m2 = mat.clone();
  return inplace::log_determinant_abs(m2);
}


template <typename T>
inline typename std::enable_if<internal::true_if_match_1<T, complex128_t>::value, float64_t>::type
log_determinant_abs_squared(Array<T, 2> const & mat)
{
  Array<complex128_t, 2> m2 = mat.clone();
  return inplace::log_determinant_abs_squared(m2);
}


//@}

//! \name eigensystem_symmetric
//@{

template <typename T>
inline typename std::enable_if<internal::true_if_match_1<T, float64_t>::value,
  std::tuple<Array<float64_t, 1>, Array<float64_t, 2>>
>::type
eigensystem_symmetric(Array<T, 2> const & matrix)
{
  std::tuple<Array<float64_t, 1>, Array<float64_t, 2> >
    ret(Array<float64_t, 1>(matrix.shape(0)), matrix.clone());

  auto & eigvals = std::get<0>(ret);
  auto & eigvecs = std::get<1>(ret);
  lapack_int info = inplace::eigensystem_symmetric(eigvecs, eigvals);
  return ret;
}

template <typename T>
inline typename std::enable_if<internal::true_if_match_1<T, float64_t>::value,
  std::tuple<Array<float64_t, 1>, Array<float64_t, 2>>
>::type
eigensystem_hermitian(Array<T, 2> const & matrix)
{
  std::tuple<Array<float64_t, 1>, Array<float64_t, 2> > 
    ret(Array<float64_t, 1>(matrix.shape(0)), matrix.clone());
  auto & eigvals = std::get<0>(ret);
  auto & eigvecs = std::get<1>(ret);
  lapack_int info = inplace::eigensystem_symmetric(eigvecs, eigvals);
  return ret;
}

template <typename T>
inline typename std::enable_if<internal::type_if_match_1<T, complex128_t>::value,
  std::tuple<Array<float64_t, 1>, Array<complex128_t, 2>>
>::type
eigensystem_hermitian(Array<T, 2> const & matrix)
{
  std::tuple<Array<float64_t, 1>, Array<complex128_t, 2> > 
    ret(Array<float64_t, 1>(matrix.shape(0)), matrix.clone());
  auto & eigvals = std::get<0>(ret);
  auto & eigvecs = std::get<1>(ret);
  lapack_int info = inplace::eigensystem_hermitian(eigvecs, eigvals);
  return ret;
}

} // namespace linalg
} // namespace array
} // namespace kore

#include "array_linalg_impl.h"
