#pragma once

#ifdef USE_CXX11

#ifdef USE_TR1
  #include <tr1/cstdint>
  #include <tr1/cinttypes>
  #include <tr1/type_traits>
namespace std { using namespace tr1; };
#else
  #include <cstdint>
  #include <cinttypes>
  #include <type_traits>
#endif

#include <cassert>
#include <complex>
#include <limits>
#include <stdexcept>
#include <vector>
#include <iterator>

namespace kore {

typedef std::intmax_t       intmax_t;
typedef std::int8_t         int8_t;
typedef std::int16_t        int16_t;
typedef std::int32_t        int32_t;
typedef std::int64_t        int64_t;
typedef std::int_least8_t   int_least8_t;
typedef std::int_least16_t  int_least16_t;
typedef std::int_least32_t  int_least32_t;
typedef std::int_least64_t  int_least64_t;
typedef std::int_fast8_t    int_fast8_t;
typedef std::int_fast16_t   int_fast16_t;
typedef std::int_fast32_t   int_fast32_t;
typedef std::int_fast64_t   int_fast64_t;
typedef std::intptr_t       intptr_t;

typedef std::uintmax_t      uintmax_t;
typedef std::uint8_t        uint8_t;
typedef std::uint16_t       uint16_t;
typedef std::uint32_t       uint32_t;
typedef std::uint64_t       uint64_t;
typedef std::uint_least8_t  uint_least8_t;
typedef std::uint_least16_t uint_least16_t;
typedef std::uint_least32_t uint_least32_t;
typedef std::uint_least64_t uint_least64_t;
typedef std::uint_fast8_t   uint_fast8_t;
typedef std::uint_fast16_t  uint_fast16_t;
typedef std::uint_fast32_t  uint_fast32_t;
typedef std::uint_fast64_t  uint_fast64_t;
typedef std::uintptr_t      uintptr_t;

typedef float float32_t;
typedef double float64_t;
typedef std::complex<float32_t> complex64_t;
typedef std::complex<float64_t> complex128_t;

static_assert(sizeof(float32_t)==sizeof(int32_t), "size of float32_t is not 32 bits");
static_assert(sizeof(float64_t)==sizeof(int64_t), "size of float64_t is not 64 bits");
}

#ifdef QUAD_MATH

#include <boost/multiprecision/float128.hpp>
namespace kore {
typedef boost::multiprecision::float128 float128_t;
typedef std::complex<float128_t> complex256_t;
}

#endif


// Array
namespace kore {

}

#endif
