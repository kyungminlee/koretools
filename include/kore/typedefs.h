#pragma once

#include <cstdint>
#include <cassert>
#include <complex>
#include <limits>
#include <stdexcept>
#include <vector>
#include <iterator>

namespace kore {

#if 0
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
#endif

using intmax_t = std::intmax_t;
using int8_t = std::int8_t;
using int16_t = std::int16_t;
using int32_t = std::int32_t;
using int64_t = std::int64_t;
using int_least8_t = std::int_least8_t;
using int_least16_t = std::int_least16_t;
using int_least32_t = std::int_least32_t;
using int_least64_t = std::int_least64_t;
using int_fast8_t = std::int_fast8_t;
using int_fast16_t = std::int_fast16_t;
using int_fast32_t = std::int_fast32_t;
using int_fast64_t = std::int_fast64_t;
using intptr_t = std::intptr_t;

using uintmax_t = std::uintmax_t;
using uint8_t = std::uint8_t;
using uint16_t = std::uint16_t;
using uint32_t = std::uint32_t;
using uint64_t = std::uint64_t;
using uint_least8_t = std::uint_least8_t;
using uint_least16_t = std::uint_least16_t;
using uint_least32_t = std::uint_least32_t;
using uint_least64_t = std::uint_least64_t;
using uint_fast8_t = std::uint_fast8_t;
using uint_fast16_t = std::uint_fast16_t;
using uint_fast32_t = std::uint_fast32_t;
using uint_fast64_t = std::uint_fast64_t;
using uintptr_t = std::uintptr_t;

using float32_t = float;
using float64_t = double;
using complex64_t = std::complex<float32_t>;
using complex128_t = std::complex<float64_t>;



static_assert(sizeof(float32_t) == 4u, "size of float32_t should be 4 bytes");
static_assert(sizeof(float64_t) == 8u, "size of float64_t should be 8 bytes");
static_assert(sizeof(complex64_t) == 8u, "size of complex64_t should be 8 bytes");
static_assert(sizeof(complex128_t) == 16u, "size of complex128_t should be 16 bytes");


} // namespace kore
