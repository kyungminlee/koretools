#pragma once
#include "../typedefs.h"

namespace kore {
namespace bitdict {
namespace constraint {
  
template <typename BitString, typename Index>
struct FixedDensity
{
  static_assert(sizeof(Index) >= sizeof(BitString), "Index should be large enough for BitString");
 public:
  FixedDensity(Index n, Index k) : _n(n), _k(k) {
      if(n <= 0) { throw std::domain_error("n should be positive"); }
      if(n > sizeof(BitString)*8) { throw std::domain_error("k should not exceed size of bitstring"); }
      if(k <= 0) { throw std::domain_error("k should be positive"); }
      if(k > n) { throw std::domain_error("k should not exceed n"); }
      if(k > sizeof(BitString)*8) { throw std::domain_error("k should not exceed size of bitstring"); }
  }
  std::tuple<bool, BitString, void*> operator()(BitString val) const {
    return std::make_tuple((!(val>>_n)) && (bitbox::bitcount(val) == _k), val, nullptr);
  }
 private:
  Index _n, _k;
};

template <typename BitString, typename Index>
struct FixedDensityModulus
{
  static_assert(sizeof(Index) >= sizeof(BitString), "Index should be large enough for BitString");
 public:
  FixedDensityModulus(Index n, Index k, Index m, Index r) : _n(n), _k(k), _m(m), _r(r) {
      if(n <= 0) { throw std::domain_error("n should be positive"); }
      if(n > sizeof(BitString)*8) { throw std::domain_error("k should not exceed size of bitstring"); }
      if(k <= 0) { throw std::domain_error("k should be positive"); }
      if(k > n) { throw std::domain_error("k should not exceed n"); }
      if(k > sizeof(BitString)*8) { throw std::domain_error("k should not exceed size of bitstring"); }
      if(m <= 0) { throw std::domain_error("m should be positive"); }
      if(m > n) { throw std::domain_error("m should not exceed n"); }
      if(r < 0 || r >= m) { throw std::domain_error("r should be in [0, m)"); }
  }

  std::tuple<bool, BitString, void*> operator()(BitString val) const {
    return std::make_tuple( (!(val>>_n))
			   && (bitbox::bitcount(val) == _k) && 
			   ((bitbox::bitwsum(val) % _m) == (_r % _m)),
			   val, nullptr);
  }
 private:
  Index _n, _k, _m, _r;
};



} // namespace constraint 
} // namespace bitdict
} // namespace kore
