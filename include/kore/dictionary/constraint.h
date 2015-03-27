#pragma once

// Binary string dictionary.
//
//

#include <exception>
#include <stdexcept>
#include "../bitbox/bitbox.h"

namespace kore {
namespace dictionary {

template <typename BitString>
struct Generic
{
  Generic() { }
  inline bool operator()(BitString w) { return true; }
};

template <typename BitString>
struct Any
{
 public:
  Any() { }
  inline bool operator()(BitString w) { return true; }
};

template <typename BitString>
struct FixedDensity : public Generic<BitString>
{
 public:
  FixedDensity(size_t k) : _k(k) { 
      if(k <= 0) { throw std::domain_error("k should be positive"); }
      if(k > sizeof(BitString)) { throw std::domain_error("k should not exceed size of bitstring"); }
  }  
  inline bool operator()(BitString val) const { return (bitbox::bitcount(val) == _k); }
 private:
  size_t _k;
};

template <typename BitString>
struct FixedDensityModulus : public Generic<BitString>
{
 public:
  FixedDensityModulus(size_t k, size_t m, size_t r) : _k(k), _m(m), _r(r) {
      if(k <= 0) { throw std::domain_error("k should be positive"); }
      if(k > sizeof(BitString)*8) { throw std::domain_error("k should not exceed size of bitstring"); }
      if(m <= 0) { throw std::domain_error("m should be positive"); }
      if(r < 0 || r >= m) { throw std::domain_error("r should be in [0, m)"); }
  }

  inline bool operator()(BitString val) const {
    return (bitbox::bitcount(val) == _k) && 
           ((bitbox::bitwsum(val) % _m) == (_r % _m));
  }
 private:
  size_t _k, _m, _r;
};

} // namespace dictionary
} // namespace kore
