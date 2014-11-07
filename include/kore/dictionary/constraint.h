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
#if 0
  inline bool operator==(const FixedDensityModulus & rhs) const {
    return (_k == rhs._k) && (_m == rhs._m) && (_r == rhs._r); 
  }
  inline bool operator<(const FixedDensityModulus &rhs) const { 
    if (_k != rhs._k) { 
      return _k < rhs._k; 
	} else if (_m != rhs._m) {
	  return _m < rhs._m; 
	} else if (_r != rhs._r) {
	  return _r < rhs._r; 
	} else {
	  return false;
	}
  }
#endif
 private:
  size_t _k, _m, _r;
};

} // namespace dictionary
} // namespace kore
