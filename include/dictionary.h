#pragma once
#include <exception>
#ifdef USE_TR1
  #include <tr1/unordered_map>
#else
  #include <unordered_map>
#endif
#include "bitbox.h"

namespace kore {
namespace dictionary {

template <typename BitString>
struct Generic
{
  Generic() { }
  inline bool operator()(BitString w) { return true; }
  inline bool operator==(const Generic & rhs) const { return true; }
  inline bool operator<(const Generic &rhs) const { return false; }
  inline bool operator>(const Generic &rhs) const { return false; }
};

template <typename BitString>
struct Any
{
 public:
  Any() { }
  inline bool operator()(BitString w) { return true; }
  inline bool operator==(const Any & rhs) const { return true; }
  inline bool operator<(const Any &rhs) const { return false; }
  inline bool operator>(const Any &rhs) const { return false; }
};

template <typename BitString>
struct FixedDensity : public Generic<BitString>
{
 public:
  FixedDensity(size_t k) : _k(k) { 
      if(k <= 0) { throw std::invalid_argument("k should be positive"); }
      if(k > sizeof(BitString)) { throw std::invalid_argument("k should not exceed size of bitstring"); }
  }  
  inline bool operator()(BitString val) const { return (bitbox::bitcount(val) == _k); }

  inline bool operator==(const FixedDensity & rhs) const { return (_k == rhs._k); }
  inline bool operator<(const FixedDensity &rhs) const { return _k < rhs._k; }
  inline bool operator>(const FixedDensity &rhs) const { return _k > rhs._k; }
 private:
  size_t _k;
};

template <typename BitString>
struct FixedDensityModulus : public Generic<BitString>
{
 public:
  FixedDensityModulus(size_t k, size_t m, size_t r) : _k(k), _m(m), _r(r) {
      if(k <= 0) { throw std::invalid_argument("k should be positive"); }
      if(k > sizeof(BitString)*8) { throw std::invalid_argument("k should not exceed size of bitstring"); }
      if(m <= 0) { throw std::invalid_argument("m should be positive"); }
      if(r < 0 || r >= m) { throw std::invalid_argument("r should be in [0, m)"); }
  }

  inline bool operator()(BitString val) const {
    return (bitbox::bitcount(val) == _k) && 
           ((bitbox::bitwsum(val) % _m) == (_r % _m));
  }

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
 private:
  size_t _k, _m, _r;
};

template <typename BitString, typename IndexType=int64_t>
struct Dictionary
{
  static_assert(sizeof(BitString) <= sizeof(IndexType), "IndexType should be large enough to store BitString");
 public:
  template <typename CriterionType>
  Dictionary(size_t n, CriterionType criterion): _n(n) //, _indices( (0x1L<<n), int64_t(-1L) )
  {
    BitString ONE = 0x1;
    if(n <= 0) { throw std::invalid_argument("n should be positive"); }
    if(n > sizeof(BitString)*8) { throw std::invalid_argument("n should not exceed size of bitstring"); }

    BitString val = 0;
    while( !((ONE << n) & val ) && !criterion(val) ) { val++; }
    // find first that matches criterion
    
    IndexType cnt = 0;
    while( !((ONE << n) & val )) {
      _indices[val] = cnt;
      _words.push_back(val);
      cnt++;
      do {
        val++;
      } while( !((ONE << n) & val ) && !criterion(val) );
    } // while length shorter than n
    if (cnt == 0) { throw std::invalid_argument("cnt should be at least 1"); }
  }

  size_t n()      const { return _n; }
  size_t length() const { return _n; }
  size_t size()   const { return _words.size(); }
#ifdef __clang__
  IndexType index(BitString w) const { return std::get<1>(*_indices.find(w)); }
#else
  IndexType index(BitString w) const { return _indices.at(w); }
#endif
  BitString word(size_t idx) const { return _words[idx]; }
  
  typename std::vector<BitString>::const_iterator begin()  const { return _words.begin(); }
  typename std::vector<BitString>::const_iterator end()    const { return _words.end(); }
  typename std::vector<BitString>::const_iterator cbegin() const { return _words.cbegin(); }
  typename std::vector<BitString>::const_iterator cend()   const { return _words.cend(); }

 private:
  size_t _n;
  //std::vector<IndexType> _indices;
  std::unordered_map<BitString, IndexType> _indices;
  std::vector<BitString> _words;
};

} // namespace dictionary
} // namespace kore
