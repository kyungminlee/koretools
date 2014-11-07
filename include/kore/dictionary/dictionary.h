#pragma once

// Binary string dictionary.
//
//

#include <stdexcept>
#include <exception>
#include <map>
#include "../bitbox/bitbox.h"

namespace kore {
namespace dictionary {

template <typename BitString, typename IndexType=int64_t>
struct Dictionary
{
  //static_assert(sizeof(BitString) <= sizeof(IndexType), "IndexType should be large enough to store BitString");
 public:
  template <typename CriterionType>
  Dictionary(size_t n, CriterionType criterion): _n(n) //, _indices( (0x1L<<n), int64_t(-1L) )
  {
    BitString ONE = 0x1;
    if(n <= 0) { throw std::domain_error("n should be positive"); }
    if(n > sizeof(BitString)*8) { throw std::domain_error("n should not exceed size of bitstring"); }

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
    if (cnt == 0) { throw std::domain_error("cnt should be at least 1"); }
  }

  size_t n()      const { return _n; }
  size_t length() const { return _n; }
  size_t size()   const { return _words.size(); }
  IndexType index(BitString w) const { return _indices.at(w); }
  BitString word(size_t idx) const { return _words[idx]; }
  
  typename std::vector<BitString>::const_iterator begin()  const { return _words.begin(); }
  typename std::vector<BitString>::const_iterator end()    const { return _words.end(); }
  typename std::vector<BitString>::const_iterator cbegin() const { return _words.cbegin(); }
  typename std::vector<BitString>::const_iterator cend()   const { return _words.cend(); }

 private:
  size_t _n;
  std::map<BitString, IndexType> _indices;
  std::vector<BitString> _words;
};

} // namespace dictionary
} // namespace kore
