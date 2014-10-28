#include <cassert>
#include <bitset>
#include <random>
#include <iostream>
#include <boost/utility.hpp>

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE BitboxTest

#include <boost/test/unit_test.hpp>

#include <tuple>
#include "kore/bitbox/bitbox.h"
#include "kore/bitdict/bitdict.h"

BOOST_AUTO_TEST_CASE(any)
{
  using namespace kore;
  using namespace kore::bitdict;
  int64_t length = 8;
  BitDictionary<> sdict(length, [length](uint64_t word) -> std::tuple<bool, uint64_t, void*> {
    if ((word >> length) != 0) {
      return std::make_tuple(false, word, nullptr);
    }
    else {
      return std::make_tuple(true, word, nullptr);
    }
  });
  BOOST_CHECK_EQUAL(sdict.size(), (1 << length));

  for (uint64_t wrd = 0 ; wrd < (1<<length)+10 ; ++wrd) {
    int64_t idx = sdict.index(wrd);
    if (wrd >= (1<<length)) {
      BOOST_CHECK_EQUAL(idx, -1);
    } else {
      BOOST_CHECK_NE(idx, -1);
      uint64_t new_wrd = sdict.word(idx);
      BOOST_CHECK_EQUAL(new_wrd, wrd);
    }
  }
  for (int64_t idx = 0; idx < sdict.size(); ++idx) {
    uint64_t wrd = sdict.word(idx);
    int64_t new_idx = sdict.index(wrd);
    BOOST_CHECK_EQUAL(idx, new_idx);
  }
}

BOOST_AUTO_TEST_CASE(fixed_density)
{
  using namespace kore;
  using namespace kore::bitdict;
  int64_t length = 8;
  int64_t density = 5;
  BitDictionary<> sdict(length, [length, density](uint64_t word) -> std::tuple<bool, uint64_t, void*> {
    if ((word >> length) != 0) {
      return std::make_tuple(false, word, nullptr);
    }
    if (kore::bitbox::bitcount(word) != density) {
      return std::make_tuple(false, word, nullptr);
    } 
    return std::make_tuple(true, word, nullptr);
  });
  BOOST_CHECK_LE(sdict.size(), (1 << length));

  for (uint64_t wrd = 0 ; wrd < (1<<length)+10 ; ++wrd) {
    int64_t idx = sdict.index(wrd);
    //std::cout << "word " << std::bitset<12>(wrd) << " : index " << idx << std::endl;
    if (idx == -1) {
    } else {
      uint64_t new_wrd = sdict.word(idx);
      BOOST_CHECK_EQUAL(new_wrd, wrd);
    }
  }
  for (int64_t idx = 0; idx < sdict.size(); ++idx) {
    uint64_t wrd = sdict.word(idx);
    int64_t new_idx = sdict.index(wrd);
    BOOST_CHECK_EQUAL(idx, new_idx);
  }
}



