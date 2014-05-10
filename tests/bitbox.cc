#include <cassert>
#include <bitset>
#include <random>
#include <iostream>
#include <boost/utility.hpp>

#include "bitbox.h"

#define BOOST_TEST_DYN_LINK
//#define BOOST_TEST_MAIN
#define BOOST_TEST_MODULE BitboxTest

#include <boost/test/unit_test.hpp>

BOOST_AUTO_TEST_CASE(simpler)
{
  using namespace std;
  using namespace kore::bitbox;
  mt19937_64 rangen;

  for (int i = 0 ; i < 1024 ; ++i) {
    uint64_t v = rangen();
    for (int digit = 0 ; digit < 64 ; ++digit) {
        uint64_t v2;
        v2 = kore::bitbox::bitset(v, digit, true);
        BOOST_CHECK(bitcheck(v2, digit));
        v2 = kore::bitbox::bitset(v, digit, false);
        BOOST_CHECK(!bitcheck(v2, digit));
    }
  }
}

BOOST_AUTO_TEST_CASE(naiveTest)
{
  using namespace std;
  using namespace kore::bitbox;

  mt19937_64 rangen;

  for (int i = 0 ; i < 1024 ; ++i) {
    uint64_t v = rangen();
    uint64_t m = rangen();
    BOOST_CHECK(bitcount(v,m)       == kore::bitbox::naive::bitcount(v,m));
    BOOST_CHECK(bitwsum(v,m)        == kore::bitbox::naive::bitwsum(v,m));
    BOOST_CHECK(bitcompress(v,m)    == kore::bitbox::naive::bitcompress(v,m));
    BOOST_CHECK(bitexpand(v,m)      == kore::bitbox::naive::bitexpand(v,m));
    BOOST_CHECK(partitioncount(v,m) == kore::bitbox::naive::partitioncount(v,m));
  }

  for (size_t f = 0 ; f < 64 ; ++f) {
    for (size_t l = 0 ; l < 64 ; ++l) {
      uint64_t m1 = makemask<uint64_t>(f, l);
      uint64_t m2 = kore::bitbox::naive::makemask<uint64_t>(f,l);
      BOOST_CHECK(m1==m2);
    }
  }
}
