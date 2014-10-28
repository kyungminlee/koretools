#include <cassert>
#include <bitset>
#include <random>
#include <iostream>
#include <boost/utility.hpp>

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE BitboxTest

#include <boost/test/unit_test.hpp>


#include "kore/bitbox/bitbox.h"

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
    BOOST_CHECK_EQUAL(bitcount(v,m)      , kore::bitbox::naive::bitcount(v,m)      );
    BOOST_CHECK_EQUAL(bitwsum(v,m)       , kore::bitbox::naive::bitwsum(v,m)       );
    BOOST_CHECK_EQUAL(bitcompress(v,m)   , kore::bitbox::naive::bitcompress(v,m)   );
    BOOST_CHECK_EQUAL(bitexpand(v,m)     , kore::bitbox::naive::bitexpand(v,m)     );
    BOOST_CHECK_EQUAL(partitioncount(v,m), kore::bitbox::naive::partitioncount(v,m));
  }

  for (size_t f = 0 ; f < 64 ; ++f) {
    for (size_t l = 0 ; l < 64 ; ++l) {
      BOOST_CHECK_EQUAL(makemask<uint64_t>(f, l),
              kore::bitbox::naive::makemask<uint64_t>(f,l));
    }
  }
}
