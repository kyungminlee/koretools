#include <cassert>
#include <random>
#include <boost/utility.hpp>

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE CollectorTest

#include <boost/test/unit_test.hpp>

#include "kore/collector/collector.h"

BOOST_AUTO_TEST_CASE(collector_test_double)
{
  using namespace std;
  using namespace kore;
  using namespace kore::collector;
  mt19937_64 rangen;
  uniform_real_distribution<double> distrib(0.0, 1.0);

  int64_t nr = 256;
  int64_t nc = 128;
  std::vector<double> phi_tgt1(nr, 0.0);
  std::vector<double> phi_tgt2(nr, 0.0);
  std::vector<double> phi_src(nc, 0.0);
  std::vector<double> mat(nr*nc, 0.0);
  for (int64_t ir = 0 ; ir < nr ; ++ir)
  for (int64_t ic = 0 ; ic < nc ; ++ic) {
    double v = distrib(rangen);
    mat[ir*nc+ic] = v;
  }

  DenseVectorCollector<int64_t, double> dvc(phi_src.data(), phi_tgt2.data(), nr, nc);
  
  for (int64_t ir = 0 ; ir < nr ; ++ir) {
    for (int64_t ic = 0 ; ic < nc ; ++ic) {
      phi_tgt1[ir] += mat[ir*nc+ic] * phi_src[ic];
      dvc(ir, ic, mat[ir*nc+ic]);
    }
  }
  for (int64_t ir = 0 ; ir < nr ; ++ir) {
    BOOST_CHECK_CLOSE(phi_tgt1[ir], phi_tgt2[ir], 1E-8);
  }
}

BOOST_AUTO_TEST_CASE(collector_test_int64)
{
  using namespace std;
  using namespace kore;
  using namespace kore::collector;
  mt19937_64 rangen;

  int64_t nr = 256;
  int64_t nc = 128;
  std::vector<int64_t> phi_tgt1(nr, 0.0);
  std::vector<int64_t> phi_tgt2(nr, 0.0);
  std::vector<int64_t> phi_src(nc, 0.0);
  std::vector<int64_t> mat(nr*nc, 0.0);
  for (int64_t ir = 0 ; ir < nr ; ++ir)
  for (int64_t ic = 0 ; ic < nc ; ++ic) {
    int64_t v = rangen();
    mat[ir*nc+ic] = v;
  }

  DenseVectorCollector<int64_t, int64_t> dvc(phi_src.data(), phi_tgt2.data(), nr, nc);
  
  for (int64_t ir = 0 ; ir < nr ; ++ir) {
    for (int64_t ic = 0 ; ic < nc ; ++ic) {
      phi_tgt1[ir] += mat[ir*nc+ic] * phi_src[ic];
      dvc(ir, ic, mat[ir*nc+ic]);
    }
  }
  for (int64_t ir = 0 ; ir < nr ; ++ir) {
    BOOST_CHECK_EQUAL(phi_tgt1[ir], phi_tgt2[ir]);
  }
}
