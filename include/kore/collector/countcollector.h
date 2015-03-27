#pragma once
#include <cassert>
#include "../utility/typecheck.h"

namespace kore {
namespace collector {

template <typename IndexType, typename ScalarType>
struct CountCollector
{
 public:
  IndexType const n_row;
  IndexType const n_col;
  IndexType count;

  // Assertion
  kore::static_assertion<std::numeric_limits<IndexType>::is_signed> check_signed;
  kore::static_assertion<std::numeric_limits<IndexType>::is_integer> check_integer;

  static const bool RowRequired = false;
  static const bool ColRequired = false;
  static const bool ValRequired = false;

  CountCollector(IndexType nr, IndexType nc) : n_row(nr), n_col(nc), count(0) {
  }

  void operator()(IndexType row, IndexType col, ScalarType val) { 
      assert(0 <= row && row < n_row);
      assert(0 <= col && col < n_col);
      count++; 
  } 
};


} // namespace collector
} // namespace kore

