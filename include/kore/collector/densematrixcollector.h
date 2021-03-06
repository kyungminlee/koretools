#pragma once
#include <iterator>
#include "../utility/typecheck.h"

namespace kore {
namespace collector {

template <typename IndexType, typename ScalarType>
struct DenseMatrixCollector
{
  ScalarType * const target;
  IndexType const n_row;
  IndexType const n_col;

  static_assert(std::numeric_limits<IndexType>::is_signed, "IndexType should be signed");
  static_assert(std::numeric_limits<IndexType>::is_integer, "IndexType should be integer");

  static const bool RowRequired = true;
  static const bool ColRequired = true;
  static const bool ValRequired = true;

  DenseMatrixCollector(ScalarType *the_target, IndexType nr, IndexType nc)
      : target(the_target), n_row(nr), n_col(nc)
  {
    assert(0 < nr && 0 < nc);
  }

  void operator()(IndexType row, IndexType col, ScalarType val)
  {
    assert(0 <= row && row < n_row);
    assert(0 <= col && col < n_col);
    target[row * n_col + col] += val;
  }
};

} // namespace collector
} // namespace kore

