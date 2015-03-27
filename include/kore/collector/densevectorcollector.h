#pragma once
#include <cassert>
#include <iterator>
#include "../utility/typecheck.h"

namespace kore {
namespace collector {

template <typename IndexType, typename ScalarType>
struct DenseVectorCollector
{
  ScalarType const * const source;
  ScalarType       * const target;
  IndexType  const n_row;
  IndexType  const n_col;

  static_assert(std::numeric_limits<IndexType>::is_signed, "IndexType should be signed");
  static_assert(std::numeric_limits<IndexType>::is_integer, "IndexType should be integer");

  static const bool RowRequired = true;
  static const bool ColRequired = true;
  static const bool ValRequired = true;
  
  DenseVectorCollector(ScalarType const * const the_source,
                       ScalarType       * const the_target,
                       IndexType nr, IndexType nc)
      : source(the_source), target(the_target), n_row(nr), n_col(nc)
  {
    assert(0 < nr && 0 < nc);
  }
  
  void operator()(IndexType row, IndexType col, ScalarType val)
  {
    assert(0 <= row && row < n_row);
    assert(0 <= col && col < n_col);
    target[row] += val * source[col];
  }
};


} // namespace collector
} // namespace kore

