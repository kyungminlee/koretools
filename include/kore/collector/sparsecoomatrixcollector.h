#pragma once
#include "../utility/typecheck.h"

namespace kore {
namespace collector {

template <typename IndexType, typename ScalarType,
          typename RowIterator, typename ColIterator, typename ValIterator>
struct SparseCooMatrixCollector
{
 public:
  IndexType const n_row;
  IndexType const n_col;
  RowIterator row_iterator;
  ColIterator col_iterator;
  ValIterator val_iterator;

  // Assertions
  static_assert(std::numeric_limits<IndexType>::is_signed, "IndexType should be signed");
  static_assert(std::numeric_limits<IndexType>::is_integer, "IndexType should be integer");
  static_assert((kore::is_same<typename std::iterator_traits<RowIterator>::value_type, IndexType>::value), "IndexType should equal RowIterator");
  static_assert((kore::is_same<typename std::iterator_traits<ColIterator>::value_type, IndexType>::value), "IndexType should equal ColIterator");
  static_assert((kore::is_same<typename std::iterator_traits<ValIterator>::value_type, ScalarType>::value), "IndexType should equal ValIterator");

  static const bool RowRequired = true;
  static const bool ColRequired = true;
  static const bool ValRequired = true;

  SparseCooMatrixCollector(IndexType nr, IndexType nc,
                           RowIterator ri, ColIterator ci, ValIterator vi)
      : n_row(nr), n_col(nc), row_iterator(ri), col_iterator(ci), val_iterator(vi) { }

  void operator()(IndexType row, IndexType col, ScalarType val)
  {
      assert(0 <= row && row < n_row);
      assert(0 <= col && col < n_col);
      *row_iterator++ = row;
      *col_iterator++ = col;
      *val_iterator++ = val;
  }
};

} // namespace collector
} // namespace kore

