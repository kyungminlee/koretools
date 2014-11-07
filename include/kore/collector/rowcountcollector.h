#pragma once
#include <cassert>
#include <vector>
#include "../utility/typecheck.h"

namespace kore {
namespace collector {

template <typename IndexType, typename ScalarType>
struct RowCountCollector
{
 public:
  IndexType const n_row;
  IndexType const n_col;
  IndexType row_start, row_end;
  std::vector<IndexType> row_count;
  
  kore::static_assertion<std::numeric_limits<IndexType>::is_signed> check_signed;
  kore::static_assertion<std::numeric_limits<IndexType>::is_integer> check_integer;

  static const bool RowRequired = true;
  static const bool ColRequired = false;
  static const bool ValRequired = false;
  
  RowCountCollector(IndexType nr, IndexType nc,
                    IndexType the_row_start=0, IndexType the_row_end=-1)
      : n_row(nr), n_col(nc), row_start(the_row_start), row_end(the_row_end)
  {
    assert(0 < n_row);
    assert(0 < n_col);
    if (row_end < 0) { row_end = n_row; }
    assert(row_start <= n_row);
    assert(row_start <= row_end && row_end <= n_row);
    
    row_count.resize(row_start - row_end, static_cast<IndexType>(0));
  }
  
  void operator()(IndexType row, IndexType col, ScalarType val) {
    assert(row_start <= row && row < row_end);
    row_count[row - row_start]++;
  }

}; // struct RowCountCollector

} // collector
} // kore
