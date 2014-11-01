#pragma once
#include <cassert>

namespace kore {
namespace collector {

template <typename IndexType, typename ScalarType>
struct DenseVectorCollector
{
  ScalarType const * const source;
  ScalarType       * const target;
  IndexType  const n_row;
  IndexType  const n_col;

#ifdef USE_CXX11
  static_assert( std::is_integral<IndexType>::value && std::is_signed<IndexType>::value, "IndexType should be signed integer");
#endif

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

template <typename IndexType, typename ScalarType>
struct DenseMatrixCollector
{
  ScalarType * const target;
  IndexType const n_row;
  IndexType const n_col;
  
#ifdef USE_CXX11
  static_assert( std::is_integral<IndexType>::value && std::is_signed<IndexType>::value, "IndexType should be signed integer");
#endif
  
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


template <typename IndexType, typename ScalarType>
struct SparseCooMatrixCollector
{
 public:
  IndexType const n_row;
  IndexType const n_col;
  std::vector<IndexType> rows, cols;
  std::vector<ScalarType> vals;

  SparseCooMatrixCollector(IndexType nr, IndexType nc) : n_row(nr), n_col(nc) { }
  
  void operator()(IndexType row, IndexType col, ScalarType val)
  {
    rows.push_back(row);
    cols.push_back(col);
    vals.push_back(val);
  }

  IndexType size() const { return vals.size(); }

  template <typename RowIterator, typename ColIterator, typename ValIterator>
  void fetch(RowIterator ri, ColIterator ci, ValIterator vi) const
  {
    std::copy(rows.begin(), rows.end(), ri);
    std::copy(cols.begin(), cols.end(), ci);
    std::copy(vals.begin(), vals.end(), vi);
  }

  void clear() { rows.clear(); cols.clear(); vals.clear(); }
};






} // namespace collector
} // namespace kore

