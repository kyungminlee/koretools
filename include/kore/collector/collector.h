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
struct CountCollector
{
 public:
  IndexType const n_row;
  IndexType const n_col;
  IndexType count;

  CountCollector(IndexType nr, IndexType nc) : n_row(nr), n_col(nc), count(0) { }

  void operator()(IndexType row, IndexType col, ScalarType val) { 
      assert(0 <= row && row < n_row);
      assert(0 <= col && col < n_val);
      count++; 
  } 
};

template <typename IndexType, typename ScalarType, typename RowIterator, typename ColIterator, typename ValIterator>
struct SparseCooMatrixCollector
{
 public:
  IndexType const n_row;
  IndexType const n_col;
  RowIterator row_iterator;
  ColIterator col_iterator;
  ValIterator val_iterator;

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

