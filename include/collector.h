#pragma once
namespace kore
{

  template <typename IndexType, typename ScalarType>
  struct DenseVectorCollector
  {
    ScalarType const * const source;
    ScalarType       * const target;
    IndexType  const n_row;
    IndexType  const n_col;

    DenseVectorCollector(ScalarType const * const the_source,
                         ScalarType       * const the_target,
                         IndexType nr, IndexType nc)
                         : source(the_source), target(the_target), n_row(nr), n_col(nc)
    {
    }

    void operator()(IndexType row, IndexType col, ScalarType val)
    {
      assert(0 <= row && row < n_row);
      assert(0 <= col && col < n_col);
      target[row] += val * source[col];
    }

  };

}