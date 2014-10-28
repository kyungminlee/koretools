#pragma once

namespace kore {
namespace bitbox {
  
  template <size_t Dim, typename BitString, typename Index>
    template <typename ... Args>
    BitShifter<Dim, BitString, Index>::BitShifter(Args ... args)
    {
      set_shape(args...);
      size_ = shape_[0] * stride_[0];
      mask_all_ = makemask<BitString>(0, size_-1);

      const BitString ONE = 0x1;
      for (Index d = 0; d < Dim; ++d) {
	mask_[d] = 0;
        for (Index i = 0; i < shape_[d]; ++i) {
          mask_[d] |= ONE << (i * stride_[d]);
        }
      }

    }
  
  template <size_t Dim, typename BitString, typename Index>
    template <typename ...Args> inline
    void BitShifter<Dim, BitString, Index>::set_shape(Index n, Args... args)
    {
      set_shape(args...);
      shape_[Dim - sizeof...(args)-1] = n;
      stride_[Dim - sizeof...(args)-1] = stride_[Dim-sizeof...(args)] * shape_[Dim-sizeof...(args)];
    }
  
  template <size_t Dim, typename BitString, typename Index>
    inline
    void BitShifter<Dim, BitString, Index>::set_shape(Index n)
    {
      shape_[Dim -1] = n;
      stride_[Dim-1] = 1;
    }
  
}
}
