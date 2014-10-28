#pragma once
#include <array>
#include "../typedefs.h"

//
// Collection of bit operations.
//

namespace kore {
namespace bitbox {

  template <size_t _Dim, typename _BitString, typename _Index>
    class BitShifter
  {
  public:
    static const size_t Dim = _Dim;
    typedef _BitString BitString;
    typedef _Index Index;
    
    //template <typename ... Args, size_t D>
    template <typename ... Args>
    BitShifter(Args ... args);
    
    template <typename ...Args>
    void set_shape(Index n, Args... args);
    void set_shape(Index n);
    
    template <typename ... Args>
    BitString operator()(BitString v, Index s, Args ... args) const;
    BitString operator()(BitString v, Index s) const;
    
    Index size() const { return size_; }
    BitString mask_all() const { return mask_all_; }
    const std::array<Index, Dim>& shape() const { return shape_; }
    const std::array<Index, Dim>& stride() const { return stride_; }
    const std::array<BitString, Dim>& mask() const { return mask_; }
  private:
    Index size_;
    BitString mask_all_;
    std::array<Index, Dim> shape_;
    std::array<Index, Dim> stride_;
    std::array<BitString, Dim> mask_;
  };
  
  template <typename _BitString, typename _Index> 
  class BitShifter2D
  {
  public:
    typedef _BitString BitString;
    typedef _Index Index;

    BitShifter2D(Index n1, Index n2) {
      assert(n1 > 0 && n2 > 0);

      size_ = n1*n2;
      shape_[0]  = n1;     shape_[1]  = n2;
      stride_[0] = n2;     stride_[1] = 1;

      mask_all_ = makemask<BitString>(0, size_ - 1);
      
      mask_[0] = 0; mask_[1] = 0;
      const BitString ONE = 0x1;
      for (Index d = 0; d < 2; ++d) {
        for (Index i = 0; i < shape_[d]; ++i) {
          mask_[d] |= ONE << (i * stride_[d]);
        }
      }
    }

    BitString operator()(BitString v, Index i1, Index i2) const {
      for (Index i = 0; i < shape_[1]; ++i) {
        BitString mask = mask_[0] << i;
        v = bitrotate(v, i1, mask);
      }

      for (Index i = 0; i < shape_[0]; ++i) {
        BitString mask = mask_[1] << (i * stride_[0]);
        v = bitrotate(v, i2, mask);
      }
      return v;
    }

    const std::array<Index, 2>& shape() const { return shape_; }
    const std::array<Index, 2>& stride() const { return stride_; }
    Index size() const { return size_; }

  private:
    Index size_;
    BitString mask_all_;
    std::array<Index, 2> shape_;
    std::array<Index, 2> stride_;
    std::array<BitString, 2> mask_;
  };


} // namespace bitbox
} // namespace kore
#include "bitshifter_impl.h"

