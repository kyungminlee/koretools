#include <cassert>
//#include <random>
#include <map>
#include <iostream>
#include <bitset>
#include <array>

#include <kore/kore>
//#include <kore/dictionary.h>

#include "tuple_stream.h"

typedef kore::uint64_t BitString;
typedef kore::int64_t Integer;
using namespace kore::bitbox;
const size_t BITSIZE = 32;


template <size_t Dim>
class BitShifter
{
public:
  BitShifter(std::initializer_list<size_t> shape) {
    if (shape.size() != Dim) {
      throw std::invalid_argument("BitShifter: number of arguments does not match");
    }
    std::copy(shape.begin(), shape.end(), shape_.begin());
    size_t stride = 1;
    BitString ONE = 0x1;
    for (size_t d = Dim - 1; d < Dim && d >= 0; d--) {
      stride_[d] = stride;

      BitString mask = 0x0;

      for (size_t i = 0; i < shape_[d]; ++i) {
        mask |= (ONE << (i * stride));
      }
      mask_[d] = mask;

      stride *= shape_[d];
    }
    size_ = stride;
    //set_stride<Dim - 1>(1);
    mask_all_ = makemask<BitString>(0, size_ - 1);
  }

public:
  size_t size_;
  std::array<size_t, Dim> shape_;
  std::array<size_t, Dim> stride_;
  std::array<BitString, Dim> mask_;
  BitString mask_all_;
};

class BitShifter1D
{
public:
  BitShifter1D(size_t n) : n_(n), mask_(makemask<BitString>(0, n - 1)) { 
    if (n == 0) { throw std::invalid_argument("BitShifter1D: size should be positive"); }
  }

  template <typename SizeInt>
  BitString operator()(BitString v, SizeInt i) const {
    return kore::bitbox::bitrotate(v, i, mask_);
  }
  size_t size() const { return n_; }
  BitString mask() const { return mask_; }
private:
  size_t n_;
  BitString mask_;
};


class BitShifter2D
{
public:
  BitShifter2D(size_t n1, size_t n2) : mask_all_(0){
    if (n1 == 0 || n2 == 0) {
      throw std::invalid_argument("BitShifter2D: size should be positive");
    }
    static const BitString ONE = 0x1;
    size_ = n1*n2;
    shape_[0] = n1; shape_[1] = n2;
    stride_[0] = n2; stride_[1] = 1;
    mask_all_ = makemask<BitString>(0, n1*n2 - 1);
    mask_[0] = 0; mask_[1] = 0;
    for (size_t i = 0; i < n2; ++i) { mask_[1] |= ONE << i; }
    for (size_t i = 0; i < n1; ++i) { mask_[0] |= ONE << (i*n2); }
  }

  template <typename SizeInt>
  BitString operator()(BitString v, SizeInt i1, SizeInt i2) const {
    BitString mask = mask_[1];
    for (size_t i = 0; i < shape_[0]; ++i) {
      v = bitrotate(v, i2, mask);
      mask = mask << stride_[0];
    }
    mask = mask_[0];
    for (size_t i = 0; i < shape_[1]; ++i) {
      v = bitrotate(v, i1, mask);
      mask = mask << 1;
    }
    return v;
  }
public:
  size_t size_;
  std::array<size_t, 2> stride_;
  std::array<size_t, 2> shape_;
  std::array<BitString, 2> mask_;
  BitString mask_all_;
};


class BitShifter3D
{
public:
  BitShifter3D(size_t n1, size_t n2, size_t n3) : mask_all_(0x0) {
    if (n1 == 0 || n2 == 0 || n3 == 0) {
      throw std::invalid_argument("BitShifter3D: size should be positive");
    }
    static const BitString ONE = 0x1;
    size_ = n1*n2*n3;
    shape_[0] = n1; shape_[1] = n2; shape_[2] = n3;
    stride_[0] = n2*n3; stride_[1] = n3; stride_[2] = 1;
    mask_all_ = makemask<BitString>(0, size_ - 1);
    mask_[0] = 0; mask_[1] = 0; mask_[2] = 0;
    for (size_t d = 0; d < 3; ++d) {
      for (size_t i = 0; i < shape_[d]; ++i){
        mask_[d] |= ONE << (i*stride_[d]);
      }
    }
  }

  template <typename SizeInt = ptrdiff_t>
  BitString operator()(BitString v, SizeInt i1, SizeInt i2, SizeInt i3) const {
    std::array<size_t, 3> n_loop;
    std::array<SizeInt, 3> i_loop;
    i_loop[0] = i1;
    i_loop[1] = i2;
    i_loop[2] = i3;

    for (int d = 0; d < 3; ++d) {
      n_loop = shape_;
      n_loop[d] = 1;
      for (size_t j1 = 0; j1 < n_loop[0]; ++j1) {
        for (size_t j2 = 0; j2 < n_loop[1]; ++j2) {
          for (size_t j3 = 0; j3 < n_loop[2]; ++j3) {
            BitString mask = mask_[d] << (j1 * stride_[0] + j2 * stride_[1] + j3*stride_[2]);
            v = bitrotate(v, i_loop[d], mask);
          }
        }
      }
    }

    return v;
  }
  size_t size_;
  std::array<size_t, 3> stride_;
  std::array<size_t, 3> shape_;
  std::array<BitString, 3> mask_;
  BitString mask_all_;
};


class FermionTranslationDictionary1D
{
public:
  FermionTranslationDictionary1D(size_t n) : length_(n), shifter_(n) {
    static const BitString ONE = 0x1;
    for (BitString word = 0x0; !(word & (ONE << n)); ++word) {
      auto parsed_word = parse(word);
      auto entry = std::get<0>(parsed_word);
      auto it = indices_.find(entry);
      if (it == indices_.end()) { // not
        indices_[entry] = words_.size();
        words_.push_back(entry);
        lengths_.push_back(std::get<2>(parsed_word));
      }
    }
  }

  void show()
  {
    std::cout << "=== WORD ===" << std::endl;
    for (auto w : words_) {
      std::cout << std::bitset<BITSIZE>(w) << std::endl;
    }
    std::cout << "=== LENGTH ===" << std::endl;
    for (auto w : lengths_) {
      std::cout << w << std::endl;
    }
    std::cout << std::endl;
    std::cout << "=== INDEX ===" << std::endl;
    for (auto w : indices_) {
      std::cout << std::bitset<BITSIZE>(w.first) << " : " << w.second << std::endl;
    }

    std::cout << std::endl;

  }
  
  std::tuple<BitString, Integer, Integer> parse(BitString v) const {
    Integer min_shift = 0;
    Integer shift;
    BitString min_v = v;
    // TODO: more efficient to look at divisors.
    for (shift = 1; shift < length_; ++shift) {
      BitString v2 = shifter_(v, shift);
      if (v2 == v) {
        break;
      }
      if (v2 < min_v) { min_shift = shift; min_v = v2; }
    }
    BitString v2 = shifter_(v, min_shift);
    return std::make_tuple(v2, min_shift, shift);
  }

  BitString word(Integer idx) const { return words_[idx]; }
  Integer index(BitString word) const { return indices_.at(word); }
  size_t length() const { return length_; }
  size_t size() const { return words_.size(); }
private:
  size_t length_;
  BitShifter1D shifter_;
  std::map<BitString, Integer> indices_;
  std::vector<BitString> words_;
  std::vector<Integer> lengths_;
};


int main(int argc, char** argv)
{

  FermionTranslationDictionary1D dict(5);
  //for (BitString i = 0; i < 16; ++i) {
  //  std::cout << dict.parse(i) << std::endl;
  //}
  using namespace std;
  
  dict.show();
  cout << dict.size() << endl;
  cout << dict.length() << endl;
  for (Integer i = 0 ; i < dict.size() ; ++i) {
    BitString w = dict.word(i);
    Integer i2 = dict.index(w);
    std::cout << i << " => " << std::bitset<BITSIZE>(w) << " => " << i2 << std::endl;
  }
  
  return 0;
}

int main2(int argc, char** argv)
{
  std::cout << std::bitset<BITSIZE>(bitrotate<BitString>(0x1, -1, 0x11111)) << std::endl;


  for (int i = -7; i < 10; ++i) {
    std::cout << i << " : " << bitrotate<uint64_t>(5, i, 6) << std::endl;
  }
  //size_t nx = 4, ny = 3, nz = 2;
  {
    BitShifter1D shifter(10);
    std::cout << shifter.size() << std::endl;
    std::cout << std::bitset<BITSIZE>(shifter.mask()) << std::endl;
    for (int i = -7; i < 10; ++i) {
      std::cout << std::bitset<BITSIZE>(shifter(5, i)) << " : " << i << std::endl;
    }
  }

  {
    BitShifter2D shifter(5, 4);
    std::cout << shifter.shape_ << std::endl;
    std::cout <<
      std::bitset<BITSIZE>(shifter.mask_[0]) << ", " <<
      std::bitset<BITSIZE>(shifter.mask_[1]) << std::endl;
    std::cout << std::bitset<BITSIZE>(shifter.mask_all_) << std::endl;
    for (int i = -7; i < 10; ++i) {
      std::cout << std::bitset<BITSIZE>(shifter(0x1, i, 0)) << " : " << i << std::endl;
    }
  }
  {
    BitShifter3D shifter(4,3,2);
    std::cout << shifter.shape_ << std::endl;
    std::cout << std::bitset<BITSIZE>(shifter.mask_[0]) << " : x" << std::endl;
    std::cout << std::bitset<BITSIZE>(shifter.mask_[1]) << " : y" << std::endl;
    std::cout << std::bitset<BITSIZE>(shifter.mask_[2]) << " : y" << std::endl;
    std::cout << std::bitset<BITSIZE>(shifter.mask_all_) << " : all" << std::endl;
    for (int i = -7; i < 10; ++i) {
      std::cout << std::bitset<BITSIZE>(shifter(0x6, 0, i, 0)) << " : " << i << std::endl;
    }
  }

  return 0;
}
