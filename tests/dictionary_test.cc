#include <cassert>
//#include <random>
#include <map>
#include <iostream>
#include <bitset>
#include <array>

#include <kore/kore>
//#include <kore/dictionary.h>

namespace aux
{
  template<std::size_t...> struct seq {};

  template<std::size_t N, std::size_t... Is>
  struct gen_seq : gen_seq<N - 1, N - 1, Is...> {};

  template<std::size_t... Is>
  struct gen_seq<0, Is...> : seq<Is...>{};

  template<class Ch, class Tr, class Tuple, std::size_t... Is>
  void print_tuple(std::basic_ostream<Ch, Tr>& os, Tuple const& t, seq<Is...>){
    using swallow = int[];
    (void) swallow{0, (void(os << (Is == 0 ? "" : ", ") << std::get<Is>(t)), 0)...};
  }
} // aux::

template<class Ch, class Tr, class... Args>
auto operator<<(std::basic_ostream<Ch, Tr>& os, std::tuple<Args...> const& t)
-> std::basic_ostream<Ch, Tr>&
{
  os << "(";
  aux::print_tuple(os, t, aux::gen_seq<sizeof...(Args)>());
  return os << ")";
}

template <typename T, size_t N>
std::ostream& operator<<(std::ostream& os, const std::array<T, N>& arr)
{
  os << "[";
  for (size_t i = 0; i < N; ++i) {
    os << arr[i] << ", ";
  }
  os << "]";
  return os;
}




typedef kore::uint64_t BitString;
typedef kore::int64_t Integer;
using namespace kore::bitbox;

#if 0
struct FermionDictionary
{
  FermionDictionary(Integer nx, Integer ny)
  : nx_(nx), ny_(ny)
  {
    Integer n = nx * ny;
    BitString ONE = 0x1;
    BitString val = 0;
    while (!((ONE << n) & val)) {

    }
  }

  std::tuple<BitString, Integer, Integer> reduce(BitString v)
  {
  }
  Integer nx_, ny_;
  //std::map < std::tuple<Integer, Integer>, std::map< std::tuple<BitString, Integer, Integer>, Integer> basis_to_index;
};
#endif

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

#if 0
  template <size_t Idx>
  typename std::enable_if< (Idx < sizeof...(Args)), BitString>::type
    shift_dim(BitString value, Args ... shift) const {
      BitString mask = mask_[Idx];
      for (size_t i = 0; i < std::get<Idx>(shape); ++i) {
        mask << (std::get<Idx>(stride_))// LOOP OVER THE REST!
      }
    }
#endif

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
  FermionTranslationDictionary1D(size_t n) : size_(n), shifter_(n) {
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
      std::cout << w << std::endl;
    }
    std::cout << "=== LENGTH ===" << std::endl;
    for (auto w : lengths_) {
      std::cout << w << std::endl;
    }
    std::cout << std::endl;
    std::cout << "=== INDEX ===" << std::endl;
    for (auto w : indices_) {
      std::cout << w.first << " : " << w.second << std::endl;
    }

    std::cout << std::endl;

  }
  std::tuple<BitString, Integer, Integer> parse(BitString v) const {
    Integer min_shift = 0;
    Integer shift;
    BitString min_v = v;
    // TODO: more efficient to look at divisors.
    for (shift = 1; shift < size_; ++shift) {
      BitString v2 = shifter_(v, shift);
      if (v2 == v) {
        break;
      }
      if (v2 < min_v) { min_shift = shift; min_v = v2; }
    }
    BitString v2 = shifter_(v, min_shift);
    return std::make_tuple(v2, min_shift, shift);
  }
private:
  size_t size_;
  BitShifter1D shifter_;
  std::map<BitString, Integer> indices_;
  std::vector<BitString> words_;
  std::vector<Integer> lengths_;
};


int main(int argc, char** argv)
{

  FermionTranslationDictionary1D dict(4);
  //for (BitString i = 0; i < 16; ++i) {
  //  std::cout << dict.parse(i) << std::endl;
  //}
  dict.show();
  return 0;
}

int main2(int argc, char** argv)
{
  const size_t BITSIZE = 32;
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


  //shift(1, 4, 2, 0, 0);
  return 0;
}

#if 0
BOOST_AUTO_TEST_CASE(collector_test_double)
{
  using namespace std;
  using namespace kore;
  using namespace kore::collector;
  mt19937_64 rangen;
  uniform_real_distribution<double> distrib(0.0, 1.0);

  int64_t nr = 256;
  int64_t nc = 128;
  std::vector<double> phi_tgt1(nr, 0.0);
  std::vector<double> phi_tgt2(nr, 0.0);
  std::vector<double> phi_src(nc, 0.0);
  std::vector<double> mat(nr*nc, 0.0);
  for (int64_t ir = 0 ; ir < nr ; ++ir)
  for (int64_t ic = 0 ; ic < nc ; ++ic) {
    double v = distrib(rangen);
    mat[ir*nc+ic] = v;
  }

  DenseVectorCollector<int64_t, double> dvc(phi_src.data(), phi_tgt2.data(), nr, nc);
  
  for (int64_t ir = 0 ; ir < nr ; ++ir) {
    for (int64_t ic = 0 ; ic < nc ; ++ic) {
      phi_tgt1[ir] += mat[ir*nc+ic] * phi_src[ic];
      dvc(ir, ic, mat[ir*nc+ic]);
    }
  }
  for (int64_t ir = 0 ; ir < nr ; ++ir) {
    BOOST_CHECK_CLOSE(phi_tgt1[ir], phi_tgt2[ir], 1E-8);
  }
}

BOOST_AUTO_TEST_CASE(collector_test_int64)
{
  using namespace std;
  using namespace kore;
  using namespace kore::collector;
  mt19937_64 rangen;

  int64_t nr = 256;
  int64_t nc = 128;
  std::vector<int64_t> phi_tgt1(nr, 0.0);
  std::vector<int64_t> phi_tgt2(nr, 0.0);
  std::vector<int64_t> phi_src(nc, 0.0);
  std::vector<int64_t> mat(nr*nc, 0.0);
  for (int64_t ir = 0 ; ir < nr ; ++ir)
  for (int64_t ic = 0 ; ic < nc ; ++ic) {
    int64_t v = rangen();
    mat[ir*nc+ic] = v;
  }

  DenseVectorCollector<int64_t, int64_t> dvc(phi_src.data(), phi_tgt2.data(), nr, nc);
  
  for (int64_t ir = 0 ; ir < nr ; ++ir) {
    for (int64_t ic = 0 ; ic < nc ; ++ic) {
      phi_tgt1[ir] += mat[ir*nc+ic] * phi_src[ic];
      dvc(ir, ic, mat[ir*nc+ic]);
    }
  }
  for (int64_t ir = 0 ; ir < nr ; ++ir) {
    BOOST_CHECK_EQUAL(phi_tgt1[ir], phi_tgt2[ir]);
  }
}
#endif
