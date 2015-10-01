#pragma once

#include <cassert>
#include <cstddef>
#include <iterator>
#include <array>
#include <vector>
#include <stdexcept>
#include <memory>
#include <type_traits>
#include <vector>
#include <initializer_list>
#include <iostream>

#include "../typedefs.h"

namespace kore {
namespace array {

template <std::size_t Rank> class ArrayStructure;
template <typename ValueType, std::size_t Rank> class Array;
template <typename ValueType, std::size_t Rank> class ConstArray;

template <std::size_t Rank>
std::ostream & operator<<(std::ostream & os, ArrayStructure<Rank> const & as);

template <typename ValueType, std::size_t Rank>
std::ostream & operator<<(std::ostream & os, Array<ValueType, Rank> const & arr);

template <typename ValueType, std::size_t Rank>
std::ostream & operator<<(std::ostream & os, ConstArray<ValueType, Rank> const & arr);

//! Structure of array
template <std::size_t _Rank>
class ArrayStructure
{
public:
  typedef std::size_t SizeType;
  typedef std::ptrdiff_t DifferenceType;

  static const std::size_t Rank = _Rank;

  static_assert(Rank > 0, "Rank should be greater than 1");

  ArrayStructure(ArrayStructure&& as)
    : length_(as.length_)
    , shape_(std::move(as.shape_))
    , stride_(std::move(as.stride_))
  {
  }

  ArrayStructure(const ArrayStructure& as)
    : length_(as.length_)
    , shape_(as.shape_)
    , stride_(as.stride_)
  {
  }

  //! Construct array structure with given shape.
  template <typename ... Args>
  ArrayStructure(SizeType size, Args ... args) { init_(size, args...); }

  //! Ravel multiindex to index
  //!
  template <typename ...Args>
  SizeType index(SizeType s, Args ... args) const {
    static_assert(Rank == 1 + sizeof...(Args), "Number of arguments should match Rank");
    return ravel_(s, args...);
  }

  //! Ravel multiindex to index
  //!
  template <typename ...Args>
  SizeType ravel(SizeType s, Args ... args) const {
    static_assert(Rank == 1 + sizeof...(Args), "Number of arguments should match Rank");
    return ravel_(s, args...);
  }

  //! Ravel multiindex to index
  //!
  SizeType ravel(const std::array<SizeType, Rank>& mi) const
  {
    SizeType idx = 0;
    for (SizeType d = 0u; d < Rank && d >= 0; ++d) {
      assert(0 <= mi[d] && mi[d] < shape_[d]);
      idx += stride_[d] * mi[d];
    }
    return idx;
  }

  //! Unravel index to multiindex
  //!
  //! 
  std::array<SizeType, Rank> unravel(SizeType s) const
  {
    std::array<SizeType, Rank> mi;
    assert(s < size());
    for (SizeType d = Rank - 1; d < Rank && d >= 0; --d) {
      mi[d] = s % shape_[d];
      s /= shape_[d];
    }
    return mi;
  }

  void show() const
  {
    std::cout << "length_ : " << length_ << std::endl;
    std::cout << "shape_ :";
    for (int d = 0; d < Rank; ++d) {
      std::cout << " " << shape_[d];
    }
    std::cout << std::endl;

    std::cout << "stride_ :";
    for (int d = 0; d < Rank; ++d) {
      std::cout << " " << stride_[d];
    }
    std::cout << std::endl;
  }

  SizeType length() const { return length_; }
  SizeType size() const { return length_; }
  std::array<SizeType, Rank> shape() const { return shape_; }
  std::array<SizeType, Rank> stride() const { return stride_; }

  SizeType shape(size_t d) const
  {
    assert(d < Rank);
    return shape_[d];
  }

  SizeType stride(size_t d) const
  {
    assert(d < Rank);
    return stride_[d];
  }

  friend std::ostream& operator<< <Rank>(std::ostream & os, ArrayStructure<Rank> const &);

private:
  //! Initialize array structure.
  //! \param size
  //! \param args
  template <typename ... Args>
  void init_(SizeType size, Args ... args)
  {
    init_(args...);
    static const std::size_t D = Rank - 1 - sizeof...(Args);
    static_assert(D < Rank, "D should be smaller than Rank");
    shape_[D] = size;
    stride_[D] = shape_[D + 1] * stride_[D + 1];
    length_ *= shape_[D];
  }

  //! Initialize array structure.
  //! \param size
  void init_(SizeType size)
  {
    static const std::size_t D = Rank - 1;
    static_assert(D + 1 == Rank, "last D+1 should match Rank");
    shape_[D] = size;
    stride_[D] = 1;
    length_ = size;
  }

  //! internal function for raveling multiindex to index
  //!
  template<typename ...Args>
  SizeType ravel_(SizeType s, Args... args) const
  {
    static const size_t D = Rank - 1 - sizeof...(Args);
    static_assert(D < Rank, "D should be smaller than Rank");
    assert(s < shape_[D]);
    return s * stride_[D] + ravel_(args...);
  }

  //! internal function for raveling multiindex to index
  //!
  SizeType ravel_(SizeType s) const
  {
    static const size_t D = Rank - 1;
    static_assert(D + 1 == Rank, "last D+1 should match Rank");
    assert(s < shape_[D]);
    return s * stride_[D];
  }

private:
  SizeType length_; //!< Total number of elements.
  std::array<SizeType, Rank> shape_; //!< Shape
  std::array<SizeType, Rank> stride_; //!< Stride
}; // class ArrayStructure


//! Array<ValueType, Rank>
template <typename _ValueType, std::size_t _Rank>
class Array
{
public:
  static const std::size_t Rank = _Rank;
  using ValueType = _ValueType;
  using ReferenceType = ValueType &;
  using ConstReferenceType = ValueType const &;
  using PointerType = ValueType *;
  using ConstPointerType = ValueType const *;
  using ArrayStructureType = ArrayStructure<Rank>;
  using SizeType = typename ArrayStructure<Rank>::SizeType;

  // Member Functions
  template <typename ... Args>
  explicit Array(SizeType s, Args ... args)
    : structure_(s, args...)
    , data_(new ValueType[structure_.size()], std::default_delete<ValueType[]>())
    , own_(true)
  {
  }

  explicit Array(const ArrayStructureType& structure)
    : structure_(structure)
    , data_(new ValueType[structure_.size()], std::default_delete<ValueType[]>())
    , own_(true)
  {
  }

  Array(const Array& arr)
    : structure_(arr.structure_)
    , data_(arr.data_)
    , own_(arr.own_)
  {
  }

  Array(Array&& arr)
    : structure_(std::move(arr.structure_))
    , data_(std::move(arr.data_))
    , own_(arr.own_)
  {
  }

  Array(ConstArray<ValueType, Rank> const & arr)
    : Array(arr.structure())
  {
    std::copy(arr.begin(), arr.end(), begin());
  }

  template<typename ... Args>
  Array(std::shared_ptr<ValueType>& data, bool own, Args ... args)
    : structure_(args...)
    , data_(data)
    , own_(own)
  {
  }

  /// Construct from raw pointer
  /// Use null deleter
  template <typename ... Args>
  Array(PointerType data, Args ... args)
    : structure_(args...)
    , data_(data, [](void const *) {})
    , own_(false)
  {
  }

  template <typename ... Args>
  Array<ValueType, sizeof...(Args)>
    reshape(Args ... args)
  {
    Array<ValueType, sizeof...(Args)> ret(data_, own_, args...);
#ifndef NDEBUG
    if (ret.size() != size()) throw std::length_error("Array::reshape()");
#endif
    return ret;
  }

  ConstArray<ValueType, Rank>
    constant() const
  {
    ConstArray<ValueType, Rank> ret(*this);
    return ret;
  }

  template <typename ... Args>
  ReferenceType operator()(Args ... args)
  {
    static_assert(Rank == sizeof...(Args), "Number of arguments should match Rank");
    SizeType idx = structure_.index(args...);
    return data_.get()[idx];
  }

  template <typename ... Args>
  ValueType operator()(Args ... args) const
  {
    static_assert(Rank == sizeof...(Args), "Number of arguments should match Rank");
    SizeType idx = structure_.index(args...);
    return data_.get()[idx];
  }

  template <typename ... Args>
  ReferenceType at(Args ... args)
  {
    static_assert(Rank == sizeof...(Args), "Number of arguments should match Rank");
    SizeType idx = structure_.index(args...);
    return data_.get()[idx];
  }

  template <typename ... Args>
  ValueType at(Args ... args) const
  {
    static_assert(Rank == sizeof...(Args), "Number of arguments should match Rank");
    SizeType idx = structure_.index(args...);
    return data_.get()[idx];
  }

  void fill(ConstReferenceType v) {
    std::fill(begin(), end(), v);
  }

  Array<ValueType, Rank> clone() const {
    Array<ValueType, Rank> ret(structure_);
    std::copy(cbegin(), cend(), ret.begin());
    return ret;
  }

  const ArrayStructureType& structure() const { return structure_; }

  std::shared_ptr<const ValueType> data() const { return data_; }
  bool own() const { return own_; }
  PointerType begin() { return data_.get(); }
  PointerType end() { return data_.get() + length(); }
  ConstPointerType begin() const { return data_.get(); }
  ConstPointerType end() const { return data_.get() + length(); }
  ConstPointerType cbegin() const { return data_.get(); }
  ConstPointerType cend() const { return data_.get() + length(); }

  SizeType length() const { return structure_.length(); }
  SizeType size() const { return structure_.size(); }
  std::array<SizeType, Rank> shape() const { return structure_.shape(); }
  std::array<SizeType, Rank> stride() const { return structure_.stride(); }
  SizeType shape(size_t d) const { return structure_.shape(d); }
  SizeType stride(size_t d) const { return structure_.stride(d); }

  //friend class ConstArray<ValueType, Rank>;
  friend std::ostream& operator<< <ValueType, Rank>(std::ostream & os, Array<ValueType, Rank> const &);
private:
  ArrayStructureType structure_;
  std::shared_ptr<ValueType> data_;
  bool own_;
}; // class Array



template <typename _ValueType, std::size_t _Rank>
class ConstArray
{
public:
  static const std::size_t Rank = _Rank;
  using ValueType = _ValueType;
  using ReferenceType = ValueType &;
  using ConstReferenceType = ValueType const &;
  using PointerType = ValueType *;
  using ConstPointerType = ValueType const *;
  using ArrayStructureType = ArrayStructure<Rank>;
  using SizeType = typename ArrayStructure<Rank>::SizeType;

  ConstArray(const Array<ValueType, Rank>& arr)
    : structure_(arr.structure())
    , data_(arr.data())
    , own_(arr.own())
  {
  }

  ConstArray(ConstArray&& arr)
    : structure_(std::move(arr.structure_))
    , data_(std::move(arr.data_))
    , own_(std::move(arr.own_))
  {
  }

  ConstArray(const ConstArray& arr)
    : structure_(arr.structure_)
    , data_(arr.data_)
    , own_(arr.own_)
  {
  }

  template<typename ... Args>
  ConstArray(std::shared_ptr<const ValueType>& data, bool own, Args ... args)
    : structure_(args...)
    , data_(data)
    , own_(own)
  {
  }

  /// Construct from raw pointer
  /// Use null deleter
  template <typename ... Args>
  ConstArray(PointerType data, Args ... args)
    : structure_(args...)
    , data_(const_cast<PointerType>(data), [](void const *) {})
    , own_(false)
  {
  }

  template <typename ... Args>
  ConstArray(ConstPointerType data, Args ... args)
    : structure_(args...)
    , data_(data, [](void const *) {})
    , own_(false)
  {
  }

  template <typename ... Args>
  ConstArray<ValueType, sizeof...(Args)>
    reshape(Args ... args)
  {
    ConstArray<ValueType, sizeof...(Args)> ret(data_, own_, args...);
#ifndef NDEBUG
    if (ret.size() != size()) throw std::length_error("ConstArray::reshape()");
#endif
    return ret;
  }

  template <typename ... Args>
  ValueType operator()(Args ... args) const
  {
    static_assert(Rank == sizeof...(Args), "Number of arguments should match Rank");
    SizeType idx = structure_.index(args...);
    return data_.get()[idx];
  }

  template <typename ... Args>
  ValueType at(Args ... args) const
  {
    static_assert(Rank == sizeof...(Args), "Number of arguments should match Rank");
    SizeType idx = structure_.index(args...);
    return data_.get()[idx];
  }

  Array<ValueType, Rank> clone() const {
    Array<ValueType, Rank> ret(structure_);
    std::copy(cbegin(), cend(), ret.begin());
    return ret;
  }

  const ArrayStructureType& structure() const { return structure_; }

  bool own() const { return own_; }
  //ConstPointerType data() const { return data_.get(); }
  std::shared_ptr<const ValueType> data() const { return data_; }
  ConstPointerType begin() const { return data_.get(); }
  ConstPointerType end() const { return data_.get() + length(); }
  ConstPointerType cbegin() const { return data_.get(); }
  ConstPointerType cend() const { return data_.get() + length(); }

  SizeType length() const { return structure_.length(); }
  SizeType size() const { return structure_.size(); }
  std::array<SizeType, Rank> shape() const { return structure_.shape(); }
  std::array<SizeType, Rank> stride() const { return structure_.stride(); }
  SizeType shape(size_t d) const { return structure_.shape(d); }
  SizeType stride(size_t d) const { return structure_.stride(d); }

  friend std::ostream& operator<< <ValueType, Rank>(std::ostream & os, ConstArray<ValueType, Rank> const &);
private:
  ArrayStructureType structure_;
  std::shared_ptr<const ValueType> data_;
  bool own_;
}; // class ConstArray

} // namespace array
} // namespace kore

#include "array_impl.h"
